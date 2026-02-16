from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Dict, Literal, Optional

import h5py
import numpy as np
import scipy as sp
from numpy.typing import NDArray
from scipy.sparse import csr_matrix


log = logging.getLogger(__name__)


# ------------------------------------------------------------
# Binary readers (.rhs/.sol vectors, .mat CSR matrix)
# ------------------------------------------------------------
def read_vec(filename: Path) -> tuple[int, NDArray[np.float64]]:
    """
    Read a vector file of the form:
        int64 nv
        float64[nv] values

    Returns (nv, vector), validating length.
    """
    filename = Path(filename)
    with filename.open("rb") as f:
        nv_arr = np.fromfile(f, dtype=np.int64, count=1)
        if nv_arr.size != 1:
            raise ValueError(f"{filename}: could not read vector length header")
        nv = int(nv_arr[0])

        vec = np.fromfile(f, dtype=np.float64)
        if vec.size != nv:
            raise ValueError(
                f"{filename}: expected {nv} values, found {vec.size}"
            )

        return nv, np.ascontiguousarray(vec)


def read_mat(filename: Path) -> csr_matrix:
    """
    Read a CSR matrix file of the form:
        bool    real
        int64   nrow
        int64   ncol
        int64   nnz
        float64[nnz] rval      # or 2*nnz float64 for complex (real=false)
        int64  [nnz] col
        int64  [nrow+1] rowstt

    Returns a scipy.sparse.csr_matrix (float64 or complex128).
    """
    filename = Path(filename)
    with filename.open("rb") as f:
        # Note: np.bool_ is a single byte; matches typical C/C++ bool
        real_arr = np.fromfile(f, dtype=np.bool_, count=1)
        nrow_arr = np.fromfile(f, dtype=np.int64, count=1)
        ncol_arr = np.fromfile(f, dtype=np.int64, count=1)
        nnz_arr  = np.fromfile(f, dtype=np.int64, count=1)

        if any(x.size != 1 for x in (real_arr, nrow_arr, ncol_arr, nnz_arr)):
            raise ValueError(f"{filename}: malformed header")

        real = bool(real_arr[0])
        nrow = int(nrow_arr[0])
        ncol = int(ncol_arr[0])
        nnz  = int(nnz_arr[0])

        if nrow < 0 or ncol < 0 or nnz < 0:
            raise ValueError(f"{filename}: negative dimension in header")

        if real:
            vals = np.fromfile(f, dtype=np.float64, count=nnz)
            if vals.size != nnz:
                raise ValueError(f"{filename}: expected {nnz} values, got {vals.size}")
            vals = np.ascontiguousarray(vals)
            dtype = np.float64
        else:
            # File encodes complex as [Re, Im] pairs
            reim = np.fromfile(f, dtype=np.float64, count=2 * nnz)
            if reim.size != 2 * nnz:
                raise ValueError(
                    f"{filename}: expected {2*nnz} (re,im) values, got {reim.size}"
                )
            vals = np.ascontiguousarray(reim[0::2] + 1j * reim[1::2])
            dtype = np.complex128

        col = np.fromfile(f, dtype=np.int64, count=nnz)
        if col.size != nnz:
            raise ValueError(f"{filename}: expected {nnz} col indices, got {col.size}")

        rowstt = np.fromfile(f, dtype=np.int64, count=nrow + 1)
        if rowstt.size != nrow + 1:
            raise ValueError(
                f"{filename}: expected {nrow+1} row pointers, got {rowstt.size}"
            )

    # Basic CSR sanity checks
    if rowstt[0] != 0:
        raise ValueError(f"{filename}: rowstt[0] must be 0, got {rowstt[0]}")
    if rowstt[-1] != nnz:
        raise ValueError(f"{filename}: rowstt[-1] must be nnz ({nnz}), got {rowstt[-1]}")
    if (np.diff(rowstt) < 0).any():
        raise ValueError(f"{filename}: rowstt must be nondecreasing")
    if col.size and (col.min() < 0 or col.max() >= ncol):
        raise ValueError(f"{filename}: column index out of bounds [0, {ncol})")

    return csr_matrix((vals, col, rowstt), shape=(nrow, ncol), dtype=dtype)


# ------------------------------------------------------------
# HDF5 writer
# ------------------------------------------------------------
def save_system(
    out_path: Path,
    matrix: csr_matrix,
    rhs: Optional[NDArray] = None,
    sol: Optional[NDArray] = None,
    *,
    overwrite: bool = False,
    compression: Literal["gzip", "lzf", None] = "gzip",
    compression_opts: Optional[int] = 4,
) -> Path:
    """
    Persist a CSR matrix and optional vectors in HDF5.

    Layout:
      /matrix/{data,indices,indptr}  (attrs: shape, format, nnz, dtype, index_dtype, indptr_dtype, index_base=0)
      /rhs/values
      /sol/values
    """
    if not isinstance(out_path, Path):
        raise ValueError("out_path must be a pathlib.Path")
    if not isinstance(matrix, csr_matrix):
        raise ValueError("matrix must be a scipy.sparse.csr_matrix")

    # Ensure parent directory exists (fixes your FileNotFoundError)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    m, n = matrix.shape

    def _as_1d(name: str, x: Optional[NDArray], expected_len: int) -> Optional[NDArray]:
        if x is None:
            return None
        arr = np.asarray(x)
        if arr.ndim != 1:
            raise ValueError(f"{name} must be 1-D")
        if arr.shape[0] != expected_len:
            raise ValueError(f"{name} length must be {expected_len}, got {arr.shape[0]}")
        return np.ascontiguousarray(arr)

    rhs_arr = _as_1d("rhs", rhs, m)
    sol_arr = _as_1d("sol", sol, n)

    if out_path.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing file: {out_path}")

    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")

    data = np.ascontiguousarray(matrix.data)
    indices = np.ascontiguousarray(matrix.indices)
    indptr = np.ascontiguousarray(matrix.indptr)

    ds_kwargs: Dict[str, object] = {}
    if compression is not None:
        ds_kwargs["compression"] = compression
        if compression == "gzip" and compression_opts is not None:
            ds_kwargs["compression_opts"] = int(compression_opts)
        ds_kwargs["chunks"] = True  # let h5py choose

    with h5py.File(tmp_path, "w") as f:
        f.attrs["format-version"] = "1"
        f.attrs["producer"] = "convert-mat:save_system"

        g_matrix = f.create_group("matrix")
        g_matrix.create_dataset("data", data=data, **ds_kwargs)
        g_matrix.create_dataset("indices", data=indices, **ds_kwargs)
        g_matrix.create_dataset("indptr", data=indptr, **ds_kwargs)
        g_matrix.attrs["shape"] = (int(m), int(n))
        g_matrix.attrs["format"] = "csr"
        g_matrix.attrs["nnz"] = int(data.size)
        g_matrix.attrs["dtype"] = str(data.dtype)
        g_matrix.attrs["index_dtype"] = str(indices.dtype)
        g_matrix.attrs["indptr_dtype"] = str(indptr.dtype)
        g_matrix.attrs["index_base"] = 0  # explicit: 0-based indexing

        if rhs_arr is not None:
            g_rhs = f.create_group("rhs")
            g_rhs.create_dataset("values", data=rhs_arr, **ds_kwargs)
            g_rhs.attrs["length"] = int(rhs_arr.size)
            g_rhs.attrs["dtype"] = str(rhs_arr.dtype)

        if sol_arr is not None:
            g_sol = f.create_group("sol")
            g_sol.create_dataset("values", data=sol_arr, **ds_kwargs)
            g_sol.attrs["length"] = int(sol_arr.size)
            g_sol.attrs["dtype"] = str(sol_arr.dtype)

    # Atomic move into place
    tmp_path.replace(out_path)
    return out_path


# ------------------------------------------------------------
# Case discovery
# ------------------------------------------------------------
def find_cases(
    root: Path,
    *,
    require_complete: bool = True,
    resolve_conflicts: Literal["error", "first", "newest"] = "error",
    resolve_paths: bool = True,
) -> Dict[str, Dict[str, Path]]:
    """
    Group *.mat/*.rhs/*.sol triplets by common stem.
    """
    if not root.exists():
        raise FileNotFoundError(f"Root path does not exist: {root}")
    if not root.is_dir():
        raise NotADirectoryError(f"Root path is not a directory: {root}")

    part_by_ext = {".mat": "mat", ".rhs": "rhs", ".sol": "sol"}
    cases: Dict[str, Dict[str, Path]] = {}

    for p in root.rglob("*"):
        if not p.is_file():
            continue
        part = part_by_ext.get(p.suffix.lower())
        if part is None:
            continue
        case_name = p.stem
        target = p.resolve() if resolve_paths else p

        slot = cases.setdefault(case_name, {})
        if part in slot:
            if resolve_conflicts == "error":
                raise ValueError(
                    f"Duplicate {part!r} for case {case_name!r}:\n  {slot[part]}\n  {target}"
                )
            elif resolve_conflicts == "first":
                continue
            elif resolve_conflicts == "newest":
                try:
                    if target.stat().st_mtime > slot[part].stat().st_mtime:
                        slot[part] = target
                except FileNotFoundError:
                    slot[part] = target
        else:
            slot[part] = target

    if require_complete:
        cases = {
            name: parts
            for name, parts in cases.items()
            if {"mat", "rhs", "sol"} <= set(parts.keys())
        }

    return cases


# ------------------------------------------------------------
# Conversion driver (exposes options; sensible defaults)
# ------------------------------------------------------------
def convert_cases(
    input_root: Path,
    output_root: Path,
    *,
    require_complete: bool = True,
    resolve_conflicts: Literal["error", "first", "newest"] = "error",
    resolve_paths: bool = True,
    overwrite: bool = False,
    compression: Literal["gzip", "lzf", None] = "gzip",
    compression_level: int = 4,
    dry_run: bool = False,
) -> dict:
    """
    Convert discovered cases to HDF5 files.

    Returns a summary dict with counts and any per-case errors.
    """
    input_root = Path(input_root)
    output_root = Path(output_root)
    if not input_root.exists():
        raise FileNotFoundError(f"Input root does not exist: {input_root}")
    if not input_root.is_dir():
        raise NotADirectoryError(f"Input root is not a directory: {input_root}")

    out_summary = {
        "converted": 0,
        "skipped": 0,
        "errors": {},  # case_name -> str(error)
        "total_cases": 0,
        "output_dir": str(output_root),
    }

    cases = find_cases(
        input_root,
        require_complete=require_complete,
        resolve_conflicts=resolve_conflicts,
        resolve_paths=resolve_paths,
    )
    out_summary["total_cases"] = len(cases)

    # Ensure the output directory exists (extra safety; writer also ensures per-file parent)
    output_root.mkdir(parents=True, exist_ok=True)

    for name, parts in sorted(cases.items()):
        try:
            mat_path = parts.get("mat")
            rhs_path = parts.get("rhs")
            sol_path = parts.get("sol")
            if mat_path is None:
                raise ValueError("missing .mat")
            if rhs_path is None and require_complete:
                raise ValueError("missing .rhs")
            if sol_path is None and require_complete:
                raise ValueError("missing .sol")

            log.info("Converting %s", name)

            # Read inputs
            A = read_mat(mat_path)
            rhs_vec = None
            sol_vec = None

            if rhs_path is not None:
                rhs_n, rhs_vec = read_vec(rhs_path)
                if rhs_n != A.shape[0]:
                    raise ValueError(
                        f"{name}: rhs length {rhs_n} != nrow {A.shape[0]}"
                    )
            if sol_path is not None:
                sol_n, sol_vec = read_vec(sol_path)
                if sol_n != A.shape[1]:
                    raise ValueError(
                        f"{name}: sol length {sol_n} != ncol {A.shape[1]}"
                    )

            out_file = output_root / f"{name}.h5"
            sp.sparse.save_npz(output_root / f"{name}.npz", A, compressed=True)

            if dry_run:
                log.info("[dry-run] would write %s", out_file)
                out_summary["skipped"] += 1
                continue

            save_system(
                out_path=out_file,
                matrix=A,
                rhs=rhs_vec,
                sol=sol_vec,
                overwrite=overwrite,
                compression=compression,
                compression_opts=compression_level,
            )
            out_summary["converted"] += 1

        except Exception as e:  # keep going on errors
            log.error("Failed to convert %s: %s", name, e)
            out_summary["errors"][name] = str(e)

    return out_summary


# ------------------------------------------------------------
# CLI
# ------------------------------------------------------------
def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Convert .mat/.rhs/.sol systems to HDF5")
    p.add_argument("-i", "--input", type=Path, required=True, help="Input root directory")
    p.add_argument("-o", "--output", type=Path, required=True, help="Output root directory")

    p.add_argument("--allow-partial", action="store_true",
                   help="Include cases that are missing .rhs or .sol")
    p.add_argument("--conflict", choices=["error", "first", "newest"], default="error",
                   help="How to resolve duplicate parts of the same case")
    p.add_argument("--no-resolve", action="store_true",
                   help="Do not resolve paths (keep as discovered)")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing .h5 files")
    p.add_argument("--compression", choices=["gzip", "lzf", "none"], default="gzip",
                   help="HDF5 dataset compression")
    p.add_argument("--level", type=int, default=4,
                   help="Compression level (gzip only)")
    p.add_argument("--dry-run", action="store_true", help="Scan and read, but don't write")
    p.add_argument("-v", "--verbose", action="count", default=0,
                   help="Increase verbosity (-v, -vv)")

    return p.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    h5_conf = h5py.get_config()
    h5_conf.complex_names = ('real', 'imag')
    args = _parse_args(argv)

    # Logging
    level = logging.WARNING
    if args.verbose == 1:
        level = logging.INFO
    elif args.verbose >= 2:
        level = logging.DEBUG
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")

    compression = None if args.compression == "none" else args.compression  # type: ignore

    summary = convert_cases(
        input_root=args.input,
        output_root=args.output,
        require_complete=not args.allow_partial,
        resolve_conflicts=args.conflict,           # type: ignore[arg-type]
        resolve_paths=not args.no_resolve,
        overwrite=args.overwrite,
        compression=compression,                   # type: ignore[arg-type]
        compression_level=args.level,
        dry_run=args.dry_run,
    )

    converted = summary["converted"]
    skipped = summary["skipped"]
    errors = summary["errors"]
    total = summary["total_cases"]

    log.info("Total cases discovered: %d", total)
    log.info("Converted: %d, Skipped: %d, Errors: %d", converted, skipped, len(errors))
    if errors:
        log.warning("Some cases failed:")
        for k, v in errors.items():
            log.warning("  %s -> %s", k, v)


if __name__ == "__main__":
    main()
