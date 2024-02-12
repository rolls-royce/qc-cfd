# Running L-QLES

L-QLES is an open source Python code for generating 1D, 2D and 3D Laplacian
operators and associated Poisson equations and their classical solutions.
The Laplacians are created using a finite volume discretisation
on Cartesian lattice meshes.
The main feature of L-QLES is the ability to _tune_ the
the mesh and, hence, the Laplacian to include the following features:

* Non-uniform mesh distributions,
* Multiple boundary condition types,
* Arbitrary mesh indexing.

## Launching L-QLES

The general syntax is:

`````
l-qles.py -i <input file> {-c <x,y,z>} {-d} {-e} {-h} {-j} {-m} {-r} {-s}

     -i {name of input file}
     -c {x,y,z} cut slice of 3D solution to be plotted, default = x
     -d allow degnerate matrices, default = False
     -e calculate eigenvalues and condition number, default = False
     -h help menu
     -j split plots into separate windows for saving, default is single window
     -m plot matrix, default = False
     -r reorder matrix and RHS to use shell ordering of mesh, default = False
     -s plot solutons and mesh, default = False
`````

## Commad line options
The default for all options is _off_ unless otherwise stated.
* __-i__  Followed by the name of the XML input file. There is no default.
* __-c__ Followed by _x_, _y_or _z_. For 3D cases this indicates which cutting
     plane is used to show the solution if the option __-s__ is turned on.
     The cutting plane is position at the mid-point of the domain and default is an _x_
     plane.

* __-d__ Laplacians that consist entirely of
     repeating and/or Neumann boundaries are degenerate. The default is to remove the degeneracy
     by applying a Dirichlet condition at a single point in the mesh determined by the input
     variable _degfix_.
     Since other authors have used the degenerate form, this option does not apply the
     degeneracy fix so that like with like comparisons can be made. Note that L-QLES cannot
     solve a degenerate Poisson equation and no solution file is stored in these cases.
     Other files have __d_ appended to their case name.

* __-e__ Calculate the eigenvalues and condition number of the Laplacian. This scales
     with $O(N^3)$ where $N$ is the dimension of the matrix and so can only be used with
     small matrices. 

* __-h__ Display help menu.

* __-j__ The default for the __-m__ and __-s__ plotting options is to
     create figures with 2 plots per pane. This option plots each figure separately which
     may be useful when preparing reports.

* __-m__ Plot the matrix values and sparsity pattern using Matplotlib.

* __-r__ Applies shell reordering described. If this
     option is on, all files have __r_  appended to their case name.
     If the dengeneracy and reordering options are both on, then __d_r_ is appended to
     the case name.

* __-s__ Plot the mesh and contours of the solution variable.



## Input file format 

The mesh and boundary conditions are set via a single XML input file.

`````
<?xml version="1.0" encoding="UTF-8"?>
<laplace>
  <case name="l1d_16_dd" dimension="1" force="1.0"></case>
  <mesh direction="x">
    <length>1.0</length>
    <ntotal>16</ntotal>
    <nclust>6</nclust>
    <cltype>2</cltype
    <cratio>1.2</cratio>
    <btype>D, D</btype>
    <bvalue>0.0, 0.0</bvalue>
    <degfix>8</degfix>
  </mesh>
</laplace>
`````

The above listing shows the input file for a 1D mesh.
2D and 3D meshes are created by changing _dimension_ in the third line and
adding the equivalent _mesh_ sections for the "y" and "z" directions.

## Output files
Note the Laplacian, $L$, is normalised to have $||L||_{max}=1$ with the
same scale factor applied to the RHS state to ensure that the
solution state corresponds to the original problem. This does not
mean that the RHS state is normalised.

L-QLES outputs 2 types of files: Python and C/C++ compatible binary files:

* __Laplacian__ This is stored using the compressed sparse row format.
    The name of the Python file is _casename_mat.npz_ and the
    C/C++ binary file is _casename_mat.bin_.

* __RHS vector__ The right-hand side vector contains the boundary values and
    the bulk inhomogeneous force term.
    The name of the Python file is _casename_rhs.npy_ and the
    C/C++ binary file is _casename_rhs.bin_.

* __Solution vector__If the Laplacian is not degenerate, the solution
    vector is output.
    The name of the Python file is _casename_sol.npy_ and the
    C/C++ binary file is _casename_sol.bin_.

* __Reordering matrix__ If the reordering option has been used then
    the column-wise permutation matrix  is written
    to the files _casename_ord.npz_ and _casename_ord.bin_.
    This matrix is needed to recover the
    solution to the original Laplacian from the solution to the reordered one.

If __-r__  and/or __-d__ options have been used, the case name is
amended as described above.

