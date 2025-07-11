/*************************************************************************************************/
/*                                                                                               */
/* Copyright (c) 2022 Rolls-Royce plc                                                            */
/*                                                                                               */
/* Redistribution and use in source and binary forms, with or without modification, are          */
/* permitted provided that the following conditions are met:                                     */
/*                                                                                               */
/* 1. Redistributions of source code must retain the above copyright notice, this list of        */
/*    conditions and the following disclaimer.                                                   */
/* 2. Redistributions in binary form must reproduce the above copyright notice, this list of     */
/*    conditions and the following disclaimer in the documentation and/or other materials        */
/*    provided with the distribution.                                                            */
/* 3. Neither the name of the copyright holder nor the names of its contributors may be used to  */
/*    endorse or promote products derived from this software without specific prior written      */
/*    permission.                                                                                */
/*                                                                                               */
/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS   */
/* OR IMPLIED  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF              */
/* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE    */
/* COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,     */
/* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE */
/* GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED    */
/* AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING     */
/* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED  */
/* OF THE POSSIBILITY OF SUCH DAMAGE.                                                            */
/*                                                                                               */
/*************************************************************************************************/
/**
 * \file    nozzle.h
 * \brief   Include file for nozzle typedefs and structures
 * \details Uses fixed width integer types to ensure interoperablility with Python
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    1st Dec 2022
 */

#ifndef NOZZLE_H
#define NOZZLE_H

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <limits.h>

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************
 *     Max and min functions                                                   *
 *******************************************************************************/
/**
 * \def MIN
 * Definition of MIN function
 * \def MAX
 * Definition of MAX function
 */
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))

/*******************************************************************************
 *     Definitions                                                             *
 *******************************************************************************/
/**
 * \def RHO_UPWIND_1ST
 * Moore and Moore 1st order density upwinding scheme
 * \def RHO_UPWIND_2ND
 * Moore and Moore 2nd order density upwinding scheme
 * \def RHO_UPWIND_3PT
 * Moore and Moore 3 point density upwinding scheme
 * \def RHO_UPWIND_JJG
 * Upwinding scheme from McGuirk and Page (AIAA J., 28, 10, Oct 1990)
 */
#define RHO_UPWIND_1ST 1
#define RHO_UPWIND_2ND 2
#define RHO_UPWIND_3PT 3
#define RHO_UPWIND_JJG 4


/*******************************************************************************
 *     Typdefs                                                                 *
 *******************************************************************************/
/**
 * \struct SETUP
 * \brief  Structure to hold setup parameters.
 *
 * \param  casename  casename which also provides root name for ouput files
 *                   in \ref FILES
 * \param  compress  false/true = incompressible/compressible flow
 * \param  couple    false/true = simple/coupled solver
 *
 * \author Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date   1st Dec 2022
 */

typedef struct
{
        char   *casename;
        bool   compress, couple;
} SETUP;

/**
 * \struct MESH
 * \brief  Structure to hold mesh parameters.
 * \details Only used for incompressible case where mesh is based on a circle
 *          passing through (0, ain), (0.5, ain-da), (1, ain).
 *          The mesh for incompressible case is set using the isentropic gas
 *          relations and the inlet Mach number,
 *          Note that number of mesh points is stored the \ref ARRAYS structure.
 *
 * \param  ain       inlet area for incompressible meshes. For compressible flows,
 *                   ain can be used as a scaling factor on the isentropic inlet area.
 * \param  da        area contraction/expansion at mid-point used for incompressible meshes
 *
 * \author Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date   1st Dec 2022
 */

typedef struct
{
        double ain, da;
} MESH;


/**
 * \struct FLOW
 * \brief  Structure to hold flow parameters. Unless stated variables are used by
 *         both incompressible and compressible cases. All variables are in 
 *         non-dimensional units.
 *
 * \param  rho       constant density for incompressible cases.
 * \param  ptin      inlet stagnation pressure
 * \param  ttin      inlet stagnation temperature (placeholder - currently hardwired to 1.0)
 * \param  min       inlet Mach number for compressible cases
 * \param  pex       exit static pressure
 * \param  gamma     ratio of specific heats (usually 1.4 for a ideal gas)
 * \param  cp        specific heat capacity at constant pressure (usually 3.5 in
 *                   non-dimensional units)
 * \param  rgas      ideal gas constant (1.0 in non-dimensional units)
 * \param  upwind    density upwinding scheme for compressible flows. Use 1, 2 and 3 in
 *                   input file for \ref RHO_UPWIND_1ST, \ref RHO_UPWIND_2ND and
 *                   \ref RHO_UPWIND_3PT respectively.
 *
 * \author Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date   1st Dec 2022
 */

typedef struct
{
        double  rho, ptin, ttin, min, pex;
        double  gamma, cp, rgas;
        int32_t upwind;
} FLOW;

/**
 * \struct SOLVER
 * \brief  Structure to solver parameters
 *
 * \param  niters    number of outer iterations
 * \param  relaxu    velocity relaxation factor for SIMPLE solver
 * \param  relaxp    pressure relaxation factor for SIMPLE solver
 * \param  relaxc    relaxation factor for coupled solver (same for all equations)
 * \param  presim    number of simple iterations prior to coupled solver
 * \param  rrampc    number of coupled iterations over which to ramp the relaxation
 *                   factor. The latter increases linear from 1/10th its value to
 *                   the full value. Usually only needed for compressible cases.
 * \param  toleru    convergence tolerance for SIMPLE momentum equation
 * \param  tolerp    convergence tolerance for SIMPLE pressure correction equation
 * \param  tolerc    convergence coupled equations
 *
 * \author Leigh Lapworth <leigh.lapworth@btinternet.com>
 * \date   5th Dec 2022
 */

typedef struct
{
        int32_t niters;
        double  relaxu, relaxp, relaxc;
        int32_t presim, rrampc;
        double  toleru, tolerp, tolerc;
} SOLVER;

/**
 * \struct FILES 
 * \brief   Structure to hold file names derived from case name
 * \details All files have length NAME_MAX from the system file <limits.h>
 *
 * \param  histfile  name of convergence history file
 * \param  solnfile  name of solution file 
 *
 * \author Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date   5th Dec 2022
 */

typedef struct
{
        char histfile[NAME_MAX];
        char solnfile[NAME_MAX];
} FILES;

/**
 * \struct ARRAYS
 * \brief  Structure to hold arrays
 *
 * \param  nx        number if axial stations in the nozzle mesh
 * \param  a         flow area at station i
 * \param  rho       density at station i
 * \param  u         velocity at station i
 * \param  p         pressure at station i for colocated grids and midway between
 *                   stations for staggered grids
 * \param  t         temperature at station i
 * \param  M         Mach number at station i
 * \param  apu       used to hold momentum equation coefficients to be used in
 *                   the pressure correction equation.
 * \param  peff      used to hold effective pressure for density upwinding
 * \param  a0,a1,a2  coefficients in upwinding density scheme
 *
 * \author Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date   1st Dec 2022
 */

typedef struct
{
        int32_t nx;
        double  *a;
        double  *rho;
        double  *u;
        double  *p;
        double  *t;
        double  *M;
        double  *apu;
        double  *peff, *a0, *a1, *a2;
} ARRAYS;

/**
 * \struct HISTORY
 * \brief  Structure to hold convergence history for current iteration
 *
 * \param  iter      counter for current iteration (counts from 0)
 * \param  dr_rms    RMS of changes to rho (compressible flow)
 * \param  du_rms    RMS of changes to u
 * \param  dp_rms    RMS of changes to p
 * \param  resr_rms  RMS of residuals upwinded gas law (compressible flow)
 * \param  resu_rms  RMS of residuals for momentum equation
 * \param  resp_rms  RMS of residuals for pressure correction equation
 *
 * \author Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date   1st Dec 2022
 */

typedef struct
{
        int32_t iter;
        double  dr_rms,   du_rms,   dp_rms;
        double  resr_rms, resu_rms, resp_rms;
} HISTORY;

/**
 * \struct NOZZLE
 * \brief   Master structure to hold all nozzle data
 * \details Used to minimise the number of arguments to be passed between modules,
 *          particularly across C-Python interface.
 *
 * \param  setup     \ref SETUP data
 * \param  mesh      \ref MESH data
 * \param  flow      \ref FLOW data
 * \param  solver    \ref SOLVER data
 * \param  arrays    \ref ARRAYS data
 * \param  files     \ref FILES data
 *
 * \author Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date   1st Dec 2022
 */

typedef struct
{
        SETUP   *setup;
        MESH    *mesh;
        FLOW    *flow;
        SOLVER  *solver;
        FILES   *files;
        ARRAYS  *arrays;
        HISTORY *history;
} NOZZLE;

/**
 * \def NEW_NOZZLE
 * Pseudo constructor - default values for initialising \ref NOZZLE structure
 */

#define NEW_NOZZLE  {NULL, NULL, NULL, NULL, NULL, NULL, NULL}

/**
 * \struct  SMATRIX 
 * \brief   Structure to hold sparse matrix
 * \details The compressed sparse row (CSR) structure is used to store the matrix.
 *          The structure also stores the right hand side vector for solving Ax=b.
 *          Storing b enables easier passing of data across the C-Python interface
 *          as the coefficients and the RHS vector are usually evaluated in the
 *          same routine to ensure a consistent treatment.
 *
 * \param   nr  number of rows
 * \param   nc  number of columns
 * \param   nnz number of non-zero entries
 * \param   col column index of each matrix entry (counted from zero)
 * \param   rst index of the start of each row (counted from zero)
 * \param   val real valued matrix entries
 * \param   b   right hand side vector for solving Ax=b
 *
 * \author Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date   20th Dec 2022
 */

typedef struct
{
        int32_t nr;
        int32_t nc;
        int32_t nnz;
        int32_t *col;
        int32_t *rst;
        double  *val;
        double  *b;
} SMATRIX;

/**
 * \def NEW_SMAT
 * Pseudo constructor - default values for initialising \ref SMATRIX structure
 */

#define NEW_SMAT {0, 0, 0, NULL, NULL, NULL, NULL}



#ifdef __cplusplus
}
#endif

#endif  /* NOZZLE_H */

