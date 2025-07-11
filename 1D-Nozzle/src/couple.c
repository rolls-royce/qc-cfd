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
 * \file    couple.c
 * \brief   Routines for coupled implicit solver
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    29th Dec 2022
 */

/******************************************************************************
 *  Include files                                                             *
 ******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include <proto.h>

/******************************************************************************
 *  couple_update - update rho,u,p vectors from solution of coupled system    *
 ******************************************************************************/
/**
 * \brief   Update rho,u,p vectors using x from the solution of coupled Ax=b
 * \details Used for updating primitive variables after coupled solver
 *          For compressible cases dx contains (du, dp, drho), 
 *          For incompressible cases dx contains (du, dp),
 *
 * \param[in,out]
 *          nozzle  \ref NOZZLE structure containing rho,u,p
 * \param[in]
 *          dx      vector containing solution to coupled linear system
 *                  A.dx = (b - A.x).
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    30th Dec 2022
 */
int  couple_update(NOZZLE *nozzle, double *dx)
{
     SETUP   *setup   = nozzle->setup;
     SOLVER  *solver  = nozzle->solver;
     ARRAYS  *arrays  = nozzle->arrays;
     HISTORY *history = nozzle->history;

     bool    compress = setup->compress;
     int32_t nx       = arrays->nx;
     double  *rho     = arrays->rho;
     double  *u       = arrays->u;
     double  *p       = arrays->p;
     double  relax    = solver->relaxc;
     double  rramp    = solver->rrampc;

     int32_t       i, roff, uoff, poff;
     static double rfac, dfac;
     static bool   first = true;
     double        du_rms, dp_rms, dr_rms;
     int           rc;

/*   check for NaN
     ------------- */
     if (compress) rc = nan_check("dx", 3*nx, dx);
     else          rc = nan_check("dx", 2*nx, dx);

     if(rc){
       error_traceback("couple_update", __FILE__, __LINE__);
       return 1;
     }

/*   set relaxation factor ramp
     -------------------------- */
     if (rramp > 0){
       if (first){
         rfac  = 0.1 * relax;
         dfac  = relax/(double) rramp;
         first = false;
       }
       else{
         rfac = MIN(1.0, rfac + dfac);
       }
       relax = rfac*relax;
     }

/*   offsets
     ------- */
     uoff = 0;
     poff = uoff + nx;
     roff = poff + nx;

/*   calculate RMS changes
     --------------------- */
     du_rms = 0.0;
     dp_rms = 0.0;

     for(i=0; i<nx; i++){
        du_rms += dx[uoff+i] * dx[uoff+i];
        dp_rms += dx[poff+i] * dx[poff+i];
     }
     du_rms = sqrt(du_rms/(nx-1));
     dp_rms = sqrt(dp_rms/(nx-1));

     history->du_rms = du_rms;
     history->dp_rms = dp_rms;

     if(compress){
       dr_rms = 0.0;
       for(i=0; i<nx; i++){
          dr_rms += dx[roff+i] * dx[roff+i];
       }
       dr_rms = sqrt(dr_rms/(nx-1));

       history->dr_rms = dr_rms;
     }

/*   update fields
     ------------- */
     if (compress) {
       uoff = 0;
       poff = uoff + nx;
       roff = poff + nx;
       for(i=0; i<nx; i++){
         u[i]   += relax*dx[uoff+i];
         p[i]   += relax*dx[poff+i];
         rho[i] += relax*dx[roff+i];
       }
     }
     else{
       uoff = 0;
       poff = uoff + nx;
       for(i=0; i<nx; i++){
         u[i] += dx[uoff+i];
         p[i] += dx[poff+i];
       }
     }

     return 0;
}

/******************************************************************************
 *  couple_mat - create sparse matrix for 1D nozzle coupled solver            *
 ******************************************************************************/
/**
 * \brief   Create coupled matrix and resiual vector for 1D nozzle.
 * \details Calls \ref couple_mat_i for incompressible nozzle and
 *          \ref couple_mat_c for compressible nozzle.
 *          Note, right hand side residual is contained within the
 *          \ref SMATRIX structure.
 *
 * \param[in]
 *          nozzle  \ref NOZZLE structure containing nozzle flow and geometry
 * \param[out]
 *          S       \ref SMATRIX containing linearised coupled system and residual
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    29th Dec 2022
 */

int  couple_mat(NOZZLE *nozzle, SMATRIX *S)
{
     SETUP *setup = nozzle->setup;
     int   rc;

     if (setup->compress) rc = couple_mat_c(nozzle, S);
     else                 rc = couple_mat_i(nozzle, S);

     return rc;
}

/******************************************************************************
 *  couple_mat_c - create sparse matrix for compressible coupled solver       *
 ******************************************************************************/
/**
 * \brief   Create coupled matrix and resiual vector for compressible 1D
 *          nozzle.
 * \details Uses only primitive variables from the \ref NOZZLE structure.
 *
 * \param[in]
 *          nozzle  \ref NOZZLE structure containing cavity flow and geometry
 * \param[out]
 *          S       \ref SMATRIX containing linearised coupled system
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    29th Dec 2022
 */

int  couple_mat_c(NOZZLE *nozzle, SMATRIX *S)
{
     SETUP   *setup   = nozzle->setup;
     SOLVER  *solver  = nozzle->solver;
     FLOW    *flow    = nozzle->flow;
     ARRAYS  *arrays  = nozzle->arrays;
     HISTORY *history = nozzle->history;

     int32_t nx    = arrays->nx;
     double  *rho  = arrays->rho;
     double  *u    = arrays->u;
     double  *p    = arrays->p;
     double  *a    = arrays->a;
     double  *t    = arrays->t;
     double  *M    = arrays->M;
     double  *peff = arrays->peff;
     double  *a0up = arrays->a0;
     double  *a1up = arrays->a1;
     double  *a2up = arrays->a2;

     double  gamma = flow->gamma;
     double  cp    = flow->cp;
     double  rgas  = flow->rgas;
     double  ttin  = flow->ttin;

     int32_t *col, *rst;
     double  *val, *b, *rhou;

     int32_t i, nnz, iz, nrc, irow, icol, uoff, poff, roff;
     double  a0, a1, t0, resu, resp, resr;
     bool   rup2ndo = false;          // use 2nd order upwinding coeffs in matrix - less stable

/*   free incoming matrix
     -------------------- */
     sparse_free(S);

/*   set temperature and get effective pressure
     ------------------------------------------ */
     noz_temp_mach(nozzle);
     noz_peff(nozzle);

/*   set rho*u
     --------- */
     rhou    = (double*) malloc(nx * sizeof(double));
     rhou[0] = rho[0]*u[0];
     for (i=1; i<nx; i++)
       rhou[i] = 0.5*(rho[i]*u[i] + rho[i-1]*u[i-1]);

/*   set total number of rows and columns
     ------------------------------------ */
     nrc  = 3*nx;
     uoff = 0;
     poff = uoff + nx;
     roff = poff + nx;

/*   allocate memory for coupled coefficient matrix S
     ------------------------------------------------ */
     nnz  = 4*(nx-1) + 2;               // momentum
     nnz += 4*(nx-1) + 1;               // continuity
     nnz += 4* nx - 1;                  // gas law (1st order upwinded density)
     if(rup2ndo) nnz += nx - 2;         // second order upwinding
     iz   = 0;

     S->nr  = nrc;
     S->nc  = nrc;
     S->nnz = nnz;
     S->col = (int32_t*) malloc(nnz     * sizeof(int32_t));
     S->rst = (int32_t*) malloc((nrc+1) * sizeof(int32_t));
     S->val = (double*)  malloc(nnz     * sizeof(double));
     S->b   = (double*)  malloc(nrc     * sizeof(double));

     rst = S->rst;
     col = S->col;
     val = S->val;
     b   = S->b;

/*   inlet bc
     -------- */
     t0        = pow(p[0],(gamma-1.0)/gamma);
     irow      = 0;
     rst[irow] = iz;

     col[iz] = uoff + irow;
     val[iz] = u[0]/cp;
     iz ++;

     col[iz] = poff + irow;
     val[iz] = (1.0 - 1.0/gamma)*pow(p[0],-1.0/gamma);
     iz ++;
     b[uoff + irow] = ttin - t0 - 0.5*u[0]*u[0]/cp;

/*   momentum - Gupta's convective volume discretisation 
     --------------------------------------------------- */
     for (i=1; i<nx; i++){
       irow      = i;
       rst[irow] = iz;

       col[iz] = uoff + irow - 1;
       val[iz] = -rhou[i];
       iz ++;

       col[iz] = uoff + irow;
       val[iz] = rhou[i];
       iz ++;

       col[iz] = poff + irow - 1;
       val[iz] = -1;
       iz ++;

       col[iz] = poff + irow ;
       val[iz] = 1;
       iz ++;

/*     residual
       -------- */
       b[uoff + irow] = -rhou[i]*(u[i] - u[i-1]) - (p[i] - p[i-1]);
     }

/*   continuity
     ---------- */
     for (i=0; i<nx-1; i++) {
       irow             = i;
       rst[poff + irow] = iz;

       col[iz] = uoff + irow;
       val[iz] = a[i]*rho[i];
       iz ++;

       col[iz] = uoff + irow+1;
       val[iz] = -a[i+1]*rho[i+1];
       iz ++;

       col[iz] = roff + irow;
       val[iz] = a[i]*u[i];
       iz ++;

       col[iz] = roff + irow+1;
       val[iz] = -a[i+1]*u[i+1];
       iz ++;

/*     residual
       -------- */
       b[poff + irow] = -(a[i]*rho[i]*u[i] - a[i+1]*rho[i+1]*u[i+1]);
     }

/*   oulet bc
     -------- */
     irow             = nx-1;
     rst[poff + irow] = iz;

     col[iz] = poff + irow;
     val[iz] = a[nx-1]*rho[nx-1];
     iz ++;
     b[poff + irow] = 0.0;

/*   gas law - inlet
     --------------- */
     i = 0;
     irow             = i;
     rst[roff + irow] = iz;

     col[iz] = uoff + irow;
     val[iz] = -rgas*rho[i]*u[i]/cp;
     iz ++;

     col[iz] = poff + irow;
     val[iz] = -1;
     iz ++;

     col[iz] = roff + irow;
     val[iz] = rgas*t[i];
     iz ++;
     b[roff + irow] = (p[i] - rgas*t[i]*rho[i]);

/*   gas law - interior 1st and 2nd order (if rup2ndo is true, this is less stable)
     ------------------------------------------------------------------------------ */
     for (i=1; i<nx; i++) {
       a0 = a0up[i];
       if(rup2ndo) a1 = a1up[i];
       else        a1 = 0.0;

       irow             = i;
       rst[roff + irow] = iz;

       col[iz] = uoff + irow;
       val[iz] = -rgas*rho[i]*u[i]/cp;
       iz ++;

       if(rup2ndo && i>1){
         col[iz] = poff + irow - 2;
         val[iz] = 0.5*a1;
         iz ++;
       }
       col[iz] = poff + irow - 1;
       val[iz] = -(1.0 - a0);
       iz ++;
       col[iz] = poff + irow;
       val[iz] = -a0 - 0.5*a1;
       iz ++;

       col[iz] = roff + irow;
       val[iz] = rgas*t[i];
       iz ++;

       b[roff + irow] = (peff[i] - rgas*t[i]*rho[i]);
     }

/*   terminate row list
     ------------------ */
     rst[nrc] = iz;

/*   free local arrays
     ----------------- */
     free(rhou);

/*   check
     ----- */
     if(iz != nnz){
       error_message("couple_mat_c",  __FILE__, __LINE__,
                     "mis-match in non-zeros: iz = %d, nnz = %d", iz, nnz);
       return 1;
     }

/*   calculate RMS residuals
     ----------------------- */
     resu = 0.0;
     resp = 0.0;
     resr = 0.0;

     for(i=0; i<nx; i++){
       resu += b[uoff+i] * b[uoff+i];
       resp += b[poff+i] * b[poff+i];
       resr += b[poff+i] * b[poff+i];
     }
     resu = sqrt(resu/(nx-1));
     resp = sqrt(resp/(nx-1));
     resr = sqrt(resr/(nx-1));

     history->resu_rms = resu;
     history->resp_rms = resp;
     history->resr_rms = resr;

     return 0;
}

/******************************************************************************
 *  couple_mat_i - create sparse matrix for incompressible coupled solver     *
 ******************************************************************************/
/**
 * \brief   Create coupled matrix and resiual vector for incompressible 1D
 *          nozzle.
 * \details Uses only primitive variables from the \ref NOZZLE structure.
 *
 * \param[in]
 *          nozzle  \ref NOZZLE structure containing cavity flow and geometry
 * \param[out]
 *          S       \ref SMATRIX containing linearised coupled system
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    29th Dec 2022
 */

int  couple_mat_i(NOZZLE *nozzle, SMATRIX *S)
{
     SETUP   *setup   = nozzle->setup;
     SOLVER  *solver  = nozzle->solver;
     FLOW    *flow    = nozzle->flow;
     ARRAYS  *arrays  = nozzle->arrays;
     HISTORY *history = nozzle->history;

     int32_t nx       = arrays->nx;
     double  *rho     = arrays->rho;
     double  *u       = arrays->u;
     double  *p       = arrays->p;
     double  *a       = arrays->a;
     double  ptin     = flow->ptin;

     int32_t *col, *rst;
     double  *val, *b, *rhou;

     int32_t i, nnz, iz, nrc, irow, icol, uoff, poff;
     double  resu, resp;
     int     rc;

/*   free incoming matrix
     -------------------- */
     sparse_free(S);

/*   set rho*u
     --------- */
     rhou    = (double*) malloc(nx * sizeof(double));
     rhou[0] = rho[0]*u[0];
     for (i=1; i<nx; i++)
       rhou[i] = 0.5*(rho[i]*u[i] + rho[i-1]*u[i-1]);

/*   set total number of rows and columns
     ------------------------------------ */
     nrc  = 2*nx;
     uoff = 0;
     poff = uoff + nx;

/*   allocate memory for coupled coefficient matrix S
     ------------------------------------------------ */
     nnz  = 4*(nx-1) + 2;               // interior udu/dx, dp/dx = 4 + inlet bc
     nnz += 2*(nx-1) + 1;               // interior udu = 4 + outlet bc
     iz   = 0;                          // initialise counter

     S->nr  = nrc;
     S->nc  = nrc;
     S->nnz = nnz;
     S->col = (int32_t*) malloc(nnz     * sizeof(int32_t));
     S->rst = (int32_t*) malloc((nrc+1) * sizeof(int32_t));
     S->val = (double*)  malloc(nnz     * sizeof(double));
     S->b   = (double*)  malloc(nrc     * sizeof(double));

     rst = S->rst;
     col = S->col;
     val = S->val;
     b   = S->b;

/*   inlet bc
     -------- */
     irow      = 0;
     rst[irow] = iz;

     col[iz] = uoff + irow;
     val[iz] = rhou[0];
     iz ++;

     col[iz] = poff + irow;
     val[iz] = 1;
     iz ++;
     b[uoff + irow] = ptin - p[0] - 0.5*rho[0]*u[0]*u[0];

/*   momentum - use Gupta discretisation for compatibility with SIMPLE
     ----------------------------------------------------------------- */
     for (i=1; i<nx; i++){
       irow      = i;
       rst[irow] = iz;

       col[iz] = uoff + irow - 1;
       val[iz] = -rhou[i];
       iz ++;

       col[iz] = uoff + irow;
       val[iz] = rhou[i];
       iz ++;

       col[iz] = poff + irow - 1;
       val[iz] = -1;
       iz ++;

       col[iz] = poff + irow ;
       val[iz] = 1;
       iz ++;

/*     residual
       -------- */
       b[uoff + irow] = -(rhou[i]*u[i] -rhou[i]*u[i-1] + p[i] - p[i-1]);
     }

/*   continuity
     ---------- */
     for (i=0; i<nx-1; i++) {
       irow             = i;
       rst[poff + irow] = iz;

       col[iz] = uoff + irow;
       val[iz] = a[i]*rho[i];
       iz ++;

       col[iz] = uoff + irow+1;
       val[iz] = -a[i+1]*rho[i+1];
       iz ++;

/*     residual
       -------- */
       b[poff + irow] = -(a[i]*rho[i]*u[i] - a[i+1]*rho[i+1]*u[i+1]);
     }

/*   oulet bc
     -------- */
     irow             = nx-1;
     rst[poff + irow] = iz;

     col[iz] = poff + irow;
     val[iz] = a[nx-1]*rho[nx-1];
     iz ++;
     b[poff + irow] = 0.0;

/*   terminate row list
     ------------------ */
     rst[nrc] = iz;

/*   free local arrays
     ----------------- */
     free(rhou);

/*   check
     ----- */
     if(iz != nnz){
       error_message("couple_mat_i",  __FILE__, __LINE__,
                     "mis-match in non-zeros: iz = %d, nnz = %d", iz, nnz);
       return 1;
     }

/*   calculate RMS residuals
     ----------------------- */
     resu = 0.0;
     resp = 0.0;

     for(i=0; i<nx; i++){
       resu += b[i]    * b[i];
       resp += b[nx+i] * b[nx+i];
     }
     resu = sqrt(resu/(nx-1));
     resp = sqrt(resp/(nx-1));

     history->resu_rms = resu;
     history->resp_rms = resp;

     return 0;
}

/******************************************************************************
 *  couple_run - run coupled implicit solver                                  *
 ******************************************************************************/
/**
 * \brief   Apply coupled implicit solver
 * \details Self contained solver that can be called with just a \ref NOZZLE
 *          instance
 *
 * \param[in]
 *          nozzle    \ref NOZZLE instance
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    29th Dec 2022
 */

int  couple_run(NOZZLE *nozzle)
{
     SETUP   *setup   = nozzle->setup;
     SOLVER  *solver  = nozzle->solver;
     FILES   *files   = nozzle->files;
     HISTORY *history = nozzle->history;

     SMATRIX S = NEW_SMAT;

     int32_t niters = solver->niters;
     int32_t iter;
     int32_t nsimple = solver->presim;
     double  *dx;
     int     rc;

/*   check setup
     ----------- */
     if(!setup->couple){
       error_message("couple_run", __FILE__, __LINE__, "%s", "solver set for SIMPLE");
       return 1;
     }

/*   for compressible flow perform initial SIMPLE iterations for stability
     --------------------------------------------------------------------- */
     if(setup->compress && nsimple>0){
       setup->couple  = 0;
       solver->niters = nsimple;

       simple_run(nozzle);

       setup->couple  = 1;
       solver->niters = niters - nsimple;
       niters = solver->niters;
     }

/*   solve coupled system
     -------------------- */
     for (iter=1; iter<=niters; iter++){
       if((rc = couple_mat(nozzle, &S))){
         error_traceback("couple_run", __FILE__, __LINE__);
         return 1;
       }

       if((rc = matsol_LU(&S, &dx))){
         error_traceback("couple_run", __FILE__, __LINE__);
         return 1;
       }

       if((rc = couple_update(nozzle, dx))){
         error_traceback("couple_run", __FILE__, __LINE__);
         return 1;

       }
       sparse_free(&S);
       free(dx);

       hist_print(nozzle);
       if(hist_check_conv(nozzle)) break;
       history->iter++;
     }

/*   save solution using casename
     ---------------------------- */
     print_flow(nozzle, files->solnfile);

     return 0;
}

