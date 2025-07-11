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
 * \file    simple.c
 * \brief   Routines for SIMPLE pressure correction solver
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    8th Dec 2022
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
 *  simple_msolve - solve momentum equation                                   *
 ******************************************************************************/
/**
 * \brief   Solve momentum equation
 * \details Solve momentum using simple marching (hyperbolic equation), 
 *          discretisation using colocated grid with opposed differencing.
 *          This follows "Development of a New Shock Capturing Formula for
 *          Pressure Correction Methods" by Ajay K. Gupta, see
 *          \href https://vtechworks.lib.vt.edu/bitstream/handle/10919/46310/LD5655.V855_1993.G867.pdf?sequence=1&isAllowed=y
 *          Gupta's convective volume integral is also used as the more usual 
 *          conservation discretisation mostly unstable. Similarly, the Jacobi
 *          solver used by Gupta is also used as Gauss-Seidel is not stable for
 *          compressible case.
 *
 * \param[in, out]
 *          nozzle  \ref NOZZLE structure to hold nozzle parameters.
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    8th Dec 2022
 */

int  simple_msolve(NOZZLE *nozzle)
{
     SOLVER  *solver  = nozzle->solver;
     ARRAYS  *arrays  = nozzle->arrays;
     HISTORY *history = nozzle->history;
     int32_t  nx       = arrays->nx;
     double   *rho     = arrays->rho;
     double   *u       = arrays->u;
     double   *p       = arrays->p;
     double   *a       = arrays->a;
     double   *apu     = arrays->apu;
     double   relaxu   = solver->relaxu;
     int32_t  i;
     double   *resu, *du;
     double   du_rms, resu_rms;
     bool     debug = false;

/*   alloc local arrays
     ------------------ */
     resu  = (double*) malloc(nx*sizeof(double));
     du    = (double*) malloc(nx*sizeof(double));
     memset(resu, 0, nx*sizeof(double));
     memset(du,   0, nx*sizeof(double));

/*   conservative discretisation
     --------------------------- */
//   for (i=0; i<nx; i++) apu[i] = a[i]*rho[i]*u[i]; 
//
//   for (i=1; i<nx; i++){
//     resu[i] = -(apu[i]*u[i] - apu[i-1]*u[i-1]) - (a[i]*p[i] - a[i-1]*p[i-1]);
//     du[i]   = resu[i]/apu[i];
//   }

/*   Gupta's convective volume discretisation and Jacobi solver
     ---------------------------------------------------------- */
     du[0] = 0.0;
     for (i=1; i<nx; i++)
       apu[i] = 0.5*(rho[i]*u[i] + rho[i-1]*u[i-1]);
   
     for (i=1; i<nx; i++){
       resu[i] = -apu[i]*(u[i] - u[i-1]) - (p[i] - p[i-1]);
       du[i]   = resu[i]/apu[i];
     }

/*   solve for u
     ----------- */
     for (i=1; i<nx; i++) u[i] = u[i] + relaxu*du[i];

/*   debug
     ----- */
     if(debug){
       printf(" i,  area,      rho,       p,          apu,       resu,      u          du\n");
       for (i=0; i<nx; i++){
         printf("%2d, % lf, % lf, % lf,  % lf, % lf, % lf, % lf\n",
                 i, a[i], rho[i], p[i], apu[i], resu[i], u[i], du[i]);
       }
     }

/*   rms residual
     ------------ */
     du_rms   = 0.0;
     resu_rms = 0.0;
     for (i=1; i<nx; i++){
       du_rms   += du[i]  *du[i];
       resu_rms += resu[i]*resu[i];
     }
     du_rms   = sqrt(du_rms/(nx-1));
     resu_rms = sqrt(resu_rms/(nx-1));

/*   save residuals
     -------------- */
     history->du_rms   = du_rms;
     history->resu_rms = resu_rms;

/*   free local arrays
     ----------------- */
     free(resu);
     free(du);

     return 0;
}

/******************************************************************************
 *  simple_pcoeff - get coefficients for pressure corrrection equation        *
 ******************************************************************************/
/**
 * \brief   Get coefficients for pressure corrrection equation
 * \details Returns matrix coefficients aw, ap and ae, and continuity residual.
 *          Use \ref simple_pcmat to convert into a sparse matrix.
 *          Calling routine must allocate memory for coeffs and residual, all
 *          arrays to be doubles of dimension nx.
 *          As well as setting cellwise continuity residuals, the RMS residual
 *          is stored in the nozzle's \ref HISTORY structure.
 *
 * \param[in]
 *          nozzle  \ref NOZZLE structure holding nozzle parameters and arrays.
 * \param[in, out]
 *          aw      coefficent for point i-1 in pc equation for i.
 * \param[in, out]
 *          ap      coefficent for point i in pc equation for i.
 * \param[in, out]
 *          ae      coefficent for point i+1 in pc equation for i.
 * \param[in, out]
 *          bb      continuity residual
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    20th Dec 2022
 */

int  simple_pcoeff(NOZZLE *nozzle,
                   double *aw, double *ap, double *ae, double *bb)
{
     SETUP   *setup   = nozzle->setup;
     FLOW    *flow    = nozzle->flow;
     ARRAYS  *arrays  = nozzle->arrays;
     HISTORY *history = nozzle->history;

     int32_t nx       = arrays->nx;
     double  *rho     = arrays->rho;
     double  *u       = arrays->u;
     double  *t       = arrays->t;
     double  *a       = arrays->a;
     double  *apu     = arrays->apu;

     double  gamma    = flow->gamma;
     double  cp       = flow->cp;
     double  rgas     = flow->rgas;

     int32_t i;
     double  dp_rms, resp_rms;
     bool    debug = false;

/*   cellwise continuity errors
     -------------------------- */
     resp_rms = 0.0;
     for (i=0; i<nx-1; i++) {
       bb[i] = -(a[i+1]*rho[i+1]*u[i+1] - a[i]*rho[i]*u[i]);
       resp_rms += bb[i]*bb[i];
     }
     bb[nx-1] = 0.0;
     history->resp_rms = sqrt(resp_rms/(nx-1));

/*   set pressure correction coefficients
     ------------------------------------ */
     for (i=0; i<nx-1; i++) {
       if (i > 0) {
         aw[i] = -a[i]*rho[i]/apu[i];
         ap[i] = -aw[i];
       }
       else {
         aw[0] =  0.0;
         if (setup->compress) ap[0] = a[0]*cp / (rgas*u[0]);      // linearised inlet bc
         else                 ap[0] = a[0] /u[0];
       }
       ae[i] = -a[i+1]*rho[i+1]/apu[i+1];
       ap[i] =  ap[i] - ae[i];
     }

/*   compressibility corrections
     --------------------------- */
     if (setup->compress) {
       for (i=0; i<nx-1; i++) {
         if (i > 0) {
           aw[i] = aw[i] - a[i]*u[i]/(rgas*t[i]);
         }
         else {
            ap[i] = ap[i] - a[i]*u[i]/(rgas*t[i]);
         }
         ap[i] = ap[i] + a[i+1]*u[i+1]/(rgas*t[i+1]);
       }
     }

/*   debug information
     ----------------- */
     if (debug){
       for (i=0; i<nx-1; i++)
         printf("Station = %d, \taap = % f \t% f \t% f \tbb = % f\r\n",
                           i,aw[i],ap[i],ae[i],bb[i]);
     }


     return 0;
}

/******************************************************************************
 *  simple_psolve - solve pressure corrrection equation                       *
 ******************************************************************************/
/**
 * \brief   Solve pressure corrrection equation
 * \details Solve pressure corrrection equation using TDMA algorithm,
 *          discretisation using colocated grid with opposed
 *          differencing.
 *
 * \param[in, out]
 *          nozzle  \ref NOZZLE structure to hold nozzle parameters.
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    8th Dec 2022
 */

int  simple_psolve(NOZZLE *nozzle)
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
     double  *t       = arrays->t;
     double  *a       = arrays->a;
     double  *apu     = arrays->apu;

     double  gamma    = flow->gamma;
     double  cp       = flow->cp;
     double  rgas     = flow->rgas;
     double  relaxp   = solver->relaxp;

     int32_t i;
     double  *bb, *aw, *ap, *ae, *dp;
     double  dp_rms, resp_rms;
     bool    debug = false;

/*   alloc local arrays
     ------------------ */
     bb = (double*) malloc(nx*sizeof(double));
     aw = (double*) malloc(nx*sizeof(double));
     ap = (double*) malloc(nx*sizeof(double));
     ae = (double*) malloc(nx*sizeof(double));
     dp = (double*) malloc(nx*sizeof(double));

/*   get pc coefficients and continuity residual
     ------------------------------------------- */
     simple_pcoeff(nozzle, aw, ap, ae, bb);

/*   solve pressure correction equations using TDMA algorithm
     -------------------------------------------------------- */
     ae[0] = ae[0]/ap[0];
     bb[0] = bb[0]/ap[0];
     for (i=1; i<nx-1; i++)
     {
       ap[i] = ap[i] - aw[i]*ae[i-1];
       ae[i] = ae[i]/ap[i];
       bb[i] = (bb[i] - aw[i]*bb[i-1])/ap[i];
     }

     dp[nx-1] = 0;
     for (i=nx-2; i>=0; i--)
     {
       dp[i] = bb[i] - ae[i]*dp[i+1];
     }

/*   Update flow u & p
     ----------------- */
     dp_rms = 0.0;
     for (i=0; i<nx; i++)
     {
       p[i] = p[i] + dp[i];
       if (i > 0) u[i] = u[i] - (dp[i] - dp[i-1])/apu[i];
       dp_rms += dp[i]*dp[i];
     }
     history->dp_rms = sqrt(dp_rms/(nx-1));

/*   update inlet velocity
     --------------------- */
     noz_set_bcs(nozzle);

/*   debug information
     ----------------- */
     if (debug){
       for (i=0; i<nx-1; i++) printf("Station = %d, \tdp = % f\r\n", i,dp[i]);
     }

/*   free local arrays
     ----------------- */
     free(bb);
     free(aw);
     free(ap);
     free(ae);
     free(dp);

     return 0;
}

/******************************************************************************
 *  simple_run - run pressure corrrection solver                              *
 ******************************************************************************/
/**
 * \brief   Apply SIMPLE solver
 * \details Self contained solver that can be called with just a \ref NOZZLE
 *          instance
 *
 * \param[in]
 *          nozzle    \ref NOZZLE instance
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    8th Dec 2022
 */

int  simple_run(NOZZLE *nozzle)
{
     SETUP   *setup   = nozzle->setup;
     SOLVER  *solver  = nozzle->solver;
     HISTORY *history = nozzle->history;
     FILES   *files   = nozzle->files;

     int32_t niters = solver->niters;
     int32_t iter;

/*   check setup
     ----------- */
     if(setup->couple){
       error_message("simple_run", __FILE__, __LINE__, "%s", "solver not set for SIMPLE");
       return 1;
     }

/*   compressible flow
     ----------------- */
     if(setup->compress){
       for(iter=0; iter<niters; iter++){
         noz_set_bcs(nozzle);
         noz_temp_mach(nozzle);
         noz_density(nozzle);
         simple_msolve(nozzle);
         simple_psolve(nozzle);

         hist_print(nozzle);
         if(hist_check_conv(nozzle)) break;
         history->iter++;
       }
     }

/*   incompressible flow
     ------------------- */
     else{
       for(iter=0; iter<niters; iter++){
         noz_set_bcs(nozzle);
         simple_msolve(nozzle);
         simple_psolve(nozzle);

         hist_print(nozzle);
         if(hist_check_conv(nozzle)) break;
         history->iter++;
       }
     }

/*   save solution using casename
     ---------------------------- */
     print_flow(nozzle, files->solnfile);

     return 0;
}

/******************************************************************************
 *  simple_pcmat - return pressure corrrection matrix                         *
 ******************************************************************************/
/**
 * \brief   Get sparse matrix for pressure corrrection equation
 * \details Coefficients of PC equation are returned in \ref SMATRIX structure.
 *          Calling routine should pass pointer to an empty sparse matrix
 *          initialised using \ref NEW_SMAT.
 *          The continuity residuals are also returned within the structure.
 * \n       This routine performs the first half of the SIMPLE iteration.
 *
 * \param[in, out]
 *          nozzle  \ref NOZZLE structure to hold nozzle parameters.
 * \param[in, out]
 *          S       \ref SMATRIX structure to hold pc equations.
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    20th Dec 2022
 */

int  simple_pcmat(NOZZLE *nozzle, SMATRIX *S)
{
     SETUP  *setup  = nozzle->setup;
     ARRAYS *arrays = nozzle->arrays;
     int32_t nx     = arrays->nx;
     double  *bb, *aw, *ap, *ae;
     int32_t nnz, iz, i;
     int32_t *col, *rst;
     double  *val;

/*   perform first part of SIMPLE iteration
     -------------------------------------- */
     noz_set_bcs(nozzle);
     if(setup->compress){
       noz_temp_mach(nozzle);
       noz_density(nozzle);
       simple_msolve(nozzle);
     }
     else{
       simple_msolve(nozzle);
     }

/*   alloc local arrays
     ------------------ */
     bb = (double*) malloc(nx*sizeof(double));
     aw = (double*) malloc(nx*sizeof(double));
     ap = (double*) malloc(nx*sizeof(double));
     ae = (double*) malloc(nx*sizeof(double));

/*   get pc coefficients and continuity residual
     ------------------------------------------- */
     simple_pcoeff(nozzle, aw, ap, ae, bb);

/*   free incoming matrix data
     ------------------------- */
     sparse_free(S);

/*   allocate memory for sparse matrix
     --------------------------------- */
     nnz = 3*(nx - 1);
     col = (int32_t*) malloc(nnz    * sizeof(int32_t));
     rst = (int32_t*) malloc((nx+1) * sizeof(int32_t));
     val = (double*)  malloc(nnz    * sizeof(double));
     
/*   set sparse matrix
     ----------------- */
     S->nr  = nx;
     S->nc  = nx;
     S->nnz = nnz;
     S->col = col;
     S->rst = rst;
     S->val = val;
     S->b   = bb;

/*   populate matrix
     --------------- */
     iz = 0;
     rst[0]  = iz;
     col[iz] = 0;
     val[iz] = ap[0];
     iz ++;
     col[iz] = 1;
     val[iz] = ae[0];
     iz ++;

     for (i=1; i<nx-1; i++){
       rst[i]  = iz;
       col[iz] = i-1;
       val[iz] = aw[i];

       iz ++;
       col[iz] = i;
       val[iz] = ap[i];
       iz ++;

       col[iz] = i+1;
       val[iz] = ae[i];
       iz ++;
     }

     rst[nx-1] = iz;
     col[iz]   = nx-1;
     val[iz]   = 1.0;
     iz ++;
     rst[nx] = iz;

/*   free local arrays
     ----------------- */
     free(aw);
     free(ap);
     free(ae);

/*   check
     ----- */
     if(iz != nnz){
       error_message("simple_pcmat",  __FILE__, __LINE__,
                     "mis-match in non-zeros: iz = %d, nnz = %d", iz, nnz);
       return 1;
     }

     return 0;
}

/******************************************************************************
 *  simple_pcupd - update flow field from external solution to pc matrix      *
 ******************************************************************************/
/**
 * \brief   Update flow field from external solution to pressure correction matrix
 * \details Updates pressure and velocity field from pressure corrections.
 *
 * \param[in, out]
 *          nozzle  \ref NOZZLE structure to hold nozzle parameters.
 * \param[in, out]
 *          dp      solution to pressure correction equation
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    24th Dec 2022
 */

int  simple_pcupd(NOZZLE *nozzle, double *dp)
{
     SETUP  *setup    = nozzle->setup;
     ARRAYS *arrays   = nozzle->arrays;
     HISTORY *history = nozzle->history;

     int32_t nx   = arrays->nx;
     double  *u   = arrays->u;
     double  *p   = arrays->p;
     double  *apu = arrays->apu;

     int32_t i;
     double  dp_rms, resp_rms;
     bool    debug = false;

/*   Update flow u & p
     ----------------- */
     dp_rms = 0.0;
     for (i=0; i<nx; i++)
     {
       p[i] = p[i] + dp[i];
       if (i > 0) u[i] = u[i] - (dp[i] - dp[i-1])/apu[i];
       dp_rms += dp[i]*dp[i];
     }
     history->dp_rms = sqrt(dp_rms/(nx-1));

/*   update inlet velocity
     --------------------- */
     noz_set_bcs(nozzle);

/*   debug information
     ----------------- */
     if (debug){
       for (i=0; i<nx-1; i++) printf("Station = %d, \tdp = % le\r\n", i,dp[i]);
     }

     return 0;
}

