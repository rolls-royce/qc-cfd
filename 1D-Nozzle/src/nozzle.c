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
 * \file    nozzle.c
 * \brief   Routines for handing \ref NOZZLE structure
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    1st Dec 2022
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
#include <unistd.h>       // gives access command to check file exists

#include <ezxml.h>
#include <proto.h>

/******************************************************************************
 *  noz_create - create an empty nozzle structure                             *
 ******************************************************************************/
/**
 * \brief   Create an empty nozzle structure
 * \details Initialises all values to default values and all pointers to NULL
 *
 * \param[in, out]
 *          nozzle  \ref NOZZLE structure to be created
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    1st Dec 2022
 */

int  noz_create(NOZZLE *nozzle)
{
     SETUP   *setup;
     MESH    *mesh;
     FLOW    *flow;
     SOLVER  *solver;
     FILES   *files;
     ARRAYS  *arrays;
     HISTORY *history;

/*   allocate structures
     ------------------- */
     nozzle->setup   = (SETUP*)   malloc(sizeof(SETUP));
     nozzle->mesh    = (MESH*)    malloc(sizeof(MESH));
     nozzle->flow    = (FLOW*)    malloc(sizeof(FLOW));
     nozzle->solver  = (SOLVER*)  malloc(sizeof(SOLVER));
     nozzle->files   = (FILES*)   malloc(sizeof(FILES));
     nozzle->arrays  = (ARRAYS*)  malloc(sizeof(ARRAYS));
     nozzle->history = (HISTORY*) malloc(sizeof(HISTORY));

     setup   = nozzle->setup;
     mesh    = nozzle->mesh;
     flow    = nozzle->flow;
     solver  = nozzle->solver;
     files   = nozzle->files;
     arrays  = nozzle->arrays;
     history = nozzle->history;

/*   initialise setup
     ---------------- */
     setup->casename = NULL;
     setup->compress = false;
     setup->couple   = false;

/*   initialise mesh
     --------------- */
     arrays->nx =  0;
     mesh->ain  =  1.0;
     mesh->da   = -0.1;

/*   initialise flow
     --------------- */
     flow->rho    = 1.0;
     flow->ptin   = 1.5;
     flow->ttin   = 1.0;
     flow->min    = 0.8;
     flow->pex    = 0.95;
     flow->gamma  = 1.4;
     flow->cp     = 3.5;
     flow->rgas   = 1.0;
     flow->upwind = 2;

/*   initialise solver
     ----------------- */
     solver->niters = 0;
     solver->relaxu = 1.0;
     solver->relaxp = 1.0;
     solver->relaxc = 1.0;
     solver->presim = 0;
     solver->rrampc = 0;
     solver->toleru = 1.0e-12;
     solver->tolerp = 1.0e-12;
     solver->tolerc = 1.0e-12;

/*   initialise files
     ---------------- */
     files->histfile[0] = '\0';
     files->solnfile[0] = '\0';

/*   initialise arrays
     ----------------- */
     arrays->a    = NULL;
     arrays->rho  = NULL;
     arrays->u    = NULL;
     arrays->p    = NULL;
     arrays->t    = NULL;
     arrays->M    = NULL;
     arrays->apu  = NULL;
     arrays->peff = NULL;
     arrays->a0   = NULL;
     arrays->a1   = NULL;
     arrays->a2   = NULL;

/*   initialise history
     ------------------ */
     history->iter     = 0;
     history->dr_rms   = 0.0;
     history->du_rms   = 0.0;
     history->dp_rms   = 0.0;
     history->resr_rms = 0.0;
     history->resu_rms = 0.0;
     history->resp_rms = 0.0;

     return 0;
}

/******************************************************************************
 *  noz_free - free nozzle structure                                          *
 ******************************************************************************/
/**
 * \brief   Free nozzle structure
 * \details Frees all dynamically allocated memory associated with the structure
 *
 * \param[in, out]
 *          nozzle  \ref NOZZLE structure to be freed
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    2nd Dec 2022
 */

void noz_free(NOZZLE *nozzle)
{
     SETUP  *setup  = nozzle->setup;
     ARRAYS *arrays = nozzle->arrays;

/*   free memory within structures
     ----------------------------- */
     if(setup->casename) free(setup->casename);
     setup->casename = NULL;

     if(arrays->a)    free(arrays->a);
     if(arrays->rho)  free(arrays->rho);
     if(arrays->u)    free(arrays->u);
     if(arrays->p)    free(arrays->p);
     if(arrays->t)    free(arrays->t);
     if(arrays->M)    free(arrays->M);
     if(arrays->apu)  free(arrays->apu);
     if(arrays->peff) free(arrays->peff);
     if(arrays->a0)   free(arrays->a0);
     if(arrays->a1)   free(arrays->a1);
     if(arrays->a2)   free(arrays->a2);

     arrays->a    = NULL;
     arrays->rho  = NULL;
     arrays->u    = NULL;
     arrays->p    = NULL;
     arrays->t    = NULL;
     arrays->M    = NULL;
     arrays->apu  = NULL;
     arrays->peff = NULL;
     arrays->a0   = NULL;
     arrays->a1   = NULL;
     arrays->a2   = NULL;

/*   free structures
     --------------- */
     free(nozzle->setup);
     free(nozzle->mesh);
     free(nozzle->flow);
     free(nozzle->solver);
     free(nozzle->files);
     free(nozzle->arrays);
     free(nozzle->history);
}

/******************************************************************************
 *  noz_files - set file names using casename and options                     *
 ******************************************************************************/
/**
 * \brief   Set file names using casename and options
 * \details File names have the form {casename}_{nx}_{xx}{yy}.{type}.
 *          Where {nx} = the number of axial stations using 3 digits;
 *                {xx} = co or in to denote compressible or incompressible flow;
 *                {yy} = co or si to denote coupled or SIMPLE solver; and, 
 *                {type} identifies the data stored on the file.
 *          For example, "nozzle_032_insi.hist" is the convergence history file
 *          for casename "nozzle" with incompressible flow using the SIMPLE
 *          solver on a mesh with 32 stations.
 *
 * \param[in]
 *          nozzle  \ref NOZZLE structure to be printed
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    7th Dec 2022
 */

void noz_files(NOZZLE *nozzle)
{
     SETUP  *setup  = nozzle->setup;
     ARRAYS *arrays = nozzle->arrays;
     FILES  *files  = nozzle->files;
     char   *cname  = setup->casename;
     int32_t nx     = arrays->nx;
     char   *comp, *coup;

/*   set case dependent suffices
     --------------------------- */
     if(setup->compress) comp = strdup("co");
     else                comp = strdup("in");

     if(setup->couple)   coup = strdup("co");
     else                coup = strdup("si");

/*   set file names
     -------------- */
     snprintf(files->histfile, NAME_MAX, "%s_%03d_%s%s.hst", cname, nx, comp, coup);
     snprintf(files->solnfile, NAME_MAX, "%s_%03d_%s%s.sol", cname, nx, comp, coup);

/*   tidy up
     ------- */
     free(comp);
     free(coup);
}

/******************************************************************************
 *  noz_initialise - initialise nozzle mesh and flow field                    *
 ******************************************************************************/
/**
 * \brief   Initialise nozzle mesh and flow field
 * \details Initialises either con-di nozzle for compressible flow or
 *          a constant curvature channel for incompressible flow.
 *          Initial flow field is set using either isentropic relations for
 *          compressible flows or Bernoulli's equation for incompessible flows.
 *
 * \param[in,out]
 *          nozzle \ref NOZZLE structure to contain nozzle data
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    7th Dec 2022
 */

int  noz_initialise(NOZZLE *nozzle)
{
     SETUP  *setup  = nozzle->setup;
     MESH   *mesh   = nozzle->mesh;
     FLOW   *flow   = nozzle->flow;
     SOLVER *solver = nozzle->solver;
     ARRAYS *arrays = nozzle->arrays;

     int32_t nx, i;
     double  ain, da, rin, ptin, ttin, Min, dM, pex, gamma, cp, rgas;
     double  *a, *rho, *u, *p, *t, *M, *apu;
     double  f, Mi, ui, pi, ti, ue, pe, te, power;
     double  x0, y0, r, x;

/*   allocate memory
     --------------- */
     nx = arrays->nx;
     arrays->a    = (double*) malloc(nx*sizeof(double));
     arrays->rho  = (double*) malloc(nx*sizeof(double));
     arrays->u    = (double*) malloc(nx*sizeof(double));
     arrays->p    = (double*) malloc(nx*sizeof(double));
     arrays->t    = (double*) malloc(nx*sizeof(double));
     arrays->M    = (double*) malloc(nx*sizeof(double));
     arrays->apu  = (double*) malloc(nx*sizeof(double));
     if(setup->compress){
       arrays->peff = (double*) malloc(nx*sizeof(double));
       arrays->a0   = (double*) malloc(nx*sizeof(double));
       arrays->a1   = (double*) malloc(nx*sizeof(double));
       arrays->a2   = (double*) malloc(nx*sizeof(double));
     }

     a   = arrays->a;
     rho = arrays->rho;
     u   = arrays->u;
     p   = arrays->p;
     t   = arrays->t;
     M   = arrays->M;

/*   get nozzle parameters
     --------------------- */
     ain  = mesh->ain;
     da   = mesh->da;
     rin  = flow->rho;
     ptin = flow->ptin;
     ttin = flow->ttin;
     pex  = flow->pex;

/*   generate compressible mesh and flow field
     ----------------------------------------- */
     if(setup->compress){
       Min   = flow->min;
       pex   = flow->pex;
       gamma = flow->gamma;
       cp    = flow->cp;
       rgas  = flow->rgas;
       dM    = 1.0/(double) (nx-1);
       power = 0.5*(gamma + 1.0)/(gamma - 1.0);

/*     set nozzle areas
       ---------------- */
       for (i=0; i<nx; i++){
         Mi   = Min + i*dM;
         a[i] = (2 + (gamma-1)*Mi*Mi) / (gamma+1);
         a[i] = pow(a[i], power) / Mi;
         a[i] = ain*a[i];                      // use ain to scale area and hence mat coeffs
       }

/*     set inlet and exit values
       ------------------------- */
       pi = pow(1.0 + 0.5*(gamma-1)*Min*Min, -gamma/(gamma-1));
       ti = pow(pi,(gamma-1)/gamma);
       ui = sqrt(2.0*cp*(ttin - ti));

       pe = pex;
       te = pow(pe,(gamma-1)/gamma);
       ue = sqrt(0.5*cp*(ttin - te));

/*     interpolate flow variables
       --------------------------*/
       for (i=0; i<nx; i++){
         f      = (float) i / (float) (nx-1);
         p[i]   = pi + (pe - pi)*f;
         u[i]   = ui + (ue - ui)*f;
         t[i]   = ttin - 0.5*u[i]*u[i]/cp;
         M[i]   = u[i]/sqrt(gamma*rgas*t[i]);
         rho[i] = p[i]/(rgas*t[i]);
       }
     }

/*   generate incompressible mesh and flow field
     ------------------------------------------- */
     else{

/*     solve for circle passing through (0, ain), (0.5, f*ain), (1, ain)
       ----------------------------------------------------------------- */
       x0 = 0.5;
       if (da != 0){
         f  = 1.0 + da;
         y0 = (f + 1.0)*ain - 0.25/(ain*(f - 1.0));
         y0 = 0.5*y0;
         r  = MAX(f*ain - y0, y0 - f*ain);
         if(da > 0) f =  1.0;
         else       f = -1.0;
         printf("circle has origin at (%le, %le) and radius %le\n", x0, y0, r);
       }
       else{
         y0 = ain;
         r  = 1.0;
         f  = 0.0;
       }
/*     set areas and flow properties
       ----------------------------- */
       for (i=0; i<nx; i++){
         x      = (double) i / (double) (nx-1);
         a[i]   = y0 + f*sqrt(r*r -((x-x0)*(x-x0)));
         rho[i] = rin;
         p[i]   = pex;
         u[i]   = sqrt(MAX(0.0, 2*(ptin-p[i])/rho[i]));
         t[i]   = 1.0;
         M[i]   = 0.0;
       }
     }

     return 0;
}

/******************************************************************************
 *  noz_continuity - adjust velocity to satisfy continuity                    *
 ******************************************************************************/
/**
 * \brief   Adjust velocity to satisfy continuity
 * \details Update velocities to match inlet flow rate
 *
 * \param[in, out]
 *          nozzle  \ref NOZZLE structure to hold nozzle parameters.
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    7th Dec 2022
 */

int  noz_continuity(NOZZLE *nozzle)
{
     ARRAYS *arrays = nozzle->arrays;
     int32_t nx   = arrays->nx;
     double  *a   = arrays->a;
     double  *rho = arrays->rho;
     double  *u   = arrays->u;
     double  flowrate;
     int32_t i;

/*   set u to match inlet flow rate
     ------------------------------ */
     flowrate = a[0] * rho[0] * u[0];
     for (i=1; i<nx; i++)
       u[i] = flowrate / (a[i] * rho[i]);

     return 0;
}

/******************************************************************************
 *  noz_set_bcs - adjust velocity to satisfy continuity                       *
 ******************************************************************************/
/**
 * \brief   Set inlet velocity and temperature
 * \details Uses isentropic relations for compressible flow and Bernoulli's
 *          equation for incompressible flow.
 *
 * \param[in, out]
 *          nozzle  \ref NOZZLE structure to hold nozzle parameters.
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    7th Dec 2022
 */

int  noz_set_bcs(NOZZLE *nozzle)
{
     SETUP   *setup  = nozzle->setup;
     SOLVER  *solver = nozzle->solver;
     FLOW    *flow   = nozzle->flow;
     ARRAYS  *arrays = nozzle->arrays;

     int32_t nx    = arrays->nx;
     double  *rho  = arrays->rho;
     double  *u    = arrays->u;
     double  *p    = arrays->p;
     double  *t    = arrays->t;

     double  gamma = flow->gamma;
     double  cp    = flow->cp;
     double  rgas  = flow->rgas;
     double  ptin  = flow->ptin;
     double  ttin  = flow->ttin;
     double  pex   = flow->pex;
     double  u2;

/*   exit pressure
     ------------- */
     p[nx-1] = pex;

/*   inlet - compressible
     -------------------- */
     if (setup->compress){
       t[0] = pow(p[0],(gamma-1)/gamma);
       u[0] = pow(2.0*cp*(ttin - t[0]),0.5);
     }

/*   inlet - incompressible
     ---------------------- */
     else{
       u2   = 2.0*(ptin - p[0])/rho[0];
       u[0] = sqrt(MAX(0.0, u2));
       t[0] = p[0]/(rgas * rho[0]);
     }

     return 0;
}

/******************************************************************************
 *  noz_temp_mach - set temperature and Mach number                           *
 ******************************************************************************/
/**
 * \brief   Set temperature and Mach number
 * \details Set temperature and Mach number from velocity and reference total
 *          temperature (=1 in non-dimensional units)
 *
 * \param[in, out]
 *          nozzle  \ref NOZZLE structure to hold nozzle parameters.
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    7th Dec 2022
 */

int  noz_temp_mach(NOZZLE *nozzle)
{
     SETUP   *setup  = nozzle->setup;
     FLOW    *flow   = nozzle->flow;
     ARRAYS  *arrays = nozzle->arrays;

     int32_t nx    = arrays->nx;
     double  *u    = arrays->u;
     double  *t    = arrays->t;
     double  *M    = arrays->M;

     double  gamma = flow->gamma;
     double  cp    = flow->cp;
     double  rgas  = flow->rgas;
     double  ttin  = flow->ttin;
     int32_t i;

/*   set temperature and Mach number
     ------------------------------- */
     if(setup->compress){
       for (i=0; i<nx; i++) {
         t[i] = MAX(0.01, ttin - 0.5*u[i]*u[i]/cp);
         M[i] = u[i]/sqrt(gamma*rgas*t[i]);
       }
     }
     else{
       for (i=0; i<nx; i++) {
         t[i] = 1.0;
         M[i] = 0.0;
       }
     }

     return 0;
}

/******************************************************************************
 *  noz_peff - set effective pressure for upwinded density                    *
 ******************************************************************************/
/**
 * \brief   Set effective pressure for upwinded density
 * \details Upwind density using one of the Moore and Moore schemes.
 *          Also stores upwinding coefficients for implicit coupled scheme.
 *
 * \param[in, out]
 *          nozzle  \ref NOZZLE structure to hold nozzle parameters.
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    7th Dec 2022
 */
int  noz_peff(NOZZLE *nozzle)
{
     SETUP  *setup  = nozzle->setup;
     FLOW   *flow   = nozzle->flow;
     ARRAYS *arrays = nozzle->arrays;

     int32_t nx     = arrays->nx;
     double  *p     = arrays->p;
     double  *M     = arrays->M;
     double  *peff  = arrays->peff;
     double  *a0up  = arrays->a0;
     double  *a1up  = arrays->a1;
     double  *a2up  = arrays->a2;
     int     upwind = flow->upwind;

     int32_t i;
     double  MaxM, a0, a1, a2, c1, c2, k1, k2, m1, m2;

/*   check for compressibility
     ------------------------- */
     if (!setup->compress){
       for (i=0; i<nx; i++) peff[i] = p[i];
       return 0;
     }

/*   inlet density from gas law
     -------------------------- */
     peff[0] = p[0];
     a0up[0] = 0.0;
     a1up[0] = 0.0;
     a2up[0] = 0.0;

/*   apply upwinding
     --------------- */
     switch(upwind)
     {

/*     1st order upwinding
       ------------------- */
       case RHO_UPWIND_1ST:
            for (i=1; i<nx; i++)
            {
              MaxM = MAX(M[i], M[i-1]);
              if (MaxM < 2.0)
                 a0 = MIN(1.0, (0.8/3)*(4.0/(MaxM*MaxM) - 1));
              else
                 a0 = 0.0;

              peff[i] = p[i-1] + a0*(p[i] - p[i-1]);
              peff[i] = MAX(peff[i], MIN(p[i], p[i-1]));
              peff[i] = MIN(peff[i], MAX(p[i], p[i-1]));

              a0up[i] = a0;
              a1up[i] = 0.0;
            }
            break;

/*     2nd order upwinding
       ------------------- */
       case RHO_UPWIND_2ND:
            for (i=1; i<nx; i++)
            {
              MaxM = MAX(M[i], M[i-1]);
              if (MaxM < 2.0)
              {
                a0 = MIN(1.0, (0.8/3)*(4.0/(MaxM*MaxM) - 1));
                a1 = 1.0 - a0;
              }
              else
              {
                a0 = 0.0;
                a1 = 4.0/(MaxM*MaxM);
              }
              peff[i] = p[i-1] + a0*(p[i] - p[i-1]);

              if (i>1)
              {
                c1 = (p[i] - p[i-1])*(p[i] - p[i-2]);
                if (c1 > 0.0) peff[i] = peff[i] + a1*(p[i] - p[i-2])/2.0;
                else  a1 = 0;
              }
              else a1 = 0;

              if (i>2)
              {
                a2 = 1.0 - a0 - a1;
                c2 = (p[i] - p[i-1])*(p[i] - p[i-3]);
                if (c2 > 0.0) peff[i] = peff[i] + a2*(p[i] - p[i-3])/3.0;
                else  a2 = 0;
              }
              else a2 = 0;

              peff[i] = MAX(peff[i], MIN(p[i], p[i-1]));
              peff[i] = MIN(peff[i], MAX(p[i], p[i-1]));

              a0up[i] = a0;
              a1up[i] = a1;
              a2up[i] = a2;
            }
            break;

/*     3 point scheme
       -------------- */
       case RHO_UPWIND_3PT:
            a0 = 0.0;
            a2 = 1.0;
            for (i=1; i<nx; i++)
            {
              peff[i] = p[i-1] + a0*(p[i] - p[i-1]);
              if (i>2)
              {
                peff[i] = peff[i] + a2*(p[i] - p[i-3])/3.0;
              }
              peff[i] = MAX(peff[i], MIN(p[i], p[i-1]));
              peff[i] = MIN(peff[i], MAX(p[i], p[i-1]));

              a0up[i] = a0;
              a1up[i] = 0.0;
              a2up[i] = a2;
            }
            break;

/*      McGuirk and Page scheme (values hardwired from AIAA paper)
        ---------------------------------------------------------- */
        case RHO_UPWIND_JJG:
            k1 = 1.0;
            m1 = 1.0;
            m2 = 0.7;

            for (i=1; i<nx; i++)
            {
              k2 = 0.7;
              if(p[i] > p[i-1]) k2 = 0.5;      // setting k2=0 as in <&P paper causes overshoot

              MaxM = M[i];
//            MaxM = 0.5*(M[i] + M[i-1]);
              a0 = MAX(0.0, MAX(k1*(1.0 - (m1*m1/(MaxM*MaxM))), 
                                k2*(1.0 - (m2*m2/(MaxM*MaxM)))));

              peff[i] = p[i] - a0*(p[i] - p[i-1]);

              a0up[i] = a0;
              a1up[i] = 0.0;
              a2up[i] = 0.0;
            }
            break;

/*      invalid scheme
        -------------- */
        default:
            error_message("noz_peff",  __FILE__, __LINE__, "invalid density upwinding scheme - %d", upwind);
            return 1;
     }

     return 0;
}

/******************************************************************************
 *  noz_density - set upwinded density                                        *
 ******************************************************************************/
/**
 * \brief   Set upwinded density
 * \details Calls \ref noz_peff to set effective pressure
 *
 * \param[in, out]
 *          nozzle  \ref NOZZLE structure to hold nozzle parameters.
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    7th Dec 2022
 */

int  noz_density(NOZZLE *nozzle)
{
     SETUP  *setup  = nozzle->setup;
     FLOW   *flow   = nozzle->flow;
     ARRAYS *arrays = nozzle->arrays;

     int32_t nx     = arrays->nx;
     double  *rho   = arrays->rho;
     double  *t     = arrays->t;
     double  *peff  = arrays->peff;
     double  rgas   = flow->rgas;  

     int32_t i, rc;

/*   check for compressibility
     ------------------------- */
     if (!setup->compress) return 0;

/*   set density from effective pressure
     ----------------------------------- */
     if((rc = noz_peff(nozzle))){
       error_traceback("noz_density", __FILE__, __LINE__);
       return rc;
     }

     for (i=0; i<nx; i++) rho[i] = peff[i]/(rgas*t[i]);

     return 0;
}

