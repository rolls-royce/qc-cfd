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
 * \file    py_interface.c
 * \brief   Interface for calling C routines from Python
 * \details Interface uses Python ctypes
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    12th Dec 2022
 */

/******************************************************************************
 *  Include files                                                             *
 ******************************************************************************/
#include <stdio.h>
#include <stdlib.h>

#include <nozzle.h> 
#include <proto.h> 
#include <c2py.h> 

/******************************************************************************
 *  c2py_noz_create - create an empty nozzle structure                        *
 ******************************************************************************/
/**
 * \brief   Create an empty nozzle structure
 * \details Wrapper for \ref noz_create.
 *          Returns a (void*) pointer to the calling python routine.
 *
 * \param[in, out]
 *          py_noz  \ref NOZZLE structure to hold nozzle parameters.
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    12th Dec 2022
 */

int  c2py_noz_create(void **py_noz)
{
     NOZZLE *nozzle;
     int    rc;

/*   create nozzle structure
     ----------------------- */
     nozzle = (NOZZLE*) malloc(sizeof(NOZZLE));
     if((rc = noz_create(nozzle)))
     {
       error_traceback("noz_density", __FILE__, __LINE__);
       return rc;
     }
     printf("c2py_noz_create: created empty nozzle at address = %p\n", nozzle);

/*   set return pointer
     ------------------ */
     *py_noz = (void*) nozzle;

     return rc;
}

/******************************************************************************
 *  c2py_noz_free - ree nozzle structure                                     *
 ******************************************************************************/
/**
 * \brief   Free nozzle structure
 * \details Wrapper for \ref noz_free. Frees allocated memory within \ref NOZZLE
 *          structure but not the structure itself. The latter is the
 *          responsibility of the calling routine.
 *
 * \param[in]
 *          py_noz  \ref NOZZLE structure to hold nozzle parameters.
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    27th Dec 2022
 */

void c2py_noz_free(void *py_noz)
{
     NOZZLE *nozzle = (NOZZLE*) py_noz;

/*   free nozzle
     ----------- */
     noz_free(nozzle);
}

/******************************************************************************
 *  c2py_read_input - read XML input file                                     *
 ******************************************************************************/
/**
 * \brief   Read XML input file for 1D nozzle
 * \details Wrapper for \ref read_input_xml.
 *          Call \ref c2py_noz_create prior to calling this routine in order to
 *          instantiate an empty \ref NOZZLE structure.
 *
 * \param[in]
 *          inputfile  name XML file containg the input data
 * \param[in, out]
 *          py_noz  \ref NOZZLE structure to hold nozzle parameters.
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    19th Dec 2022
 */

int  c2py_read_input(const char *inputfile, void *py_noz)
{
     NOZZLE *nozzle = (NOZZLE*) py_noz;
     int    rc;

/*   read input file
     --------------- */
     rc = read_input_xml(inputfile, nozzle);
     if(rc){
       error_traceback("c2py_read_input", __FILE__, __LINE__);
     }

     return rc;
}

/******************************************************************************
 *  c2py_noz_init - initialise nozzle                                         *
 ******************************************************************************/
/**
 * \brief   Initialise nozzle
 * \details Wrapper for \ref noz_initialise \ref couple_run.
 *          Call \ref c2py_noz_create and \ref c2py_read_input prior to this
 *          routine to instantiate nozzle and get input parameters.
 * \n       Any pre-SIMPLE iterations required by the coupled solver are
 *          performed here and viewed as part of the initialisation. This is
 *          different from the C version where the pre-SIMPLE iterations are 
 *          included in \ref couple_run.
 *
 * \param[in, out]
 *          py_noz  \ref NOZZLE structure to hold nozzle parameters.
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    19th Dec 2022
 */

int  c2py_noz_init(void *py_noz)
{
     NOZZLE *nozzle = (NOZZLE*) py_noz;
     SETUP  *setup  = nozzle->setup;
     SOLVER *solver = nozzle->solver;

     int32_t niters  = solver->niters;
     int32_t nsimple = solver->presim;
     int     rc;

/*   initialise nozzle mesh and flow field
     ------------------------------------- */
     print_setup(nozzle);
     if((rc = noz_initialise(nozzle))){
       error_traceback("c2py_noz_init", __FILE__, __LINE__);
       return 1;
     }

/*   perform pre-SIMPLE iterations for coupled solver
     ------------------------------------------------ */
     if(setup->compress && nsimple>0){
       setup->couple  = 0;
       solver->niters = nsimple;

       simple_run(nozzle);

       setup->couple  = 1;
       solver->niters = niters;
     }

     return rc;
}

/******************************************************************************
 *  c2py_run_c - run entirely using C library code                            *
 ******************************************************************************/
/**
 * \brief   Run entirely using C library code 
 * \details Wrapper for \ref simple_run and \ref couple_run
 *
 * \param[in, out]
 *          py_noz  \ref NOZZLE structure to hold nozzle parameters.
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    19th Dec 2022
 */

int  c2py_run_c(void *py_noz)
{
     NOZZLE *nozzle = (NOZZLE*) py_noz;
     SETUP  *setup   = nozzle->setup;
     int    rc;

/*   run solver
     ---------- */
     if (setup->couple){
       if((rc = couple_run(nozzle))){
         error_traceback("main", __FILE__, __LINE__);
         return rc;
       }
     }
     else{
       if((rc = simple_run(nozzle))){
         error_traceback("main", __FILE__, __LINE__);
         return rc;
       }
     }

/*   print solution
     -------------- */
     print_flow(nozzle, NULL);

     return rc;
}

/******************************************************************************
 *  c2py_get_smat - get sparse matrix for pressure correction or coupled eqns *
 ******************************************************************************/
/**
 * \brief   Get sparse matrix for pressure correction or coupled eqns
 * \details Returns sparse matrix components uding compressed row storage
 *          and right hand side residual of the equations.
 * \n       For a coupled discretisation, the complete coupled matrix and
 *          residual are retured.
 * \n       For a SIMPLE discretisation, the momentum equation is solved first
 *          to provide the terms needed for the pressure correction equation. 
 *
 * \param[in]
 *          py_noz  \ref NOZZLE structure to hold nozzle parameters.
 * \param[in,out]
 *          PS      Empty \ref SMATRIX initialised using ctypes SMATRIX class
 *                  from Python.
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    24th Dec 2022
 */

int  c2py_get_smat(void *py_noz, SMATRIX *PS)
{
     NOZZLE  *nozzle = (NOZZLE*) py_noz;
     SETUP   *setup  = nozzle->setup;
     int     rc;

/*   get matrix
     ---------- */
     if (setup->couple){
       if((rc = couple_mat(nozzle, PS))){
         error_traceback("main", __FILE__, __LINE__);
         return rc;
       }
     }
     else{
       if((rc = simple_pcmat(nozzle, PS))){
         error_traceback("main", __FILE__, __LINE__);
         return rc;
       }
     }

     return rc;
}

/******************************************************************************
 *  c2py_put_soln - store python solution and update flow field               *
 ******************************************************************************/
/**
 * \brief   Store python solution and update flow field  
 * \details Updates solution using solution of linear system and completes
 *          interation.
 *
 * \param[in, out]
 *          py_noz  \ref NOZZLE structure to hold nozzle parameters.
 * \param[in]
 *          x       Solution of Ax=b from Python
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    24th Dec 2022
 */

int  c2py_put_soln(void *py_noz, double *x)
{
     NOZZLE  *nozzle = (NOZZLE*) py_noz;
     SETUP   *setup  = nozzle->setup;
     int     rc;

/*   update flow field
     ----------------- */
     if (setup->couple){
       if((rc = couple_update(nozzle, x))){
         error_traceback("main", __FILE__, __LINE__);
         return rc;
       }
     }
     else{  
       if((rc = simple_pcupd(nozzle, x))){
         error_traceback("main", __FILE__, __LINE__);
         return rc;
       }
     }

     return rc;
}

/******************************************************************************
 *  c2py_hist_chk - print current residuals and check convergence             *
 ******************************************************************************/
/**
 * \brief   Print current residuals and check convergence
 * \details Relies on Python ctypes c_bool and stdbool.h
 *          interation. Also checks whether maximum number of iterations
 *          has been reached.
 *
 * \param[in]
 *          py_noz  \ref NOZZLE structure to hold nozzle parameters.
 * \param[out]
 *          status  True if converged, false otherwise
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    24th Dec 2022
 */

int  c2py_hist_chk(void *py_noz, bool *status)
{
     NOZZLE  *nozzle  = (NOZZLE*) py_noz;
     SOLVER  *solver  = nozzle->solver;
     HISTORY *history = nozzle->history;

     int32_t niters   = solver->niters;
     int     rc;

     hist_print(nozzle);
     if(!(*status = hist_check_conv(nozzle))) history->iter++;

     if(history->iter == niters) *status = true;

     return 0;
}

/******************************************************************************
 *  c2py_save_soln - save solution                                            *
 ******************************************************************************/
/**
 * \brief   Save flow field to solution file
 * \details Wrapper for \ref print_flow
 *
 * \param[in]
 *          py_noz  \ref NOZZLE structure to hold nozzle parameters.
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    24th Dec 2022
 */

int  c2py_save_soln(void *py_noz)
{
     NOZZLE *nozzle = (NOZZLE*) py_noz;
     FILES  *files  = nozzle->files;
     int    rc;

/*   read input file
     --------------- */
     rc = print_flow(nozzle, files->solnfile);
     if(rc){
       error_traceback("c2py_read_input", __FILE__, __LINE__);
     }

     return rc;
}

/******************************************************************************
 *  c2py_print_soln - print solution to screen                                *
 ******************************************************************************/
/**
 * \brief   Print flow field to screen
 * \details Wrapper for \ref print_flow
 *
 * \param[in]
 *          py_noz  \ref NOZZLE structure to hold nozzle parameters.
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    24th Dec 2022
 */

int  c2py_print_soln(void *py_noz)
{
     NOZZLE *nozzle = (NOZZLE*) py_noz;
     int    rc;

/*   read input file
     --------------- */
     rc = print_flow(nozzle, NULL);
     if(rc){
       error_traceback("c2py_read_input", __FILE__, __LINE__);
     }

     return rc;
}
