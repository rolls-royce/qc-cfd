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
 * \file    history.c
 * \brief   Routines for convergence histories
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
 *  hist_print - print current iteration                                     *
 ******************************************************************************/
/**
 * \brief   Print convergence data for current iteration
 * \details Wrapper routine for \ref hist_output to write convergence history
 *          to both the screen and the convergence history file specified in
 *          \ref FILES.
 *
 * \param[in]
 *          nozzle   \ref NOZZLE structure containing all nozzle data
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    8th Dec 2022
 */

int  hist_print(NOZZLE *nozzle)
{
     SETUP   *setup   = nozzle->setup;
     FILES   *files   = nozzle->files;
     HISTORY *history = nozzle->history;

     int32_t i = history->iter;
     static  bool first = true;
     int     rc;

     if((rc = hist_output(nozzle, NULL, first))){
       error_traceback("hist_print", __FILE__, __LINE__);
       return rc;
     }

     if((rc = hist_output(nozzle, files->histfile, first))){
       error_traceback("hist_print", __FILE__, __LINE__);
       return rc;
     }

     first = false;
     return 0;
}

/******************************************************************************
 *  hist_output - output convergence history for plotting in gnuplot          *
 ******************************************************************************/
/**
 * \brief   Output convergence history for plotting in gnuplot 
 * \details Pass filename or NULL to write to stdout
 *
 * \param[in]
 *          nozzle   \ref NOZZLE structure needed for flow regime
 * \param[in]
 *          filename name of convergence history file or NULL to write to stdout
 * \param[in]
 *          first    indicator to write first line as header
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    8th Dec 2022
 */

int  hist_output(NOZZLE *nozzle, char *filename, bool first)
{
     SETUP   *setup   = nozzle->setup;
     HISTORY *history = nozzle->history;

     int32_t i = history->iter;
     FILE    *fptr = NULL;

/*   open file if specified
     ---------------------- */
     if(filename){
       if(first) fptr = fopen(filename, "w");
       else      fptr = fopen(filename, "a");

       if(!fptr){
         error_message("hist_output", __FILE__, __LINE__, 
                       "error opening output file: %s\n",filename);
         return 1;
       }
     }
     else
       fptr = stdout;

/*   print file header on first open
     ------------------------------- */
     if(first && fptr){
       if(setup->compress){
         fprintf(fptr, " \"iter\",    \"du_rms\",        \"dp_rms\",        \"dr_rms\",        \"resu_rms\",      \"resp_rms\",      \"resr_rms\"\n");
       }
       else{
         fprintf(fptr, " \"iter\",    \"du_rms\",        \"dp_rms\",        \"resu_rms\",      \"resp_rms\"\n");
       }
     }

/*   print current iteration
     ----------------------- */
     if(setup->compress){
       fprintf(fptr, "%8d, % 12.8le, % 12.8le, % 12.8le, % 12.8le, % 12.8le, % 12.8le\n",
               i+1, history->du_rms,   history->dp_rms, history->dr_rms,
                    history->resu_rms, history->resp_rms, history->resr_rms);

     }
     else{
       fprintf(fptr, "%8d, % 12.8le, % 12.8le, % 12.8le, % 12.8le\n",
               i+1, history->du_rms,   history->dp_rms, history->resu_rms, history->resp_rms);

     }

/*   close file
     ---------- */
     if(filename) fclose(fptr);

     return 0;
}

/******************************************************************************
 *  hist_check_conv - check for convergence                                   *
 ******************************************************************************/
/**
 * \brief   Check for convergence
 * \details Check for convergence
 *
 * \param[in]
 *          nozzle   \ref NOZZLE structure needed for flow regime
 *
 * \return  true/false for convergence achieved
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    8th Dec 2022
 */

bool hist_check_conv(NOZZLE *nozzle)
{
     SETUP   *setup   = nozzle->setup;
     SOLVER  *solver  = nozzle->solver;
     HISTORY *history = nozzle->history;

     double  toleru   = solver->toleru;
     double  tolerp   = solver->tolerp;
     double  tolerc   = solver->tolerc;

     double  dr_rms   = history->dr_rms;
     double  du_rms   = history->du_rms;
     double  dp_rms   = history->dp_rms;

     bool pass = false;

/*   check convergence
     ----------------- */
     if(setup->couple){
        if (du_rms < tolerc && dp_rms < tolerc) pass = true;
     }
     else{
        if (du_rms < toleru && dp_rms < tolerp) pass = true;
     }

     return pass;
}

