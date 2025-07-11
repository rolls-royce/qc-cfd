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
 * \file    io.c
 * \brief   I/O routines for 1-dimensional nozzle
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
#include <unistd.h>       // gives access command to check file exists

#include <ezxml.h>
#include <proto.h>

/******************************************************************************
 *  read_input_xml - read XML input file                                      *
 ******************************************************************************/
/**
 * \brief   Read XML input file for 1D nozzle
 * \details XML file is read using the ezxml library
 *          Only data relevant required by the setup section is read.
 *          The XML structure is designed to allow multiple flow and solver options 
 *          to be stored in a single file to allow testing via a simple change in
 *          the set-up section. This adds to the verbosity of the file but reduces
 *          the number of files to be stored.
 *
 * \param[in]
 *          inputfile  name XML file containg the input data
 * \param[out]
 *          nozzle     \ref NOZZLE structure to contain input data
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    1st Dec 2022
 */

int  read_input_xml(const char *inputfile, NOZZLE *nozzle)
{
     char   routine[] = "read_input_xml";
     SETUP  *setup  = nozzle->setup;
     MESH   *mesh   = nozzle->mesh;
     FLOW   *flow   = nozzle->flow;
     SOLVER *solver = nozzle->solver;
     ARRAYS *arrays = nozzle->arrays;

     ezxml_t xfile, xsetup, xmesh, xinco, xcomp;
     ezxml_t xreg, xbc, xgas, xmeth;
     ezxml_t xdata;
     int     idata;
     char    *xstr;

/*   check input file exists
     ----------------------- */
     if(access(inputfile, F_OK ) != 0){
       error_message(routine, __FILE__, __LINE__, 
                     "input file %s does not exit", inputfile);
       return 1;
     }

/*   load input xml file
     ------------------- */
     printf("reading input from file: %s\n", inputfile);
     xfile = ezxml_parse_file(inputfile);

/*   debug echo file
     --------------- */
     if(false){
       ezxml_pretty(xfile);
       xstr = ezxml_toxml(xfile);
       printf("%s\n",xstr);
     }

/*   load main sections
     ------------------ */
     xsetup  = ezxml_child(xfile, "setup");
     xmesh   = ezxml_child(xfile, "mesh");
     xcomp   = ezxml_child(xfile, "compressible");
     xinco   = ezxml_child(xfile, "incompressible");

/*   read setup data
     --------------- */
     xdata = ezxml_child(xsetup, "casename");
     if(xdata){
       setup->casename = strdup(xdata->txt);
     }
     else{
       error_message(routine, __FILE__, __LINE__, "%s", "casename must be specified");
       return 1;
     }

     if((xdata = ezxml_child(xsetup, "compressible"))){
       sscanf(xdata->txt, "%d", &idata);
       setup->compress = (bool) idata;
     }
     if((xdata = ezxml_child(xsetup, "couple"))){
       sscanf(xdata->txt, "%d", &idata);
       setup->couple   = (bool) idata;
     }

/*   read mesh data
     -------------- */
     xdata = ezxml_child(xmesh, "nx");
     if(xdata) sscanf(xdata->txt, "%d", &arrays->nx);

     xdata = ezxml_child(xmesh, "ain");
     if(xdata) sscanf(xdata->txt, "%lf", &mesh->ain);

     if(!setup->compress){
       xdata = ezxml_child(xmesh, "da");
       if(xdata) sscanf(xdata->txt, "%lf", &mesh->da);
     }

/*   read flow data
     -------------- */
     if(setup->compress) xreg = ezxml_child(xcomp, "flow");
     else                xreg = ezxml_child(xinco, "flow");

     xbc   = ezxml_child(xreg, "bc");
     xgas  = ezxml_child(xreg, "gaslaw");

     if(setup->compress){
       xdata = ezxml_child(xbc, "ptin");
       if(xdata) sscanf(xdata->txt, "%lf", &flow->ptin);

       xdata = ezxml_child(xbc, "min");
       if(xdata) sscanf(xdata->txt, "%lf", &flow->min);

       xdata = ezxml_child(xbc, "pex");
       if(xdata) sscanf(xdata->txt, "%lf", &flow->pex);

       xdata = ezxml_child(xgas, "gamma");
       if(xdata) sscanf(xdata->txt, "%lf", &flow->gamma);

       xdata = ezxml_child(xgas, "rgas");
       if(xdata) sscanf(xdata->txt, "%lf", &flow->rgas);

       xdata = ezxml_child(xgas, "cp");
       if(xdata) sscanf(xdata->txt, "%lf", &flow->cp);

       xdata = ezxml_child(xgas, "upwind");
       if(xdata) sscanf(xdata->txt, "%d", &flow->upwind);
     }
     else{
       xdata = ezxml_child(xbc, "ptin");
       if(xdata) sscanf(xdata->txt, "%lf", &flow->ptin);

       xdata = ezxml_child(xbc, "pex");
       if(xdata) sscanf(xdata->txt, "%lf", &flow->pex);

       xdata = ezxml_child(xreg, "rho");
       if(xdata) sscanf(xdata->txt, "%lf", &flow->rho);
     }

/*   read solver data
     ---------------- */
     if(setup->compress) xreg = ezxml_child(xcomp, "solver");
     else                xreg = ezxml_child(xinco, "solver");

     if(setup->couple){
       xmeth = ezxml_child(xreg, "couple");
       xdata = ezxml_child(xmeth, "niters");
       if(xdata) sscanf(xdata->txt, "%d", &solver->niters);

       xdata = ezxml_child(xmeth, "presim");
       if(xdata) sscanf(xdata->txt, "%d", &solver->presim);

       xdata = ezxml_child(xmeth, "relaxc");
       if(xdata) sscanf(xdata->txt, "%lf", &solver->relaxc);

       xdata = ezxml_child(xmeth, "rrampc");
       if(xdata) sscanf(xdata->txt, "%d", &solver->rrampc);

       xdata = ezxml_child(xmeth, "tolerc");
       if(xdata) sscanf(xdata->txt, "%lf", &solver->tolerc);
     }
     else{
       xmeth = ezxml_child(xreg, "simple");
       xdata = ezxml_child(xmeth, "niters");
       if(xdata) sscanf(xdata->txt, "%d", &solver->niters);

       xdata = ezxml_child(xmeth, "relaxu");
       if(xdata) sscanf(xdata->txt, "%lf", &solver->relaxu);

       xdata = ezxml_child(xmeth, "relaxp");
       if(xdata) sscanf(xdata->txt, "%lf", &solver->relaxp);

       xdata = ezxml_child(xmeth, "toleru");
       if(xdata) sscanf(xdata->txt, "%lf", &solver->toleru);

       xdata = ezxml_child(xmeth, "tolerp");
       if(xdata) sscanf(xdata->txt, "%lf", &solver->tolerp);
     }

/*   set file names
     -------------- */
     noz_files(nozzle);

/*   free local memory
     ----------------- */
     ezxml_free(xfile);

     return 0;
}

/******************************************************************************
 *  print_setup - print nozzle structure to stdout                            *
 ******************************************************************************/
/**
 * \brief   Dump nozzle structure
 * \details Prints nozzle structure based on set-up options to stdout
 *
 * \param[in]
 *          nozzle  \ref NOZZLE structure to be printed
 *
 * \return  error code
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    2nd Dec 2022
 */

void print_setup(NOZZLE *nozzle)
{
     SETUP  *setup  = nozzle->setup;
     MESH   *mesh   = nozzle->mesh;
     FLOW   *flow   = nozzle->flow;
     SOLVER *solver = nozzle->solver;
     FILES  *files  = nozzle->files;
     ARRAYS *arrays = nozzle->arrays;

/*   banner
     ------ */
     printf("\nNozzle setup and options\n");
     printf(  "------------------------\n");

/*   setup
     ----- */
     printf("setup:\n");
     printf("\tcasename    = %s\n",  setup->casename);
     printf("\tflow regime = ");
     setup->compress > 0 ? printf("compressible\n") : printf("incompressible\n");
     printf("\tsolver      = ");
     setup->couple   > 0 ? printf("coupled\n")      : printf("SIMPLE\n");
     printf("\n");

/*   mesh
     ---- */
     printf("mesh:\n");
     printf("\n");
     printf("\tnumber of stations     = % d\n", arrays->nx);
     if(setup->compress){
       printf("\tarea scaling factor    = % lf\n", mesh->ain);
     }
     else{
       printf("\tinlet area             = % lf\n", mesh->ain);
       printf("\tarea contraction ratio = % lf\n", mesh->da);
     }
     printf("\n");

/*   flow
     ---- */
     printf("flow:\n");
     if(setup->compress){
       printf("\tinlet total pressure            = %lf\n", flow->ptin);
       printf("\tinlet total temperature         = %lf\n", flow->ttin);
       printf("\tinlet Mach number               = %lf\n", flow->min);
       printf("\texit static pressure            = %lf\n", flow->pex);
       printf("\tratio of specific heats (gamma) = %lf\n", flow->gamma);
       printf("\tspecific heat capacity (cp)     = %lf\n", flow->cp);
       printf("\tidealgas constant (R)           = %lf\n", flow->rgas);
       printf("\tdensity upwinding scheme        = %d\n",  flow->upwind);
     }
     else{
       printf("\tconstant density     = %lf\n", flow->rho);
       printf("\tinlet total pressure = %lf\n", flow->ptin);
       printf("\texit static pressure = %lf\n", flow->pex);
     }
     printf("\n");

/*   solver
     ------ */
     printf("solver:\n");
     printf("\tnumber of iterations    = %d\n", solver->niters);
     if(setup->couple){
       printf("\trelaxation factor       = %lf\n", solver->relaxc);
       printf("\tpre-SIMPLE iterations   = %d\n",  solver->presim);
       printf("\trelaxation factor ramp  = %d\n",  solver->rrampc);
       printf("\tconvergence tolerence   = %le\n", solver->tolerc);
     }
     else{
       printf("\tu-relaxation factor     = %lf\n", solver->relaxu);
       printf("\tp-relaxation factor     = %lf\n", solver->relaxp);
       printf("\tu-convergence tolerence = %le\n", solver->toleru);
       printf("\tp-convergence tolerence = %le\n", solver->tolerp);
     }
     printf("\n");

/*   files
     ----- */
     printf("files:\n");
     printf("\tconvergence history file = %s\n", files->histfile);
     printf("\n");

}

/******************************************************************************
 *  print_flow - output nozzle solution for plotting in gnuplot               *
 ******************************************************************************/
/**
 * \brief   Output nozzle solution for plotting in gnuplot 
 * \details Output nozzle solution for plotting in gnuplot, set filename to
 *          NULL to print to stdout.
 *
 * \param[in]
 *          nozzle   \ref NOZZLE structure to hold nozzle parameters.
 * \param[in]
 *          filename name of file for output or NULL to print to stdout
 *
 * \return  error code: 0 if successful, 1 if failure occurred
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    7th Dec 2022
 */

int  print_flow(NOZZLE *nozzle, char *filename)
{
     char    routine[] = "printf_flow";
     SETUP   *setup  = nozzle->setup;
     ARRAYS  *arrays = nozzle->arrays;
     int32_t nx     = arrays->nx;
     FILE    *fptr;
     int32_t i;
     double  x;

/*   open file if specified
     ---------------------- */
     if(filename){
       if(!(fptr = fopen(filename, "w"))){
         error_message(routine, __FILE__, __LINE__, "error opening output file: %s",filename);
         return 1;
       }
       printf("writing output to file: %s\n",filename);
     }
     else
       fptr = stdout;

/*   print column headers
     -------------------- */
     if(setup->compress){
       fprintf(fptr, " \"index\",   \"x\",              \"area\",           \"density\",        \"velocity\",       \"pressure\",      \"temperature\",     \"Mach number\"\n");
     }
     else{
       fprintf(fptr, " \"index\",   \"x\",              \"area\",           \"density\",        \"velocity\",       \"pressure\"\n");
     }

/*   print nozzle solution
     --------------------- */
     if(setup->compress){
       for(i=0; i<nx; i++){
         x = (double) i/ (double) (nx - 1);
         fprintf(fptr, "%8d, % 16.12lf, % 16.12lf, % 16.12lf, % 16.12lf, % 16.12lf, % 16.12lf, % 16.12lf\n",
                 i+1, x, arrays->a[i], arrays->rho[i], arrays->u[i],
                         arrays->p[i], arrays->t[i],   arrays->M[i]);
       }
     }
     else{
       for(i=0; i<nx; i++){
         x = (double) i/ (double) (nx - 1);
         fprintf(fptr, "%8d, % 16.12lf, % 16.12lf, % 16.12lf, % 16.12lf, % 16.12lf\n",
                 i+1, x, arrays->a[i], arrays->rho[i], arrays->u[i], arrays->p[i]);
       }
     }

/*   close file
     ---------- */
     if(filename) fclose(fptr);

     return 0;
}

