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
 * \file    main.c
 * \brief   Main routine for 1-dimensional nozzle
 * \details This is used to create an executable for the classical solvers.
 *          Hybrid quantum-classical solvers must provide their own C or Python
 *          main routine and load the shared object library for the nozzle solver.
 *          This routine requires the name of the input file to be provided on the
 *          command line. Hybrid solvers are responsible for setting the name of
 *          the input file.
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

#include <proto.h> 

/******************************************************************************
 *  Main routine - \cond and \endcond excludes main from Doxgen processing    *
 ******************************************************************************/
/// \cond
int  main(int argc, char **argv)
{
     int    arg, rc;
     bool   verbose = false;
     char   *inputfile = NULL;
     NOZZLE nozzle;
     SETUP  *setup;

/*   banner
     ------ */
     printf("\n***********************************************************\n");
     printf(  "*                                                         *\n");
     printf(  "*    1D Nozzle - Incompressible and Compressible Flow     *\n");
     printf(  "*                                                         *\n");
     printf(  "***********************************************************\n");

/*   process command line arguments
     ------------------------------ */
     if(argc > 0)
     {
       for(arg=1; arg<argc; arg++)
       {

/*       help
         ---- */
         if(!strncmp(argv[arg], "-h", 2)){
           printf("\nSyntax:\n\n");
           printf("%s [-h] [-i input file] [-v level]\n\n",argv[0]);
           printf("    -h  this help menu\n");
           printf("    -i  name of input file (default = input.xml)\n");
           printf("    -v  verbosity level (default = 0)\n");
           printf("\n");
           exit(0);
         }

/*       process flags
         ------------- */
         if(!strncmp(argv[arg], "-i", 2)) inputfile = strdup(argv[arg+1]);
         if(!strncmp(argv[arg], "-v", 2)) verbose   = true;
       }
     }

/*   check input file
     ---------------- */
     if(!inputfile){
       error_message("main", __FILE__, __LINE__, "%s", "no input file specified");
       exit(1);
     }

/*   create empty nozzle structure
     ----------------------------- */
     rc = noz_create(&nozzle);

/*   read and echo input file using XML
     ---------------------------------- */
     rc = read_input_xml(inputfile, &nozzle);
     free(inputfile);
     if(rc){
       error_traceback("main", __FILE__, __LINE__);
       exit(1);
     }
     print_setup(&nozzle);

/*   initialise nozzle mesh and flow field
     ------------------------------------- */
     noz_initialise(&nozzle);
     print_flow(&nozzle, NULL);

/*   launch solver
     ------------- */
     setup = nozzle.setup;
     if(setup->couple) rc = couple_run(&nozzle);
     else              rc = simple_run(&nozzle);

     if(rc){
       error_traceback("main", __FILE__, __LINE__);
       exit(1);
     }

     print_flow(&nozzle, NULL);

/*   tidy-up
     ------- */
     noz_free(&nozzle);
}
/// \endcond
