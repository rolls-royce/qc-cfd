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
 * \file   proto.h
 * \brief  C prototypes for 1D nozzle solver
 * \author Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date   1st Dec 2022
 */

#ifndef PROTO_H
#define PROTO_H

#include <stdio.h>
#include <stdbool.h>

#include <nozzle.h>

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************
 *     Prototypes for couple.c                                                 *
 *******************************************************************************/
extern int  couple_update             (NOZZLE*, double*);
extern int  couple_mat                (NOZZLE*, SMATRIX*);
extern int  couple_mat_c              (NOZZLE*, SMATRIX*);
extern int  couple_mat_i              (NOZZLE*, SMATRIX*);
extern int  couple_run                (NOZZLE*);

/*******************************************************************************
 *     Prototypes for error.c                                                  *
 *******************************************************************************/
extern int  nan_check                 (const char*, int, double*);
extern void error_message             (const char*, const char*, const int, char*, ...);
extern void error_traceback           (const char*, const char*, const int);

/*******************************************************************************
 *     Prototypes for gsl_interface.c                                          *
 *******************************************************************************/
extern int  matsol_LU                 (SMATRIX*, double**);

/*******************************************************************************
 *     Prototypes for history.c                                                *
 *******************************************************************************/
extern int  hist_print                (NOZZLE*);
extern int  hist_output               (NOZZLE*, char*, bool);
extern bool hist_check_conv           (NOZZLE*);

/*******************************************************************************
 *     Prototypes for io.c                                                     *
 *******************************************************************************/
extern int  read_input_xml            (const char*, NOZZLE*);
extern void print_setup               (NOZZLE*);
extern int  print_flow                (NOZZLE*, char*);

/*******************************************************************************
 *     Prototypes for nozzle.c                                                 *
 *******************************************************************************/
extern int  noz_create                (NOZZLE*);
extern void noz_free                  (NOZZLE*);
extern void noz_files                 (NOZZLE*);
extern int  noz_initialise            (NOZZLE*);
extern int  noz_continuity            (NOZZLE*);
extern int  noz_set_bcs               (NOZZLE*);
extern int  noz_temp_mach             (NOZZLE*);
extern int  noz_peff                  (NOZZLE*);
extern int  noz_density               (NOZZLE*);

/*******************************************************************************
 *     Prototypes for simple.c                                                 *
 *******************************************************************************/
extern int  simple_msolve             (NOZZLE*);
extern int  simple_pcoeff             (NOZZLE*, double*, double*, double*, double*);
extern int  simple_psolve             (NOZZLE*);
extern int  simple_run                (NOZZLE*);
extern int  simple_pcmat              (NOZZLE*, SMATRIX*);
extern int  simple_pcupd              (NOZZLE*, double*);

/*******************************************************************************
 *     Prototypes for sparse.c                                                 *
 *******************************************************************************/
extern void sparse_free               (SMATRIX*);
extern void sparse_print              (char*, SMATRIX);
extern void sparse_matvec             (SMATRIX*, double*, double*);


#ifdef __cplusplus
}
#endif

#endif  /* PROTO_H */

