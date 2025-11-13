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
 * \file   c2py.h
 * \brief  Include file for interfacing C code to Python
 * \author Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date   19th Dec 2022
 */

#ifndef C2PY_H
#define C2PY_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <nozzle.h>

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
 *  Prototypes for calling from Python (py_interface.c)                       *
 ******************************************************************************/
#ifdef _MSC_VER
    #define EXPORT_SYMBOL __declspec(dllexport)
#else
    #define EXPORT_SYMBOL
#endif

EXPORT_SYMBOL int  c2py_noz_create (void**);
EXPORT_SYMBOL void c2py_noz_free   (void*);
EXPORT_SYMBOL int  c2py_read_input (const char*, void*);
EXPORT_SYMBOL int  c2py_noz_init   (void*);
EXPORT_SYMBOL int  c2py_run_c      (void*);
EXPORT_SYMBOL int  c2py_get_smat   (void*, SMATRIX*);
EXPORT_SYMBOL int  c2py_put_soln   (void*, double*);
EXPORT_SYMBOL int  c2py_hist_chk   (void*, bool*);
EXPORT_SYMBOL int  c2py_save_soln  (void*);
EXPORT_SYMBOL int  c2py_print_soln (void*);
EXPORT_SYMBOL int  c2py_set_bcs    (void*);


#ifdef __cplusplus
}
#endif

#endif  /* C2PY_H */

