#!/usr/bin/python3

#################################################################################################
#                                                                                               #
# Copyright 2022 Rolls-Royce plc                                                                #
#                                                                                               #
# Redistribution and use in source and binary forms, with or without modification, are          #
# permitted provided that the following conditions are met:                                     #
#                                                                                               #
# 1. Redistributions of source code must retain the above copyright notice, this list of        #
#    conditions and the following disclaimer.                                                   #
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of     #
#    conditions and the following disclaimer in the documentation and/or other materials        #
#    provided with the distribution.                                                            #
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to  #
#    endorse or promote products derived from this software without specific prior written      #
#    permission.                                                                                #
#                                                                                               #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS   #
# OR IMPLIED  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF              #
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE    #
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,     #
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE #
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED    #
# AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING     #
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED  #
# OF THE POSSIBILITY OF SUCH DAMAGE.                                                            #
#                                                                                               #
#################################################################################################

##
# @file    run_pysolve.py
# @brief   Hybrid Python-C code using C matrix assembly and Python matrix solver
# @details Sparse matrices are passed across the interface from the C assembly
#          routines. These can be for either SIMPLE of coupled discretisations,
#          Python is unware of the matrix type and simply solves Ax=b and
#          returns x to the C library.
#          In this context, hybrid refers a classical C-Python nozzle solver.
#          This provides the template for replacing the Python matrix solver
#          with an emulated quantum linear equation solver.
# @author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
# @date    24th Dec 2022
##

import sys, getopt
from   ctypes import *
import sys
import pathlib
import numpy as np
from   scipy import sparse
from   scipy.sparse.linalg import spsolve

#     Python ctypes sparse matrix structure - this must match the
#     C \ref SMATRIX structure
#     -----------------------------------------------------------
class SMATRIX(Structure):
      _fields_ = [("nr",  c_int32),
                  ("nc",  c_int32),
                  ("nnz", c_int32),
                  ("col", POINTER(c_int32)),
                  ("rst", POINTER(c_int32)),
                  ("val", POINTER(c_double)),
                  ("b",   POINTER(c_double))
                 ]
      _defaults_ = {"nr" : 0,  
                    "nc" : 0,  
                    "nnz": 0,
                    "col": None, 
                    "rst": None, 
                    "val": None, 
                    "b"  : None}

###########################################################
#   line number macro                                     #
###########################################################
def LINE():
    return sys._getframe(1).f_lineno

###########################################################
#   exit on error                                         #
###########################################################
def error_exit(error, line):
    message = error + ' in ' + str(sys.argv[0]) + ':' + str(line)
    sys.exit(str(message))

###########################################################
#   get command line arguments                            #
###########################################################
def read_args():
    inputfile = ''
    argv = sys.argv[1:]

    try:
       opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
       print('Input error:\n\t%s -h -i <inputfile>' % str(sys.argv[0]))
       sys.exit(2)
    for opt, arg in opts:
       if opt == '-h':
          print('Usage:\n\t%s -h -i <inputfile>' % str(sys.argv[0]))
          sys.exit()
       elif opt in ("-i", "--ifile"):
          inputfile = arg

    return inputfile

###########################################################
#   main routine                                          #
###########################################################
if __name__ == "__main__":
    libname = pathlib.Path().absolute()
    print("libname: ", libname)

#   load shared object library using ctypes
#   ---------------------------------------
    if sys.platform.startswith("win"):
       c_lib = CDLL(libname / "nozzle.dll")
    else:
       c_lib = CDLL(libname / "../lib/libnozzle.so")

#   get name of input file
#   ----------------------
    inputfile = read_args()
    if not inputfile:
       error_exit("No input file specified", LINE())

#   create empty nozzle structure
#   -----------------------------
    noz_cptr = c_void_p()
    rc = c_lib.c2py_noz_create(byref(noz_cptr))
    if rc:
       error_exit("Traceback error", LINE())

#   read input file
#   ---------------
    i_cptr = inputfile.encode('utf-8')
    c_lib.c2py_read_input.argtypes = [c_char_p, c_void_p]
    rc = c_lib.c2py_read_input(i_cptr, noz_cptr)
    if rc:
       error_exit("Traceback error", LINE())

#   initialise nozzle
#   -----------------
    rc = c_lib.c2py_noz_init(noz_cptr)
    if rc:
       error_exit("Traceback error", LINE())

#   run hybrid C-Python nozzle solver
#   ---------------------------------
    S = SMATRIX()
    while(True):
       rc = c_lib.c2py_get_smat(noz_cptr, byref(S))
       if rc:
          error_exit("Traceback error", LINE())

       col = np.ctypeslib.as_array(S.col, shape=(S.nnz,))
       rst = np.ctypeslib.as_array(S.rst, shape=(S.nr+1,))
       val = np.ctypeslib.as_array(S.val, shape=(S.nnz,))
       b   = np.ctypeslib.as_array(S.b,   shape=(S.nr,))

       SP = sparse.csr_matrix((val, col, rst),shape=(S.nr, S.nc), dtype=np.double)
       x  = spsolve(SP, b)

       xc = x.ctypes.data_as(POINTER(c_double))
       rc = c_lib.c2py_put_soln(noz_cptr, xc)
       if rc:
          error_exit("Traceback error", LINE())

       done = c_bool()
       c_lib.c2py_hist_chk(noz_cptr, byref(done))       # checks both convergence and max iters
       if done:
          break

#   print and save solution
#   ----------------------- */
    c_lib.c2py_print_soln(noz_cptr)
    c_lib.c2py_save_soln(noz_cptr)

#   free C memory associated with nozzle
#   ------------------------------------
    c_lib.c2py_noz_free(noz_cptr)

#   exit message
#   ------------
    print("\nSuccessful completion\n")
