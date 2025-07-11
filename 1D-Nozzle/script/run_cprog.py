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
# @file    run_cprog.py
# @brief   Simple python interface to get input file and run C code
# @details Tests use of Python ctype interface for receiving and passing
#          pointers to a C \ref NOZZLE structure and processing return
#          codes. Use of ctypes functionality uses ctypes prefix to 
#          aid understanding.
# @author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
# @date    19th Dec 2022
##

import sys, getopt
import ctypes
import sys
import pathlib
import time

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
       c_lib = ctypes.CDLL(libname / "nozzle.dll")
    else:
       c_lib = ctypes.CDLL(libname / "../lib/libnozzle.so")

#   get name of input file
#   ----------------------
    inputfile = read_args()
    if not inputfile:
       error_exit("No input file specified", LINE())

#   create empty nozzle structure
#   -----------------------------
    noz_cptr = ctypes.c_void_p()
    rc = c_lib.c2py_noz_create(ctypes.byref(noz_cptr))
    if rc:
       error_exit("Traceback error", LINE())

#   read input file
#   ---------------
    i_cptr = inputfile.encode('utf-8')
    c_lib.c2py_read_input.argtypes = [ctypes.c_char_p, ctypes.c_void_p]
    rc = c_lib.c2py_read_input(i_cptr, noz_cptr)
    if rc:
       error_exit("Traceback error", LINE())

#   initialise nozzle
#   -----------------
    rc = c_lib.c2py_noz_init(noz_cptr)
    if rc:
       error_exit("Traceback error", LINE())

#   run entirely using C library
#   ----------------------------
    t1 = time.perf_counter()
    rc = c_lib.c2py_run_c(noz_cptr)
    t2 = time.perf_counter()
    print(f"\nC solver took {t2-t1:0.4f} seconds")
    if rc:
       error_exit("Traceback error", LINE())

#   free C memory associated with nozzle
#   ------------------------------------
    c_lib.c2py_noz_free(noz_cptr)

    print("\nSuccessful completion\n")
