#!/usr/bin/python3

# Copyright 2022 Rolls-Royce plc

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import sys, getopt
import csv
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

###########################################################
#   get command line arguments                            #
###########################################################
def read_args(argv):
    mfile = ''
    bfile = ''
    xfile = ''
    try:
       opts, args = getopt.getopt(argv,"hm:b:x:",["mfile=","bfile=","xfile="])
    except getopt.GetoptError:
       print ('plot-mat.py -m <matfile> -b <rhsfile> -x <solfile>')
       sys.exit(2)
    for opt, arg in opts:
       if opt == '-h':
          print ('plot-mat.py -m <matfile> -b <rhsfile> -x <solfile>')
          sys.exit()
       elif opt in ("-m", "--mfile"):
          mfile = arg
       elif opt in ("-b", "--bfile"):
          bfile = arg
       elif opt in ("-x", "--xfile"):
          xfile = arg

    return mfile, bfile, xfile

###########################################################
#   read vector file                                      #
###########################################################
def read_vec(filename):

#   read vector data

    file = open(filename, "rb")

    nv = np.fromfile(file, dtype=np.long, count=1)
    v  = np.fromfile(file, dtype=np.double)
    file.close()

    print("vector:");
    print(v)

    return nv, v

###########################################################
#   read matrix file                                      #
###########################################################
def read_mat(filename):

#   read matrix data

    file = open(filename, "rb")

    real = np.fromfile(file, dtype=np.bool, count=1)
    nrow = np.fromfile(file, dtype=np.long, count=1)
    ncol = np.fromfile(file, dtype=np.long, count=1)
    nonz = np.fromfile(file, dtype=np.long, count=1)

    nr  = nrow[0]
    nc  = ncol[0]
    nnz = nonz[0]
    print("matrix:")
    print(real, nrow, ncol, nnz)

    rval = np.fromfile(file, dtype=np.double, count=nnz)
    print(rval)

    col  = np.fromfile(file, dtype=np.long, count=nnz)
    print(col)

    rowstt = np.fromfile(file, dtype=np.long, count=nr+1)
    print(rowstt)

    file.close()

#   create sparse matrix
    S = sparse.csr_matrix((rval, col, rowstt),shape=(nr, nc), dtype=np.double)

    return S

###########################################################
#   main routine                                          #
###########################################################
def main(argv):

    nb = 0
    nx = 0

#   get filenames
    mfile, bfile, xfile = read_args(argv)

#   read vectors
    if(bfile): nb, b = read_vec(bfile)
    if(xfile): nx, x = read_vec(xfile)

#   read matrix
    if(mfile): S = read_mat(mfile)

#   plot sparsity pattern of matrix
    if(mfile):
       plt.spy(S, markersize=2)
       plt.show()

#   plot vectors
    if(bfile and xfile):
       fig, ax = plt.subplots()
       plt.title("Right hand side and solution vectors")
       plt.xlabel("Index")
       plt.ylabel("Normalised pressure")
       ax.plot(b, label = 'RHS',  linestyle = 'solid',  color = 'black')
       ax.plot(x, label = 'Soln', linestyle = 'dashed', color = 'black')
       ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
       ax.legend()
#      ax.yaxis.set_label_coords(.035, .5)
       plt.grid()
       plt.show()
       return

    if(bfile):
       plt.title("Right hand side vector")
       plt.xlabel("index")
       plt.ylabel("normalised pressure")
       plt.plot(b)
       plt.grid()
       plt.show()

    if(xfile):
       plt.title("Solution vector")
       plt.xlabel("index")
       plt.ylabel("normalised pressure")
       plt.plot(x)
       plt.grid()
       plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
