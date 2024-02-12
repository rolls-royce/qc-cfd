#!/usr/bin/python3

#################################################################################################
#                                                                                               #
# Code for generating 1D 2d and 3D Laplacian operators with representative boundary conditions  #
# for testing Quantum Linear Equation Solvers                                                   #
#                                                                                               #
# Copyright 2024 Rolls-Royce plc                                                                #
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

import os.path
import sys, getopt
import numpy as np
from numpy import linalg as lin
from scipy import sparse
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

from mesh    import parse_meshfile, generate_mesh
from matvec  import matvec_1d,  matvec_2d,  matvec_3d
from plot    import plotsol_1d, plotsol_2d, plotsol_3d, plotmat
from reorder import reorder
from save    import case_save_npz, case_save_bin

###########################################################
#   read command line arguments                           #
###########################################################
def read_args(argv):
    inputfile=''
    mplot = False
    splot = False
    cut3d = "x"
    degen = False
    eigen = False
    order = False
    psplt = False
    
    try:
       opts, args = getopt.getopt(argv,"hmsdejri:c:",["i=","c="])
    except getopt.GetoptError:
       print ('\nusage error, use:')
       print ('\tlaplace.py -h for help\n')
       sys.exit(2)

    for opt, arg in opts:
       if opt == '-h':
          print ('\nusage:')
          print ('\tl-qles.py -i <input file> {-c <x,y,z>} {-d} {-e} {-h} {-j} {-m} {-r} {-s}\n')
          print ('\t\t-i {name of input file}')
          print ('\t\t-c {x,y,z} cut slice of 3D solution to be plotted, default = x')
          print ('\t\t-d allow degnerate matrices, default = False')
          print ('\t\t-e calculate eigenvalues and condition number, default = False')
          print ('\t\t-h help menu')
          print ('\t\t-j split plots into separate windows for saving, default is single window')
          print ('\t\t-m plot matrix, default = False')
          print ('\t\t-r reorder matrix and RHS to use shell ordering of mesh, default = False')
          print ('\t\t-s plot solutons and mesh, default = False')
          sys.exit()
       elif opt in ("-c", "--c"):
          cut3d = arg
       elif opt in ("-d", "--d"):
          degen = True
       elif opt in ("-e", "--e"):
          eigen = True
       elif opt in ("-i", "--i"):
          inputfile = arg
       elif opt in ("-m", "--m"):
          mplot = True
       elif opt in ("-j", "--j"):
          psplt = True
       elif opt in ("-r", "--r"):
          order = True
       elif opt in ("-s", "--s"):
          splot = True

    if not os.path.isfile(inputfile):
       print('\nfile', inputfile, 'does not exist\n') 
       sys.exit(3)

    return inputfile, mplot, splot, cut3d, degen, eigen, order, psplt


###########################################################
#   main routine                                          #
###########################################################
def laplace(argv):

#   read input file 
    inputfile, mplot, splot, cut, degen, eigen, order, psplt = read_args(argv)
    casename, ndims, rdict = parse_meshfile(inputfile)
#   print("rdict:\n", rdict)

#   generate mesh coordinates
    if ndims > 0: x = generate_mesh(rdict['x'])
    if ndims > 1: y = generate_mesh(rdict['y'])
    if ndims > 2: z = generate_mesh(rdict['z'])

#   generate matrix and rhs
    if ndims == 1:
       a, b = matvec_1d(x, rdict['x'], rdict['force'], degen)
    elif ndims == 2:
       a, b = matvec_2d(x, y, rdict['x'], rdict['y'], rdict['force'], degen)
    elif ndims == 3:
       a, b = matvec_3d(x, y, z, rdict['x'], rdict['y'], rdict['z'], rdict['force'], degen)

#   reorder: solve PAP^{-1} Px = Pb where P is a permutation matrix, need to permute solution later
    if order:
       if ndims == 1:
          q, a, b  = reorder(a, b, len(x), 1, 1)
       elif ndims == 2:
          q, a, b = reorder(a, b, len(x), len(y), 1)
       elif ndims == 3:
          q, a, b = reorder(a, b, len(x), len(y), len(z))
    else:
       q = np.identity(np.shape(a)[0])

#   solve (scipy sparse linalg solver is more reliable than numpy linalg lin.solve(a,b))
    SP = csr_matrix(a)
    s  = spsolve(SP, b)
    status = np.allclose(np.dot(a, s), b)
    print("solution status = ", status)

#   permute solution
    if order:
       s = np.matmul(q, s)

#   save python npz and C binary files
    case_save_npz(a, b, s, q, status, degen, order, casename)
    case_save_bin(a, b, s, q, status, degen, order, casename)

#   eigen analysis - use symmetrised Hernmitian matrix
    if eigen:
       print("\ncalculating eigenvalues:")
       nt    = a.shape[0]
       asym  = np.block([[np.zeros([nt,nt]),a],[np.transpose(a),np.zeros([nt,nt])]])
#      kappa = np.linalg.cond(asym)
       evals = np.linalg.eigvals(asym)
       emin  = min(np.abs(evals))
       emax  = max(np.abs(evals))
       kappa = emax/emin
       print("\tcondition number:  ",kappa)
       print("\tmin abs eigenvalue:", min(abs(evals)))
       print("\tmax abs eigenvalue:", max(abs(evals)))

#   plot
    if splot:
       if ndims == 1:
          plotsol_1d(x, s, status)
       if ndims == 2:
          plotsol_2d(x, y, s, status, psplt)
       if ndims == 3:
          plotsol_3d(x, y, z, s, cut, status, psplt)

    if mplot: plotmat(a, psplt)

###########################################################
#   call main                                             #
###########################################################
if __name__ == "__main__":
    laplace(sys.argv[1:])

