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

import numpy as np

m = lambda i, j, k, ni, nj, nk: min(i,ni-1) + ni*min(j,nj-1) + ni*nj*min(k,nk-1)

###########################################################
#   closest neighbour reordering                          #
###########################################################
def reorder(a, b, ni, nj, nk):

    na = ni*nj*nk                    # note nj and/or nk = 1 for 1D and 2D meshes
    r = np.zeros(shape=(na), dtype=np.int)
    f = np.zeros(shape=(na), dtype=np.int)
    r[0] = 0 
    f[0] = 1 

#   not the most elegant or efficient method but good enough for demos
    n = 0
    while n<na:
       f = -f
       for k in range(0, nk):
          for j in range(0, nj):
             for i in range(0, ni):
                mp = m(i, j, k, ni, nj, nk)
                if f[mp] < 0:
                   for ka in range(0,min(2,nk)):
                      for ja in range(0,min(2,nj)):
                         for ia in range(0,min(2,ni)):
                            mn =  m(i+ia, j+ja, k+ka, ni, nj, nk)
                            if r[mn] == 0:
                               r[mn] = n
                               f[mn] = 1
                               n = n+1
                   f[mp] = 0
#   print(r)

#   invert mapping
    ma = np.zeros(shape=(na), dtype=np.int)
    n = 0
    for i in range(0, na):
       ma[r[i]] = n
       n = n+1
#   print("mapping:\n", ma)

#   hand coded post permutation to a: i.e. AQ.Qx = b superseded below
#   pa = np.zeros(shape=(na, na))
#   for i in range(0, na):
#      for j in range(0, na):
#         pa[i][j] = a[i][ma[j]]

#   permutation operators - not p and q are 1-sparse and hence their own inverse
    p = f = np.zeros(shape=(na, na))
    q = f = np.zeros(shape=(na, na))
    for i in range(0, na):
       p[i][ma[i]] = 1
       q[ma[i]][i] = 1

#   apply reordering as preconditioner: i.e. PAQ.Qx = Pb   Not the usual form of preconditioning
    aq  = np.matmul(a, q)
    paq = np.matmul(p, aq)
    pb  = np.matmul(p, b)
         
    return q, paq, pb

