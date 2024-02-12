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
from   scipy.sparse import csr_matrix, save_npz

###########################################################
#   save npz and npy files x=solution, not coordinates    #
###########################################################
def case_save_npz(a, b, x, q, status, degen, order, casename):

#   append _r to casename if the matrix has been reordered
    cname = casename
    if degen: cname += '_d'
    if order: cname += '_r'

#   save matrix
    filename = cname + '_mat.npz'
    print('\nsaving sparse matrix to npz file:  ', filename)
    s = csr_matrix(a)
    save_npz(filename, s)

#   save RHS
    filename = cname + '_rhs.npy'
    print('saving RHS vector to npy file:     ', filename)
    np.save(filename, b)

#   save solution
    if status:
       filename = cname + '_sol.npy'
       print('saving solution vector to npy file:', filename)
       np.save(filename, x)

#   for reordering, QLES needs Q to get x = Qx from linear solution
    if order:
       filename = cname + '_ord.npz'
       print('saving reorder matrix to npy file :', filename)
       s = csr_matrix(q)
       save_npz(filename, s)


###########################################################
#   save binary files x=solution, not coordinates         #
###########################################################
def case_save_bin(a, b, x, q, status, degen, order, casename):

#   append _r to casename if the matrix has been reordered
    cname = casename
    if degen: cname += '_d'
    if order: cname += '_r'

#   save matrix
    filename = cname + '_mat.bin'
    print('\nsaving sparse matrix to binary file:  ', filename)

    s = csr_matrix(a)
    rank = s.shape
    nr   = np.long(rank[0])
    nc   = np.long(rank[1])
    nnz  = np.long(s.nnz)

    real = np.array([True], dtype=np.bool)
    dims = np.array([nr,nc,nnz], dtype=np.long)
    rval = np.array([s.data],    dtype=np.double)
    rstt = np.array([s.indptr],  dtype=np.long) 
    col  = np.array([s.indices], dtype=np.long)

    with open(filename, "wb") as fp:
       real.tofile(fp)
       dims.tofile(fp)
       rval.tofile(fp)
       col.tofile(fp)
       rstt.tofile(fp)

#   save RHS
    filename = cname + '_rhs.bin'
    print('saving RHS vector to binary file:     ', filename)

    nb = np.array([len(b)], dtype=np.long)
    vb = np.array([b],      dtype=np.double)

    with open(filename, "wb") as fp:
       nb.tofile(fp)
       vb.tofile(fp)

#   save solution if found
    if status:
       filename = cname + '_sol.bin'
       print('saving solution vector to binary file:', filename)
 
       nx = np.array([len(x)], dtype=np.long)
       vx = np.array([x],      dtype=np.double)
    
       with open(filename, "wb") as fp:
          nx.tofile(fp)
          vx.tofile(fp)

#   for reordering, QLES needs Q to get x = Qx from linear solution
    if order:
       filename = cname + '_ord.bin'
       print('saving reorder matrix to binary file :', filename)

       s = csr_matrix(q)
       rank = s.shape
       nr   = np.long(rank[0])
       nc   = np.long(rank[1])
       nnz  = np.long(s.nnz)

       real = np.array([True], dtype=np.bool)
       dims = np.array([nr,nc,nnz], dtype=np.long)
       rval = np.array([s.data],    dtype=np.double)
       rstt = np.array([s.indptr],  dtype=np.long)
       col  = np.array([s.indices], dtype=np.long)

       with open(filename, "wb") as fp:
          real.tofile(fp)
          dims.tofile(fp)
          rval.tofile(fp)
          col.tofile(fp)
          rstt.tofile(fp)

