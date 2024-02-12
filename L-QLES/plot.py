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
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

###########################################################
#   plot 1d solution vector                               #
###########################################################

def plotsol_1d(x, s, status):
    plt.plot(x, s, 'o', linestyle = 'solid',  color = 'blue')
    if not status: plt.title ("Solution status = False")
    plt.xlabel("Mesh coordinate")
    plt.ylabel("scalar value")
    plt.grid()
    plt.show()
    return

###########################################################
#   plot 2d solution field and mesh                       #
###########################################################

def plotsol_2d(x, y, s, status, splitp):

#   solution field
    if not splitp:
       fig = plt.figure(figsize=(12,6));
       ax = plt.subplot(121)
#      ax.set_anchor('N')
    else:
       fig = plt.figure(figsize=(6,6));
       ax = plt.subplot()
    ax.set_aspect('equal')
    u = s.reshape(len(y),len(x))
    plt.contourf(x, y, u)
    clb=plt.colorbar(orientation = 'horizontal')
    clb.set_label('Scalar values', fontsize=14);
    plt.title("Scalar Field\n", fontsize=14)
    plt.xlabel("x-coordinate", fontsize=14)
    plt.ylabel("y-coordinate", fontsize=14)
    if splitp: plt.show()

#   mesh
    if not splitp:
       ax = plt.subplot(122)
#      ax.set_anchor('N')
    else:
       fig = plt.figure(figsize=(6,6));
       ax = plt.subplot()
    ax.set_aspect('equal')

    X, Y = np.meshgrid(x, y)
    plt.plot(X,Y,  linestyle = 'solid',  color = 'blue');
    plt.plot(np.transpose(X), np.transpose(Y), linestyle = 'solid',  color = 'blue');
    plt.title("Mesh", fontsize=14)
    plt.xlabel("x-coordinate")
    plt.ylabel("y-coordinate")
    plt.axis('off')
    plt.tight_layout()

#   plot
    if not splitp:
       if not status: fig.suptitle ("Solution status = False")
    plt.show()

    return

###########################################################
#   plot 3d solution field and mesh                       #
###########################################################

def plotsol_3d(x, y, z, s, cut, status, splitp):

#   get slice
    nx = len(x)
    ny = len(y)
    nz = len(z)
    U = s.reshape(nz, ny, nx)
    X, Y, Z = np.meshgrid(z, y, x)

    if cut in ("x", "X"):
       pp = int((nx+1)/2)
       xp, yp = np.meshgrid(y, z)
       up = U[:, :, pp]
       xpl = "y-coordinate"
       ypl = "z-coordinate"
    elif cut in ("y", "Y"):
       pp = int((ny+1)/2)
       xp, yp = np.meshgrid(x, z)
       up = U[:, pp, :]
       xpl = "x-coordinate"
       ypl = "z-coordinate"
    elif cut in ("z", "Z"):
       pp = int((nz+1)/2)
       xp, yp = np.meshgrid(x, y)
       up = U[pp, :, :]
       xpl = "x-coordinate"
       ypl = "y-coordinate"
    else:
       print("\ninvalid cut: ", cut)
       return

#   solution field
    if not splitp:
       fig = plt.figure(figsize=(12,6))
       ax = plt.subplot(121)
       ax.set_anchor('N')
    else:
       fig = plt.figure(figsize=(6,6));
       ax = plt.subplot()
    ax.set_aspect('equal')

    plt.contourf(xp, yp, up)
    clb=plt.colorbar(orientation = 'horizontal')
    clb.set_label('Scalar values', fontsize=14);
    plt.title("Scalar Field\n", fontsize=14)
    plt.xlabel(xpl, fontsize=14)
    plt.ylabel(ypl, fontsize=14)
    if splitp: plt.show()

#   mesh
    if not splitp:
       ax = plt.subplot(122)
       ax.set_anchor('N')
    else:
       fig = plt.figure(figsize=(6,6));
       ax = plt.subplot()
    ax.set_aspect('equal')

    plt.plot(xp, yp,  linestyle = 'solid',  color = 'blue');
    plt.plot(np.transpose(xp), np.transpose(yp), linestyle = 'solid',  color = 'blue');
    plt.title("Mesh", fontsize=14)
    plt.xlabel(xpl)
    plt.ylabel(ypl)
    plt.axis('off')
#   plt.tight_layout()

#   plot
    if not splitp:
       if not status: fig.suptitle ("Solution status = False")
    plt.show()

    return

###########################################################
#   plot matrix                                           #
###########################################################
def plotmat(A, splitp):
    if not splitp:
       fig = plt.figure(figsize=(12,6))
       ax = plt.subplot(121)
    plt.imshow(A,interpolation='none');
    clb=plt.colorbar();
    clb.set_label('Matrix elements values');
    plt.title('Matrix values',fontsize=16)
    plt.tight_layout()
    plt.grid()
    if splitp: plt.show()

    if not splitp:
       ax = plt.subplot(122)
    SP = csr_matrix(A)
    plt.spy(SP, markersize=2)
    plt.title('Matrix sparsity pattern\n',fontsize=14)
    plt.tight_layout()
    plt.grid()
    plt.show()
    return
