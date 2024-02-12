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

mij  = lambda i, j, ni:        i + ni*j
mijk = lambda i, j, k, ni, nj: i + ni*j + ni*nj*k

###########################################################
#   check consistency of boundary conditions              #
###########################################################
def check_bcs(bcs, axis):
    valid = ['D', 'N', 'R', 'S']

    if bcs[0] not in valid:
       print('\nBoundary type', bcs[0], 'in', axis, 'direction is not valid\n')
       exit(1);
    if bcs[1] not in valid:
       print('\nBoundary type', bcs[1], 'in', axis, 'direction is not valid\n')
       exit(1);

    if (bcs[0] == 'R' and bcs[1] != 'R') or (bcs[0] != 'R' and bcs[1] == 'R'):
       print('\nMismatched repeating boundary conditions:', bcs[0], bcs[1],
             'in', axis, 'direction' '\n')
       exit(1);

###########################################################
#   generate 1D matrix and rhs, solution vectors          #
###########################################################
def matvec_1d(x, xdict, f, degen):

#   check consistency of bcs
    bcs = xdict["btype"].replace(" ", "").split(',')
    check_bcs(bcs, 'x')

    repeat = False
    if bcs[0] == 'R' and bcs[1] == 'R': repeat = True

#   initialise
    nx = len(x)
    a  = np.zeros(shape=(nx, nx))
    b  = np.zeros(shape=(nx))
    dx = np.zeros(shape=(nx+1))

    for i in range(1, nx): dx[i] = x[i] - x[i-1]
    dx[0]  = dx[1]
    dx[nx] = dx[nx-1]

#   set interior matrix entries
    for i in range(1, nx-1):
       a[i][i-1] = -1.0/dx[i]
       a[i][i+1] = -1.0/dx[i+1]
       a[i][i]   = -a[i][i-1] - a[i][i+1]

#   set boundary entries
    if repeat:
       a[0][1]    = -1.0/dx[0]
       a[0][nx-1] = -1.0/dx[nx]
       a[0][0]    = -a[0][1] - a[0][nx-1]

       a[nx-1][nx-1] = a[0][0]
       a[nx-1][nx-2] = a[0][nx-1]
       a[nx-1][0]    = a[0][1]
    else:
       if bcs[0] == 'D':
         a[0][0] = a[1][1]
       elif bcs[0] == 'N' or bcs[0] == 'S':
         a[0][0] =  a[1][1]
         a[0][1] = -a[0][0]

       if bcs[1] == 'D' or bcs[1] == 'S':
         a[nx-1][nx-1] = a[nx-2][nx-2]
       elif bcs[1] == 'N':
         a[nx-1][nx-1] =  a[nx-2][nx-2]
         a[nx-1][nx-2] = -a[nx-1][nx-1]

#   set RHS state - Symmetry is special case of Neumann with zero gradient
    bval = xdict["bvalue"].replace(" ", "").split(',')
    if bcs[0] != 'S': b[0]    = float(bval[0])*a[0][0]
    if bcs[1] != 'S': b[nx-1] = float(bval[1])*a[nx-1][nx-1]
    for i in range(1, nx-1): b[i] = f

#   if degen is off, fix row i if matrix is degenerate
    if not degen:
       if bcs[0] != 'D' and bcs[1] != 'D':
          i = int(xdict["degfix"])
          a[i][i-1] = 0.0
          a[i][i+1] = 0.0
          b[i]      = b[i]*a[i][i]

#   scale to give ||a|| = 1.0 in max norm
    amax = np.amax(a)
    a = a/amax
    b = b/amax

#   debug print
#   print(np.get_printoptions())
    with np.printoptions(precision=2, suppress=True, linewidth=100):
       print(a)
       print(b)

    return a, b

###########################################################
#   generate 2D matrix and rhs, solution vectors          #
###########################################################
def matvec_2d(x, y, xdict, ydict, f, degen):

#   check consistency of bcs
    bcx = xdict["btype"].replace(" ", "").split(',')
    bcy = ydict["btype"].replace(" ", "").split(',')
    check_bcs(bcx, 'x')
    check_bcs(bcy, 'y')

    repeatx = False
    repeaty = False
    if bcx[0] == 'R' and bcx[1] == 'R': repeatx = True
    if bcy[0] == 'R' and bcy[1] == 'R': repeaty = True

#   initialise
    nx = len(x)
    ny = len(y)
    n2 = nx*ny
    a  = np.zeros(shape=(n2, n2))
    b  = np.zeros(shape=(n2))
    dx = np.zeros(shape=(nx+1))
    dy = np.zeros(shape=(ny+1))

    for i in range(1, nx): dx[i] = x[i] - x[i-1]
    dx[0]  = dx[1]
    dx[nx] = dx[nx-1]

    for j in range(1, ny): dy[j] = y[j] - y[j-1]
    dy[0]  = dy[1]
    dy[ny] = dy[ny-1]

    bvx = xdict["bvalue"].replace(" ", "").split(',')
    bvy = ydict["bvalue"].replace(" ", "").split(',')

#   set Dirichet bcs as these take precedence
    for i in range(0, nx, nx-1):
       ib = int(i/(nx-1))
       if bcx[ib] == 'D':
          ax = dx[i]
          for j in range(0, ny):
             m = mij(i, j, nx)
             if a[m][m] == 0:
                ay = 0.5*(dy[j] + dy[j+1])
                a[m][m] = 2*ay/ax + ax/dy[j] + ax/dy[j+1]
                b[m]    = float(bvx[ib])*a[m][m] 

    for j in range(0, ny, ny-1):
       jb = int(j/(ny-1))
       if bcy[jb] == 'D':
          ay = dy[j]
          for i in range(0, nx):
             m = mij(i, j, nx)
             if a[m][m] == 0:
                ax = 0.5*(dx[i] + dx[i+1])
                a[m][m] = 2*ax/ay + ay/dx[i] + ay/dx[i+1]
                b[m]    = float(bvy[jb])*a[m][m]

#   set Neumann/Symmetry bcs next as these take precedence over repeating bcs
    for i in range(0, nx, nx-1):
       ib = int(i/(nx-1))
       ia = 1-2*ib
       if bcx[ib] == 'N' or bcx[ib] =='S':
          ax = dx[i]
          for j in range(0, ny):
             m = mij(i, j, nx)
             if a[m][m] == 0:
                ay = 0.5*(dy[j] + dy[j+1])
                a[m][m]    = 2*ay/ax + ax/dy[j] + ax/dy[j+1]
                a[m][m+ia] = -a[m][m]
                if bcx[ib] == 'N': b[m] = float(bvx[ib])*a[m][m]

    for j in range(0, ny, ny-1):
       jb = int(j/(ny-1))
       ja = 1-2*jb
       if bcy[jb] == 'N' or bcy[jb] =='S':
          ay = dy[j]
          for i in range(0, nx):
             m = mij(i, j, nx)
             if a[m][m] == 0:
                ax = 0.5*(dx[i] + dx[i+1])
                a[m][m]       = 2*ax/ay + ay/dx[i] + ay/dx[i+1]
                a[m][m+ja*nx] = -a[m][m]
                if bcy[jb] == 'N': b[m] = float(bvy[jb])*a[m][m]

#   set dx and dy for repeating bcs
    if repeatx:
       dx[0]  = dx[nx-1]
       dx[nx] = dx[1]
    if repeaty:
       dy[0]  = dy[ny-1]
       dy[ny] = dy[1]

#   set interior matrix entries
    for i in range(0, nx):
       ax = 0.5*(dx[i] + dx[i+1])
       for j in range(0, ny):
          m  = mij(i, j, nx)
          if a[m][m] == 0:
             ay = 0.5*(dy[j] + dy[j+1])
             me = mij(i+1, j,   nx)
             mw = mij(i-1, j,   nx)
             mn = mij(i,   j+1, nx)
             ms = mij(i,   j-1, nx)

             if repeatx:
                if i == 0:     mw = mij(nx-1, j, nx)
                if i == nx-1:  me = mij(0,    j, nx)
             else:
                if i == 0:     mw = m                        # should be redundant trap as only repeats unset
                if i == nx-1:  me = m
             if repeaty:
                if j == 0:     ms = mij(i, ny-1, nx)
                if j == ny-1:  mn = mij(i, 0,    nx)
             else:
                if j == 0:     ms = m                        # should be redundant trap as only repeats unset
                if j == ny-1:  mn = m

             a[m][mw] = -ay/dx[i]
             a[m][me] = -ay/dx[i+1]
             a[m][ms] = -ax/dy[j]
             a[m][mn] = -ax/dy[j+1]
             a[m][m]  = -a[m][me] - a[m][mw] - a[m][ms] - a[m][mn]
             b[m]     = f*ax*ay

#   if degen is off, fix row mij(i,j,nx)  if matrix is degenerate
    if not degen:
       if bcx[0] != 'D' and bcx[1] != 'D':
          i = int(xdict["degfix"])
          if bcy[0] != 'D' and bcy[1] != 'D':
             j = int(ydict["degfix"])
             m  = mij(i, j, nx)
             me = mij(i+1, j,   nx)
             mw = mij(i-1, j,   nx)
             mn = mij(i,   j+1, nx)
             ms = mij(i,   j-1, nx)

             a[m][mw] = 0.0
             a[m][me] = 0.0
             a[m][ms] = 0.0
             a[m][mn] = 0.0
             b[m]     = b[m]*a[m][m]


#   scale to give ||a|| = 1.0 in max norm
    amax = np.amax(a)
    a = a/amax
    b = b/amax

#   debug print
#   print(np.get_printoptions())
    with np.printoptions(precision=2, suppress=True, linewidth=100):
       print(a)

    return a, b

###########################################################
#   generate 3D matrix and rhs, solution vectors          #
###########################################################
def matvec_3d(x, y, z, xdict, ydict, zdict, f, degen):

#   check consistency of bcs
    bcx = xdict["btype"].replace(" ", "").split(',')
    bcy = ydict["btype"].replace(" ", "").split(',')
    bcz = zdict["btype"].replace(" ", "").split(',')
    check_bcs(bcx, 'x')
    check_bcs(bcy, 'y')
    check_bcs(bcz, 'z')

    repeatx = False
    repeaty = False
    repeatz = False
    if bcx[0] == 'R' and bcx[1] == 'R': repeatx = True
    if bcy[0] == 'R' and bcy[1] == 'R': repeaty = True
    if bcz[0] == 'R' and bcz[1] == 'R': repeatz = True

#   initialise
    nx = len(x)
    ny = len(y)
    nz = len(z)
    n3 = nx*ny*nz
    a  = np.zeros(shape=(n3, n3))
    b  = np.zeros(shape=(n3))
    dx = np.zeros(shape=(nx+1))
    dy = np.zeros(shape=(ny+1))
    dz = np.zeros(shape=(nz+1))

    for i in range(1, nx): dx[i] = x[i] - x[i-1]
    dx[0]  = dx[1]
    dx[nx] = dx[nx-1]

    for j in range(1, ny): dy[j] = y[j] - y[j-1]
    dy[0]  = dy[1]
    dy[ny] = dy[ny-1]

    for k in range(1, nz): dz[k] = z[k] - z[k-1]
    dz[0]  = dz[1]
    dz[nz] = dz[nz-1]

    bvx = xdict["bvalue"].replace(" ", "").split(',')
    bvy = ydict["bvalue"].replace(" ", "").split(',')
    bvz = zdict["bvalue"].replace(" ", "").split(',')

#   set Dirichet bcs as these take precedence
    for i in range(0, nx, nx-1):
       ib = int(i/(nx-1))
       if bcx[ib] == 'D':
          ax = dx[i]
          for j in range(0, ny):
             ay = 0.5*(dy[j] + dy[j+1])
             for k in range(0, nz):
                az = 0.5*(dz[k] + dz[k+1])
                m = mijk(i, j, k, nx, ny)
                if a[m][m] == 0:
                   a[m][m] = 2*ay*az/ax + ax*az/dy[j] + ax*az/dy[j+1] + ax*ay/dz[k] + ax*ay/dz[k+1] 
                   b[m]    = float(bvx[ib])*a[m][m]

    for j in range(0, ny, ny-1):
       jb = int(j/(ny-1))
       if bcy[jb] == 'D':
          ay = dy[j]
          for i in range(0, nx):
             ax = 0.5*(dx[i] + dx[i+1])
             for k in range(0, nz):
                az = 0.5*(dz[k] + dz[k+1])
                m = mijk(i, j, k, nx, ny)
                if a[m][m] == 0:
                   a[m][m] = 2*ax*az/ay + ay*az/dx[i] + ay*az/dx[i+1] + ax*ay/dz[k] + ax*ay/dz[k+1]
                   b[m]    = float(bvy[jb])*a[m][m]

    for k in range(0, nz, nz-1):
       kb = int(k/(nz-1))
       if bcz[kb] == 'D':
          az = dz[k]
          for i in range(0, nx):
             ax = 0.5*(dx[i] + dx[i+1])
             for j in range(0, ny):
                ay = 0.5*(dy[j] + dy[j+1])
                m = mijk(i, j, k, nx, ny)
                if a[m][m] == 0:
                   a[m][m] = 2*ax*ay/az + ay*az/dx[i] + ay*az/dx[i+1] + ax*az/dy[j] + ax*az/dy[j+1]
                   b[m]    = float(bvz[kb])*a[m][m]

#   set Neumann/Symmetry bcs next as these take precedence over repeating bcs
    for i in range(0, nx, nx-1):
       ib = int(i/(nx-1))
       ia = 1-2*ib
       if bcx[ib] == 'N' or bcx[ib] =='S':
          ax = dx[i]
          for j in range(0, ny):
             ay = 0.5*(dy[j] + dy[j+1])
             for k in range(0, nz):
                az = 0.5*(dz[k] + dz[k+1])
                m = mijk(i, j, k, nx, ny)
                if a[m][m] == 0:
                   a[m][m] = 2*ay*az/ax + ax*az/dy[j] + ax*az/dy[j+1] + ax*ay/dz[k] + ax*ay/dz[k+1]
                   a[m][mijk(i+ia, j, k, nx, ny)] = -a[m][m]
                   if bcx[ib] == 'N': b[m] = float(bvx[ib])*a[m][m]

    for j in range(0, ny, ny-1):
       jb = int(j/(ny-1))
       ja = 1-2*jb
       if bcy[jb] == 'N' or bcy[jb] =='S':
          ay = dy[j]
          for i in range(0, nx):
             ax = 0.5*(dx[i] + dx[i+1])
             for k in range(0, nz):
                az = 0.5*(dz[k] + dz[k+1])
                m = mijk(i, j, k, nx, ny)
                if a[m][m] == 0:
                   a[m][m] = 2*ax*az/ay + ay*az/dx[i] + ay*az/dx[i+1] + ax*ay/dz[k] + ax*ay/dz[k+1]
                   a[m][mijk(i, j+ja, k, nx, ny)] = -a[m][m]
                   if bcy[jb] == 'N': b[m] = float(bvy[jb])*a[m][m]

    for k in range(0, nz, nz-1):
       kb = int(k/(nz-1))
       ka = 1-2*kb
       if bcz[kb] == 'N' or bcz[kb] =='S':
          az = dz[k]
          for i in range(0, nx):
             ax = 0.5*(dx[i] + dx[i+1])
             for j in range(0, ny):
                ay = 0.5*(dy[j] + dy[j+1])
                m = mijk(i, j, k, nx, ny)
                if a[m][m] == 0:
                   a[m][m] = 2*ax*ay/az + ay*az/dx[i] + ay*az/dx[i+1] + ax*az/dy[j] + ax*az/dy[j+1]
                   a[m][mijk(i, j, k+ka, nx, ny)] = -a[m][m]
                   if bcy[jb] == 'N': b[m] = float(bvz[kb])*a[m][m]

#   set dx, dy, dz for repeating bcs
    if repeatx:
       dx[0]  = dx[nx-1]
       dx[nx] = dx[1]
    if repeaty:
       dy[0]  = dy[ny-1]
       dy[ny] = dy[1]
    if repeatz:
       dz[0]  = dz[nz-1]
       dz[nz] = dz[1]

#   set interior matrix entries
    for i in range(0, nx):
       ax = 0.5*(dx[i] + dx[i+1])
       for j in range(0, ny):
          ay = 0.5*(dy[j] + dy[j+1])
          for k in range(0, nz):
             m = mijk(i, j, k, nx, ny)
             if a[m][m] == 0:
                az = 0.5*(dz[k] + dz[k+1])
                me = mijk(i+1, j,   k,   nx, ny)
                mw = mijk(i-1, j,   k,   nx, ny)
                mn = mijk(i,   j+1, k,   nx, ny)
                ms = mijk(i,   j-1, k,   nx, ny)
                mu = mijk(i,   j,   k+1, nx, ny)
                md = mijk(i,   j,   k-1, nx, ny)

                if repeatx:
                   if i == 0:     mw = mijk(nx-1, j, k, nx, ny)
                   if i == nx-1:  me = mijk(0,    j, k, nx, ny)
                else:
                   if i == 0:     mw = m                        # should be redundant trap as only repeats unset
                   if i == nx-1:  me = m
                if repeaty:
                   if j == 0:     ms = mijk(i, ny-1, k, nx, ny)
                   if j == ny-1:  mn = mijk(i, 0,    k, nx, ny)
                else:
                   if j == 0:     ms = m                        # should be redundant trap as only repeats unset
                   if j == ny-1:  mn = m
                if repeatz:
                   if k == 0:     md = mijk(i, j, nz-1, nx, ny)
                   if k == nz-1:  mu = mijk(i, j,    0, nx, ny)
                else:
                   if k == 0:     md = m                        # should be redundant trap as only repeats unset
                   if k == nz-1:  mu = m

                a[m][mw] = -ay*az/dx[i]
                a[m][me] = -ay*az/dx[i+1]
                a[m][ms] = -ax*az/dy[j]
                a[m][mn] = -ax*az/dy[j+1]
                a[m][md] = -ax*ay/dz[k]
                a[m][mu] = -ax*ay/dz[k+1]
                a[m][m]  = -a[m][me] - a[m][mw] - a[m][ms] - a[m][mn] - a[m][md] - a[m][mu]
                b[m]     = f*ax*ay*az

#   if degen is off, fix row mijk(i,j,k, nx, ny)  if matrix is degenerate
    if not degen:
       if bcx[0] != 'D' and bcx[1] != 'D':
          i = int(xdict["degfix"])
          if bcy[0] != 'D' and bcy[1] != 'D':
             j = int(ydict["degfix"])
             if bcz[0] != 'D' and bcz[1] != 'D':
                k = int(zdict["degfix"])
                m = mijk(i, j, k, nx, ny)
                me = mijk(i+1, j,   k,   nx, ny)
                mw = mijk(i-1, j,   k,   nx, ny)
                mn = mijk(i,   j+1, k,   nx, ny)
                ms = mijk(i,   j-1, k,   nx, ny)
                mu = mijk(i,   j,   k+1, nx, ny)
                md = mijk(i,   j,   k-1, nx, ny)
                a[m][mw] = 0.0
                a[m][me] = 0.0
                a[m][ms] = 0.0
                a[m][mn] = 0.0
                a[m][md] = 0.0
                a[m][mu] = 0.0
                b[m]     = b[m]*a[m][m]

#   scale to give ||a|| = 1.0 in max norm
    amax = np.amax(a)
    a = a/amax
    b = b/amax

#   debug print
#   print(np.get_printoptions())
    with np.printoptions(precision=3, suppress=True, linewidth=100):
       print(a)

    return a, b

