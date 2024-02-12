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

import xml.etree.ElementTree as ET
import numpy as np

###########################################################
#   get mesh parameters from XML input file               #
###########################################################
def parse_meshfile(inputfile):

    print('\nreading input from XML file:', inputfile)
    root = ET.parse(inputfile).getroot()

    rdict = {}

    for child in root:
       if child.tag == 'case':
          casename = child.attrib["name"]
          ndims    = int(child.attrib["dimension"])
          force    = float(child.attrib["force"])
       elif child.tag == 'mesh':
          xyz = child.attrib["direction"]
          mdict = {}
          for m in child:
             mdict[m.tag] = m.text.strip()
          rdict[xyz] = mdict.copy()
          mdict.clear()

    rdict['force'] = force
    return casename, ndims, rdict

###########################################################
#   generate mesh in single coordinate direction          #
###########################################################
def generate_mesh(mdict):
    L  = int(mdict["length"])
    r  = float(mdict["cratio"])
    nt = int(mdict["ntotal"])
    nc = int(mdict["nclust"])
    ct = int(mdict["cltype"])

    if ct == 2:
       fc = 2
    else:
       fc = 1

#   solve clustering equations
    nu = nt - fc*(nc-1)
    if r == 1.0:
       d = L/(nt-1)
       D = d
       C = (nc-1)*d
    else:
       rc = pow(r, nc-1)
       f1 = (nu-1)*rc
       f2 = fc*(rc-1)/(r-1)
       d  = L/(f1+f2)
       D  = rc*d
       C  = (L - (nu-1)*D)/fc

#   debug
    if False:
       if ct == 2:
          print("nu, d, D, C, L = %d %f, %f, %f %f" % (nu, d, D, C, (nu-1)*D+2*C))
       else:
          print("nu, d, D, C, L = %d %f, %f, %f %f" % (nu, d, D, C, (nu-1)*D+C))

#   solve for coordinates
    x = np.zeros(shape=(nt))

    x[0] = 0.0
    for i in range(1, nc-1):
      x[i] = x[i-1] + d*pow(r,i-1)

    x[nc-1] = C
    for i in range(0, nu-1):
      x[nc+i] = x[nc+i-1] + D

    ns = nc+nu-1
    if ct == 2:
      for i in range(0, nc-1):
        x[ns+i] = x[ns+i-1] + d*pow(r, nc-i-2)

#   flip if cluster type = -1
    if ct == -1:
       return L - np.flip(x)
    else:
       return x

