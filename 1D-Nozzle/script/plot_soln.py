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

import sys, getopt
from os.path import exists
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

###########################################################
#   get command line arguments                            #
###########################################################
def read_args():
    filelist = []
    argv = sys.argv[1:]

    for arg in argv:
       if arg == '-h':
          print('Usage:\n\t%s -h <file1> <file2> ...' % str(sys.argv[0]))
          print('\tA maximum of 5 solution files may be provided\n')
          sys.exit()
       else:
          if exists(arg):
             filelist.append(arg)
          else:
             print('Warning: file %s does not exist' % arg)

       if len(filelist) == 5:
          print('Warning: maximum number of files reached')
          break

#   check solution files are all of the same flow type
#   --------------------------------------------------
    ni = 0
    nc = 0
    for f in filelist:
       file_CSV = open(f)
       data_CSV = csv.reader(file_CSV)
       list_CSV = list(data_CSV)
       head = list_CSV[0]

       i = 0
       for h in head:
          head[i] = h.strip(' " " ')
          i += 1

       if 'Mach number' in head:
          nc += 1
       else:
          ni += 1

    if ni>0 and nc>0:
       print('\nError: cannot mix incompressible and compressible solutions\n')
       sys.exit(1)
    if nc>0:
       compressible = True
    else:
       compressible = False

    return compressible, filelist


###########################################################
#   read flow solution                                    #
###########################################################
def read_soln(filename):

#   read data
    file_CSV = open(filename)
    data_CSV = csv.reader(file_CSV)
    list_CSV = list(data_CSV)

#   get headers
    head = list_CSV[0]

#   strip spaces and quotation marks
    i = 0
    for h in head:
       head[i] = h.strip(' " " ')
       i += 1

#   get solution
    del list_CSV[0]
    soln = np.asarray(list_CSV, dtype=np.float64)

    return head, soln


###########################################################
#   plot incompressible flow                              #
###########################################################
def plot_incomp(filelist):

    fig, axs = plt.subplots(ncols=2, nrows=1)
    ax1, ax2 = axs.flat

    for ax in (ax1, ax2):
       ax.grid()

    colours = iter(['r', 'g', 'b', 'c', 'm'])

    ax1.title.set_text('velocity')
    ax2.title.set_text('pressure')

    for f in filelist:
       c = next(colours)
       head, soln = read_soln(f)
       ax1.plot(soln[:,1], soln[:,4], label = f, linestyle = 'solid',  color = c)
       ax2.plot(soln[:,1], soln[:,5], label = f, linestyle = 'solid',  color = c)

    ax1.legend(loc="lower center")
    ax2.legend(loc="upper center")

    fig.suptitle("Flow solution", fontsize="x-large")
    plt.gcf().set_size_inches(10, 5)               # pane dimensions in inches
    plt.show()


###########################################################
#   plot compressible flow                                #
###########################################################
def plot_comp(filelist):

    fig, axs = plt.subplots(ncols=2, nrows=2)
    ax1, ax2, ax3, ax4 = axs.flat

    for ax in (ax1, ax2, ax3, ax4):
       ax.grid()

    colours = iter(['r', 'g', 'b', 'c', 'm'])

    ax1.title.set_text('velocity')
    ax2.title.set_text('pressure')
    ax3.title.set_text('area')
    ax4.title.set_text('Mach number')

    for f in filelist:
       c = next(colours)
       head, soln = read_soln(f)
       ax1.plot(soln[:,1], soln[:,4], label = f, linestyle = 'solid',  color = c)
       ax2.plot(soln[:,1], soln[:,5], label = f, linestyle = 'solid',  color = c)
       ax3.plot(soln[:,1], soln[:,2], label = f, linestyle = 'solid',  color = c)
       ax4.plot(soln[:,1], soln[:,7], label = f, linestyle = 'solid',  color = c)

    ax1.legend(loc="lower center")
    ax2.legend(loc="upper center")
    ax3.legend(loc="upper center")
    ax4.legend(loc="lower center")

    fig.suptitle("Flow solution", fontsize="x-large")
    plt.gcf().set_size_inches(10, 8)               # pane dimensions in inches
    plt.show()


###########################################################
#   main routine                                          #
###########################################################
if __name__ == "__main__":

    compressible, filelist = read_args()

    if compressible:
       plot_comp(filelist)
    else:
       plot_incomp(filelist)

