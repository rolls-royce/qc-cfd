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

###########################################################
#   get command line arguments                            #
###########################################################
def read_args():
    filelist = []
    argv = sys.argv[1:]

    for arg in argv:
       if arg == '-h':
          print('Usage:\n\t%s -h <file1> <file2> ...' % str(sys.argv[0]))
          print('\tA maximum of 5 history files may be provided\n')
          sys.exit()
       else:
          if exists(arg):
             filelist.append(arg)
          else:
             print('Warning: file %s does not exist' % arg)

       if len(filelist) == 5:
          print('Warning: maximum number of files reached')
          break

    return filelist


###########################################################
#   read convergence history                              #
###########################################################
def read_hist(filename):

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

#   get history
    del list_CSV[0]
    hist = np.asarray(list_CSV, dtype=np.float64)

    return head, hist


###########################################################
#   main routine                                          #
###########################################################
if __name__ == "__main__":

    filelist = read_args()

    fig, axs = plt.subplots(ncols=2, nrows=2)
    ax1, ax2, ax3, ax4 = axs.flat

    for ax in (ax1, ax2, ax3, ax4):
       ax.semilogy()
       ax.grid()
    
    colours = iter(['r', 'g', 'b', 'c', 'm'])

    ax1.title.set_text('du_rms')
    ax2.title.set_text('dp_rms')
    ax3.title.set_text('resu_rms')
    ax4.title.set_text('resp_rms')

    for f in filelist:
       c = next(colours)
       head, hist = read_hist(f)
       ax1.plot(hist[:,1], label = f, linestyle = 'solid',  color = c)
       ax2.plot(hist[:,2], label = f, linestyle = 'solid',  color = c)
       if 'dr_rms' in head:                        # compressible history
          ax3.plot(hist[:,4], label = f, linestyle = 'solid',  color = c)
          ax4.plot(hist[:,5], label = f, linestyle = 'solid',  color = c)
       else:
          ax3.plot(hist[:,3], label = f, linestyle = 'solid',  color = c)
          ax4.plot(hist[:,4], label = f, linestyle = 'solid',  color = c)

    ax1.legend(loc="upper right")
    ax2.legend(loc="upper right")
    ax3.legend(loc="upper right")
    ax4.legend(loc="upper right")

    fig.suptitle("Convergence histories", fontsize="x-large")
    plt.gcf().set_size_inches(10, 8)               # pane dimensions in inches
    plt.show()

