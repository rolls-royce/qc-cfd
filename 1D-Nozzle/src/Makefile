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

BIN_DIR  = ../bin
LIB_DIR  = ../lib/
PROGRAM  = $(BIN_DIR)/nozzle
LIBRARY  = $(LIB_DIR)libnozzle.so

C_SRC    = couple.c \
           error.c \
           ezxml.c \
           history.c \
           io.c \
           matsol.c \
           nozzle.c \
           simple.c \
           sparse.c

BIN_SRC  = $(C_SRC) main.c
LIB_SRC  = $(C_SRC) py_interface.c

CC       = cc
MKDIR    = mkdir -p
CFLAGS   = -O2 -fPIC
INCS     = -I../include
LIBS     = -lm

BIN_OBJS = $(BIN_SRC:.c=.o)
LIB_OBJS = $(LIB_SRC:.c=.o)

prog: $(PROGRAM)
lib:  $(LIBRARY)
all:  prog lib

$(PROGRAM): $(BIN_OBJS)
	$(MKDIR) $(BIN_DIR)
	$(CC) $(BIN_OBJS) $(INCS) $(LIBS) -o $@

$(LIBRARY): $(LIB_OBJS)
	$(MKDIR) $(LIB_DIR)
	$(CC) $(LIB_OBJS) $(INCS) -shared -o $@

.c.o:
	$(CC) -c $(CFLAGS) $(INCS) $<

clean:
	rm -f *.o $(PROGRAM) $(LIBRARY)
