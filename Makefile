# Except as specified in the OpenAccess terms of use of Cadence or Silicon
# Integration Initiative, this material may not be copied, modified,
# re-published, uploaded, executed, or distributed in any way, in any medium,
# in whole or in part, without prior written permission from Cadence.
#
#                Copyright 2002-2005 Cadence Design Systems, Inc.
#                           All Rights Reserved.
#
#  $Author: pdsim $
#  $Revision: #21 $
#  $Date: 2006/11/15 $
# ******************************************************************************
# ******************************************************************************


#include ../macro.defs
all : DEFAULT
PROG = WERSim
#thread_lib =  /u/home/puneet/shaodiwa/MeRAM/switch_cuda/MeRam_Cuda_2/pthread_library
#thread_o = $(thread_lib)/common.o $(thread_lib)/thread.o
NVCC = /u/local/cuda/5.0/bin/nvcc
DEFAULT: $(PROG)
	#$(MAKE) -f $(thread_lib)/Makefile
CUPROG = LLG 

CCFLAG = -I/u/local/cuda/5.0/include/ -arch=sm_13 --optimize 2 #-G -g#/u/local/cuda/current/include/ #-G -g
CCPATH   = g++
 


# Compile the application code
objects = $(PROG).o 


$(PROG).o: $(PROG).cu LLG.cu NDR_Solver.cu
	$(NVCC) $(CCFLAG) -o $(PROG).o \
	 -c $(PROG).cu




# Link the executable

$(PROG): $(PROG).o 
	$(NVCC) $(CCFLAG)  -o $(PROG)  $(PROG).o
          



#$(PROG): $(PROG).o  $(OA_LIB_LIST)
#	$(CCPATH) $(DEBUG) $(CXXOPTS) -o $(PROG)  $(PROG).o \
#	 $(CCLNFLAG) 


clean: 
	@/bin/rm  -rf $(PROG).o 
	@/bin/rm  -rf $(PROG)
	@/bin/rm  -rf $(CUPROG).o
