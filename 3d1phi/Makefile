# Copyright (c) 1998 Lawrence Livermore National Security, LLC and other
# HYPRE Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

########################################################################
# Compiler and external dependences
########################################################################
CC        = mpicxx
F77       = mpif77
CXX       = mpicxx
F90       = mpif90
CUF       = mpixlcuf
LINK_CC   = ${CXX}
LINK_CXX  = ${CXX}
LINK_FC   = ${F77}
LINK_CUF  = ${CUF}

XL_DIR=$(dir $(shell which xlc))..

HYPRE_DIR = ../hypre

########################################################################
# CUDA
########################################################################
ifeq ($(use_cuda), 1)
   CUDA_INCL = -I${CUDA_HOME}/include
   CUDA_LIBS = -L${CUDA_HOME}/lib64 -lcudart -lcusparse -lcurand -lstdc++ -L$(XL_DIR)/xlC/16.1.1/lib -libmc++
   CUDA_ARCH = -gencode arch=compute_60,code=sm_60 -gencode arch=compute_70,code=sm_70
   NVCC_LDFLAGS = -ccbin=${CXX} ${CUDA_ARCH}
   COPTS_CUDA = -DHYPRE_EXAMPLE_USING_CUDA
   FOPTS_CUDA = -DHYPRE_EXAMPLE_USING_CUDA -qsuppress=cmpmsg
endif

ifeq ($(use_cuf), 1)
   CUDA_LIBS = -L${CUDA_HOME}/lib64 -lcudart -lcusparse -lcurand -lstdc++ -L$(XL_DIR)/xlC/16.1.1/lib -libmc++
endif

########################################################################
# Device OMP
########################################################################
ifeq ($(use_domp), 1)
   DOMP_LIBS = -L${CUDA_HOME}/lib64 -lcudart -lcusparse -lcurand -lstdc++ -L$(XL_DIR)/xlC/16.1.1/lib -libmc++
   COPTS_DOMP = -DHYPRE_EXAMPLE_USING_DEVICE_OMP
   FOPTS_DOMP = -qoffload -W@,-v -qsmp=omp -qinfo=omperrtrace -DHYPRE_EXAMPLE_USING_DEVICE_OMP
   LOPTS_DOMP = -qoffload -W@,-v -qsmp=omp
endif

########################################################################
# Compiling and linking options
########################################################################
CINCLUDES = -I$(HYPRE_DIR)/include $(CUDA_INCL)
#CDEFS = -DHYPRE_EXVIS
CDEFS =
COPTS = -O2  -g  $(COPTS_CUDA) $(COPTS_DOMP)
FOPTS = -g $(FOPTS_CUDA) $(FOPTS_DOMP)
CFLAGS = $(COPTS) $(CINCLUDES) $(CDEFS)
FINCLUDES = $(CINCLUDES)
FFLAGS = $(FOPTS) $(FINCLUDES)
CUFFLAGS = -qcuda
CXXOPTS = $(COPTS) -Wno-deprecated
CXXINCLUDES = $(CINCLUDES) -I..
CXXDEFS = $(CDEFS)
IFLAGS_BXX =
CXXFLAGS  = $(CXXOPTS) $(CXXINCLUDES) $(CXXDEFS) $(IFLAGS_BXX)
IF90FLAGS =
F90FLAGS = $(FFLAGS) $(IF90FLAGS)

LINKOPTS = $(LOPTS_CUDA) $(LOPTS_DOMP)
LIBS = -L$(HYPRE_DIR)/lib -lHYPRE -lm $(CUDA_LIBS) $(DOMP_LIBS)
LFLAGS = $(LINKOPTS) $(LIBS)
LFLAGS_B =\
 -L${HYPRE_DIR}/lib\
 -lbHYPREClient-C\
 -lbHYPREClient-CX\
 -lbHYPREClient-F\
 -lbHYPRE\
 -lsidl -ldl -lxml2
LFLAGS77 = $(LFLAGS)
LFLAGS90 =

########################################################################
# Rules for compiling the source files
########################################################################
.SUFFIXES: .c .f .cuf .cxx

%o:%c
	$(CC) $(CXXFLAGS) -c $<

########################################################################
# List of all programs to be compiled
########################################################################
ALLPROGS = ex5 \
           
BIGINTPROGS = ex5big ex15big
FORTRANPROGS = ex5f ex12f
MAXDIMPROGS = ex17 ex18
COMPLEXPROGS = ex18comp

all: $(ALLPROGS)

default: all

gpu: all

bigint: $(BIGINTPROGS)

fortran: $(FORTRANPROGS)

maxdim: $(MAXDIMPROGS)

complex: $(COMPLEXPROGS)
########################################################################
# Example 1
########################################################################
ex5: first_step.o poisson.o io.o io_2d2phc.o math_utils.o matrix_solver.o mem.o utils.o utils_2d2phc.o second_step.o third_step.o ex5.o calc_matrix.o first_step_rhs.o second_step_rhs.o
	$(LINK_CC) -O3  $^ -o $@ $(LFLAGS) -lmpi

########################################################################
# Clean up
########################################################################
clean:
	rm *.o
distclean: clean
	rm -f $(ALLPROGS) $(ALLPROGS:=*~)
	rm -f $(BIGINTPROGS) $(BIGINTPROGS:=*~)
	rm -f $(FORTRANLPROGS) $(FORTRANPROGS:=*~)
	rm -f $(MAXDIMPROGS) $(MAXDIMPROGS:=*~)
	rm -f $(COMPLEXPROGS) $(COMPLEXPROGS:=*~)
	rm -fr README*
