
  SHELL = /bin/bash

  CC     = icpc #g++ #g++-3.4 # gcc g++ g++-3.4
  
  OPT = -O3 -openmp

ifndef $(LIBDIR)
  LIBDIR = $(HOME)/lib64
endif

  I = -I. -I$(HOME)/include -I$(MATLAB)/extern/include
  L = -L. -L$(LIBDIR) -L$(MATLAB)/bin/$(ARCH) $(GSL_SHLIB) $(GSL_BLAS_SHLIB) 

  LIB  = -lmx -lmex -lmyppmd -lmwlapack -lmwblas #-larprec # $(MKL_SHLIB)

ifndef $(MEXDIR)
#  MEXDIR = $(HOME)/Documents/MATLAB/MEX
endif

ifeq ($(ARCH),maci64)
  EXTEN = mexmaci64
endif

ifeq ($(ARCH),maci)
  EXTEN = mexmaci
endif

ifeq ($(ARCH),glnxa64)
  EXTEN = mexa64
endif

ifeq ($(ARCH),glnx86)
  EXTEN = mexglx
endif

#%: %.cpp %.h
%: %.cpp
	${CC} $(OPT) -fPIC -shared -fexceptions -o $@.$(EXTEN) $@.cpp $(I) $(L) $(LIB)
	mv $@.$(EXTEN) $(MEXDIR)
