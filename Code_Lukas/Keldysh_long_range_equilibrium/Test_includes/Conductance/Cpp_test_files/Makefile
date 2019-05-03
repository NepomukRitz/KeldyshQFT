  SHELL = /bin/bash

  CC    = icpc #vtcxx #g++ -3.4 # gcc g++ g++-3.4
#  CC    = g++ #vtcxx #g++ -3.4 # gcc g++ g++-3.4
  
#OPT   =  -qopenmp -g -Wall 
#  OPT   = -qopenmp -g -Wall
#  OPT   = -g -qopenmp -Wall -O0 #maybe use O2
  OPT   =  -qopenmp -Wall -O3 #maybe use O2

  I = -I. -I$(HOME)/include -I$(HOME)/Code/Keldysh_long_range_equilibrium/Includes -I$(MATLAB)/extern/include $(GSL_INC) 
#  I = -I. -I$(HOME)/include -I$(MATLAB)/extern/include $(ACML_INC) $(MKL_INC)
  L = -L. -L$(HOME)/lib64 -L$(MATLAB)/bin/$(ARCH) $(GSL_SHLIB)
  LINK = -Wl,-rpath=$(HOME)/lib64
#  L = -L. -L$(HOME)/lib64 -L$(MATLAB)/bin/$(ARCH) $(ACML_SHLIB) $(MKL_SHLIB)

#  LIB  = -lmx -lm -lmat -lmex -lmyppmd -lmwlapack -lmwblas $(GSL_BLAS_SHLIB)
  LIB  = -lmx -lm -lmat -lmex -lmyppmd -lmwlapack -lmwblas $(GSL_BLAS_SHLIB)
#  LIB  = -lmx -lm -lmat -lmex -lmypp -lmwlapack -lmwblas $(GSL_BLAS_SHLIB)

%:  %.cpp
	$(CC) $(OPT) -o $@ $@.cpp $(I) $(L) $(LIB) $(LINK)
	mv $@ ~/bin
    

# -------------------------------------------------------------------------- #

