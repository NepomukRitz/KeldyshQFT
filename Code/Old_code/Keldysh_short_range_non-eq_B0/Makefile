
  SHELL = /bin/bash

  CC    = icpc #vtcxx #g++ -3.4 # gcc g++ g++-3.4
  
  OPT   = -qopenmp -O3 

  I = -I. -I$(HOME)/include -I$(MATLAB)/extern/include $(GSL_INC)
#  I = -I. -I$(HOME)/include -I$(MATLAB)/extern/include $(ACML_INC) $(MKL_INC)
  L = -L. -L$(HOME)/lib64 -L$(MATLAB)/bin/$(ARCH) $(GSL_SHLIB)
#  L = -L. -L$(HOME)/lib64 -L$(MATLAB)/bin/$(ARCH) $(ACML_SHLIB) $(MKL_SHLIB)

#  LIB  = -lmx -lm -lmat -lmex -lmyppmd -lmwlapack -lmwblas $(GSL_BLAS_SHLIB)
  LIB  = -lmx -lm -lmat -lmex -lmypp -lmwlapack -lmwblas $(GSL_BLAS_SHLIB)
  LINK = -Wl,-rpath=$(HOME)/lib64

dfRGK_lin:  dfRG2_general_B0_clean.cpp
	$(CC) $(OPT) -DUSE_STEADY_STATE_APPROXIMATE_FDT -DLINEAR -o $@ dfRG2_general_B0_clean.cpp $(I) $(L) $(LIB)
	cp $@ ~/DATA/TEST/
	mv $@ ~/bin/

%:  %.cpp
	$(CC) $(OPT) -o $@ $@.cpp $(I) $(L) $(LIB) $(LINK)
	#cp $@ ~/DATA/TEST/
	mv $@ ~/bin/
    
%.o:  %.c $(DEPS)
	$(CC) $(OPT) -c -o $@ $< ${I}

# -------------------------------------------------------------------------- #

