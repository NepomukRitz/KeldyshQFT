
  SHELL = /bin/bash

  CC    = mpicxx #vtcxx #g++ -3.4 # gcc g++ g++-3.4
# CC    = g++ #vtcxx #g++ -3.4 # gcc g++ g++-3.4
  
# OPT   =  -qopenmp -g -Wall -O0 -static-intel
# OPT   =  -qopenmp -O3 -static-intel
#  OPT   = -qopenmp   -O2 -Qoption,cpp,--extended_float_types -std=c++11 #-O3
  OPT   = -qopenmp -Wall -g -O0 -Qoption,cpp,--extended_float_types -std=c++11 #-O3
# OPT   = -fopenmp -Wall -g -O3  -std=c++11 #-O3
#  OPT   = -qopenmp   -O2 -g -Wall -std=c++11 #-O3
# OPT   = -g  -qopenmp -O0 -static-intel
#  OPT   = -qopenmp -O3
#  OPT   = -fopenmp -g  -O0 #-O3

PATH_TO_LIB=~/fRG_development/lib_without_matlab/

  I = -I. -I$(HOME)/fRG_development/lib_without_matlab -I$(HOME)/fRG_development/Keldysh_long_range_equilibrium/Includes -I$(HOME)/PAID_repository/PAID/include -I$EBROOTGSL/include
#  I = -I. -I$(HOME)/include -I$(MATLAB)/extern/include $(ACML_INC) $(MKL_INC)
  L = -L. -L$(HOME)/lib64  -L$EBROOTGSL/lib
#  LINK = -Wl,-rpath-link
#  L = -L. -L$(HOME)/lib64 -L$(MATLAB)/bin/$(ARCH) $(ACML_SHLIB) $(MKL_SHLIB)

#  LIB  = -lmx -lm -lmat -lmex -lmyppmd -lmwlapack -lmwblas $(GSL_BLAS_SHLIB)
  LIB  =  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lgsl -lgslcblas
%:  %.cpp
	$(CC) $(OPT) -o $@ $@.cpp $(PATH_TO_LIB)matrix.cpp $(PATH_TO_LIB)blnla.cpp $(PATH_TO_LIB)physicalpp.cpp $(PATH_TO_LIB)approxxpp.cpp $(PATH_TO_LIB)ozaki.cpp $(PATH_TO_LIB)patterson_set.cpp $(PATH_TO_LIB)basic.cpp $(I) $(L) $(LIB) 
	#mv $@ ~/fRG/HPC_Juelich/Keldysh_long_range_equilibrium/Main_program_mpi/Compiled
	mv $@ ~/bin
    

# -------------------------------------------------------------------------- #

