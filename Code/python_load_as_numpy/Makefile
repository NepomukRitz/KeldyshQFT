  SHELL = /bin/bash

  CC    = icpc #vtcxx #g++ -3.4 # gcc g++ g++-3.4
#  CC    = g++ #vtcxx #g++ -3.4 # gcc g++ g++-3.4
  
 OPT   =  -shared -O3 `python3-config --cflags` 
#  OPT   = -qopenmp -O3
#  OPT   = -fopenmp -g  -O0 #-O3

PATH_TO_LIB=~/fRG_development/lib_without_matlab/

  #I = -I. -I$(HOME)/include -I$(HOME)/Aufgeraeumt/Includes -I$EBROOTGSL/include -I/usr/local/include/python2.7
  #I = -I. -I$(HOME)/include -I$(HOME)/Aufgeraeumt/Includes -I$EBROOTGSL/include -I/usr/local/include/python3.6
  I = -I. -I$(HOME)/fRG_development/lib_without_matlab/ -I$(HOME)/fRG_development/Keldysh_long_range_equilibrium/Includes -I$EBROOTGSL/include -I/usr/local/include/python3.6
#  I = -I. -I$(HOME)/include -I$(MATLAB)/extern/include $(ACML_INC) $(MKL_INC)
  #L = -L. -L$(HOME)/lib64  -L$EBROOTGSL/lib 
  #L = -L. -L$(HOME)/lib64  -L$EBROOTGSL/lib -L/usr/local/software/jureca/Stages/2017b/software/Boost/1.65.1-ipsmpi-2017b-Python-2.7.14/lib
  L = -L. -L$(HOME)/lib64  -L$EBROOTGSL/lib -L/gpfs/software/juwels/stages/2018a/software/SciPy-Stack/2018a-intel-para-2018a-Python-3.6.5/lib/
  #LINK = -Wl,-rpath=/usr/local/software/jureca/Stages/2017b/software/Boost/1.65.1-ipsmpi-2017b-Python-2.7.14/lib
  LINK = -Wl,-rpath=/gpfs/software/juwels/stages/2018a/software/SciPy-Stack/2018a-intel-para-2018a-Python-3.6.5/lib
#  L = -L. -L$(HOME)/lib64 -L$(MATLAB)/bin/$(ARCH) $(ACML_SHLIB) $(MKL_SHLIB)

#  LIB  = -lmx -lm -lmat -lmex -lmyppmd -lmwlapack -lmwblas $(GSL_BLAS_SHLIB)
  #LIB  =  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lgsl -lgslcblas -lboost_python3 -lpython2.7
  LIB  =  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lgsl -lgslcblas -lboost_python36 
%:  %.cpp
	$(CC) $(OPT) -o $@.so -fPIC $@.cpp  $(PATH_TO_LIB)matrix.cpp $(PATH_TO_LIB)blnla.cpp $(PATH_TO_LIB)physicalpp.cpp $(PATH_TO_LIB)approxxpp.cpp $(PATH_TO_LIB)ozaki.cpp $(PATH_TO_LIB)patterson_set.cpp $(PATH_TO_LIB)basic.cpp $(I) $(L) $(LIB) $(LINK)
	#$(CC) $(OPT) -o $@.so -fPIC $@.cpp $(I) $(L) $(LIB) 
	#mv $@ ~/Aufgeraeumt/Test_Includes/P_bubble_central_zero_mag/Compiled
    

# -------------------------------------------------------------------------- #

