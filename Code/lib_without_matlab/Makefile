
  SHELL = /bin/bash

  CC     = icpc
  

  OPT = -g #$(OPT2) $(OPT1) $(OPT0)

  I = -I.  $(GSL_INC) $(MKL_INC)
  
  LIB  = $(MKL_SHLIB) $(GSL_SHLIB) $(GSL_BLAS_SHLIB)

 .SUFFIXES: .c .o .a .cpp .cc .h .cc

 objects =  matrix.o blnla.o physicalpp.o approxxpp.o
 objects2 = ozaki.o patterson_set.o basic.o
default: c


dlib:
	make -f Makefile libmyppd.so
	make -f Makefile clean
    
lib:
	make  -f Makefile libmypp.so
	make  -f Makefile clean
	make  -f Makefile libmyppd.so
	make  -f Makefile clean

fast:
	make  -f Makefile libmypp.so

libmypp.so: $(objects)  $(objects2)
	$(CC) $(OPT) -shared -o $@ $(objects) $(objects2) $(L) $(I) $(LIB)
	cp *.h $(HOME)/include_without_matlab
	mv $@ $(HOME)/lib64_without_matlab
   
libmyppd.so: DEBUGFLAG=-DDEBUG
libmyppd.so: $(objects) $(objects2)
	$(CC) $(OPT) -shared -o $@ $(objects) $(objects2) $(L) $(I) $(LIB)
	cp *.h $(HOME)/include_without_matlab
	mv $@ $(HOME)/lib64_without_matlab

%.o:  %.cpp %.h
	$(CC) $(OPT) -shared -fPIC $(DEBUGFLAG) -c -o $@ $< $(I) $(LIB)

clean:
	rm $(objects)

# -------------------------------------------------------------------------- #

