CXX=mpicxx
CFLAGS+= -std=c++14 -O3 -lhdf5 -lhdf5_cpp -fopenmp  
INCS+=-I/home/j/J.Thoenniss/Documents/hdf5-1.10.4/hdf5/include -I/software/opt/xenial/x86_64/openmpi/1.8.4/include
LIBS+=-L/home/j/J.Thoenniss/Documents/hdf5-1.10.4/hdf5/lib -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm

FLAGS+= 
%.o: %.cpp
	$(CXX) $(CFLAGS) $(INCS) $(LIBS) -c $<

all: main.o kagome1d.o hdf5_routines.o 
	$(CXX) -o main  main.o  kagome1d.o hdf5_routines.o $(CFLAGS) $(LIBS) $(INCS)
	@echo "compiled successfully"	
main: main.o
	$(CXX) -c  main.o $(CFLAGS) $(LIBS) $(INCS)
	@echo "---successfully compiled and linked test"

kagome: kagome1d.o
	$(CXX) -c kagome1d.o $(CFLAGS) $(LIBS) $(INCS)

hdf5: hdf5_routines.o
	$(CXX) -c hdf5_routines.o $(CFLAGS) $(LIBS) $(INCS)

.phony: clean
clean:
	rm -f *.o main  
