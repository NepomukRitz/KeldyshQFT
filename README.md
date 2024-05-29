# KeldyshQFT

**A C++ Codebase for real-frequency multiloop functional renormalization group and parquet computations for the 
single  impurity Anderson model**

This is the source code designed for the accurate calculation of real-frequency dynamical correlation functions of 
the  single-impurity Anderson model (AM) in the regime of weak to intermediate coupling as detailed in [this 
publication](https://doi.org/10.1103/PhysRevB.109.115128).

It provides a fully MPI+OpenMP parallelized implementation of the multiloop functional renormalization
group, and, independently, the solution of the self-consistent parquet equations within the parquet approximation.

## Prerequisites & Installation

The code is written in [`C++17`](https://en.cppreference.com/w/cpp/17)
and is built using [`CMake`](https://cmake.org/documentation/), demanding at least version 3.10.

To build the code, adapt the file 
[`CMakeLists.txt`](https://github.com/NepomukRitz/KeldyshDiagrammatics/blob/master_to_be_published/source/CMakeLists.txt) from  the `source/` directory to specify all compiler and linker settings and to  properly include the required
[`GSL`](https://www.gnu.org/software/gsl/doc/latex/gsl-ref.pdf), [`boost`](https://www.boost.org/),
[`Eigen3`](http://eigen.tuxfamily.org/) and [`HDF5`](https://www.hdfgroup.org/HDF5/) libraries.

Afterwards, specify all parameters and system settings in the header files included in the directory 
`source/parameters/`, 
see also the [documentation](https://nepomukritz.github.io/KeldyshQFT/parameters.html).

Finally, compile 
[`main.cpp`](https://github.com/NepomukRitz/KeldyshDiagrammatics/blob/master_to_be_published/source/main.cpp). 

## Documentation

Find the documentation [here](https://nepomukritz.github.io/KeldyshQFT/). Comments are welcome!

## License

Published under the [MIT license](LICENSE.md).

## How to cite

If you wish to cite this software, please use the metadata provided in the [citation file](CITATION.cff).