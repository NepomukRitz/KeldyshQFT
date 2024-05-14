=============
Prerequisites
=============

TODO: Add links

The code itself is written in ``C++17``  and is built using ``CMake``, demanding at least version 3.10. It requires the ``GSL``, ``boost`` and ``Eigen3`` libraries as well as the ``HDF5`` library for input and output. For parallelization, the ``OpenMP`` and ``MPI`` interfaces are used. All these libraries are hence required on the system that is used to build and run the code. Notably, we do not supply precompiled executables, that could be run directly, for several reasons: First, the code makes heavy use of preprocessor flags that must be set \emph{before} compilation and that are in part used to specify the concrete problem at hand, see (TODO: link to parameters). Second, depending on the system at hand, special compilers for the particular architecture might be available, which could optimize the code during compilation and linking, leading to a significant improvement in performance. The user should hence adapt the file ``CMakeLists.txt`` accordingly, such that the required libraries are included and linked properly and all compiler settings are set as desired.

This technical documentation supplied with the code is generated automatically using the tools ``Doxygen``, ``Sphinx``, ``Breathe`` and ``CMake``.