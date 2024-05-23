.. KeldyshDiag documentation master file, created by
   sphinx-quickstart on Wed Feb 14 16:36:10 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

====================
KeldyshDiagrammatics
====================
**A C++ Codebase for real-frequency multiloop functional renormalization group and parquet computations for the single impurity Anderson model**


This is the technical documentation for the computational
framework designed for the accurate calculation of real-frequency
dynamical correlation functions of the single-impurity Anderson model
(AM) in the regime of weak to intermediate coupling, see its `source code <https://github.com/NepomukRitz/KeldyshDiagrammatics/tree/master_to_be_published>`_.

It uses quantum field theory within the Keldysh formalism to directly access the
self-energy and dynamical susceptibilities along the real-frequency
axis, as detailed in `this publication <https://doi.org/10.1103/PhysRevB.109.115128>`_. The primary challenge we
address is the computationally demanding task of obtaining the full
three-dimensional real-frequency dependence of the four-point vertex.

This codebase provides a fully MPI+OpenMP parallelized implementation of the multiloop
functional renormalization group, and, independently, the solution
of the self-consistent parquet equations within the parquet
approximation. It leverages vectorization to handle the additional
complexity imposed by the Keldysh formalism, using optimized data
structures and highly performant integration routines.


The code itself is written in `C++17 <https://en.cppreference.com/w/cpp/17>`_
and is built using `CMake <https://cmake.org/documentation/>`_, demanding at least version 3.10. It requires the
`GSL <https://www.gnu.org/software/gsl/doc/latex/gsl-ref.pdf>`_, `boost <https://www.boost.org/>`_ and
`Eigen3 <http://eigen.tuxfamily.org/>`_ libraries as well as the `HDF5 <https://www.hdfgroup.org/HDF5/>`_ library for
input and output. For parallelization, the `OpenMP <https://www.openmp.org/>`_ and
`MPI <https://hpc.nmsu.edu/discovery/mpi/introduction/>`_ interfaces are used.
All these libraries are hence required on the system that is used to build and run the code.
Notably, we do not supply precompiled executables, that could be run directly, for several reasons:
First, the code makes heavy use of preprocessor flags that must be set *before* compilation and that are
in part used to specify the concrete problem at hand, see :doc:`parameters`. Second, depending on the system
at hand, special compilers for the particular architecture might be available, which could optimize the code during
compilation and linking, leading to a significant improvement in performance. The user should hence adapt the file
`CMakeLists.txt <https://github.com/NepomukRitz/KeldyshDiagrammatics/blob/master_to_be_published/source/CMakeLists.txt>`_
accordingly, such that the required libraries are included and linked properly and all compiler settings are set as desired.

.. toctree::
   :maxdepth: 1
   :caption: Contents

   correlationfunctions
   diagrams
   postprocessing
   io
   miscellaneous
   parameters
   algorithms

This technical documentation supplied with the code is generated automatically using the tools
`Doxygen <https://www.doxygen.nl/>`_, `Sphinx <https://www.sphinx-doc.org/>`_,
`Breathe <https://breathe.readthedocs.io/en/latest/>`_ and `CMake <https://cmake.org/documentation/>`_.

----------------

* :ref:`genindex`
