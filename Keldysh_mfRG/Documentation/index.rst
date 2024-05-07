.. KeldyshDiag documentation master file, created by
   sphinx-quickstart on Wed Feb 14 16:36:10 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

A C++ Codebase for real-frequency multiloop functional renormalization group and parquet computations for the single impurity Anderson model
============================================================================================================================================

TODO: Add introductory paragraph


.. toctree::
   :maxdepth: 2
   :caption: Contents

   prerequisites
   correlationfunctions
   diagrams
   postprocessing
   io
   miscellaneous
   parameters
   algorithms



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. [TODO: remove] Docs
.. ===================
.. 
.. .. doxygenstruct:: fRG_config
..    :members:
.. 
.. .. doxygenclass:: Hartree_Solver
..    :members:
.. 
.. .. doxygendefine:: KELDYSH_FORMALISM
.. 
.. .. doxygenfunction:: run_parquet(const fRG_config &config, const std::vector<double> &U_NRG_list, const int version, const bool overwrite_old_results)
.. 
.. .. doxygenfunction:: sopt_state_impl(State<Q>& Psi, const Bubble_Object& Pi, const State<Q>& bareState)
.. 
.. .. doxygennamespace:: utils