===============
Post processing
===============

The following functions can be invoked after a calculation has been performed.
They are not required for the actual calculations themselves but are useful to extract additional information from
their results, either as consistency checks or to infer derived quantities for later analysis.

.. doxygenfunction:: check_SE_causality(const SelfEnergy<Q>& selfEnergy)

.. doxygenfunction:: check_FDTs_selfenergy(const SelfEnergy<Q>& selfenergy, const double T, const bool verbose)

.. doxygenfunction:: sum_rule_spectrum(const State<state_datatype> &state)

.. doxygenfunction:: compute_postprocessed_susceptibilities(const std::string& filename)

.. doxygenfunction:: save_slices_through_fullvertex(const std::string& filename, const int ispin)