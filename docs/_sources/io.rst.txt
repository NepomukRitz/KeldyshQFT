===
I/O
===

.. doxygenfunction:: write_state_to_hdf(const H5std_string FILE_NAME, double Lambda_i, const int Lambda_size, const State<Q,diff>& state_in, const bool verbose=true, const bool is_converged=false)

.. doxygenfunction:: add_state_to_hdf(const H5std_string FILE_NAME, int Lambda_it, const State<Q,diff>& state_in, const bool is_converged=false, const bool verbose=true)

.. doxygenfunction:: read_state_from_hdf(const H5std_string& filename, int Lambda_it)