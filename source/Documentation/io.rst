===
I/O
===

For input and output purposes, the `HDF5 <https://www.hdfgroup.org/HDF5/>`_ file format  is used throughout. To organize the data for output,
the contents of a state to be stored are split up into different datasets that correspond, to all the asymptotic
classes of the vertex in each channel, the self-energy, the frequency grids used, as well as the most important
parameters of the calculation. The output file is then on a high level organized in terms of '':math:`\Lambda` layers'',
the idea being that each layer enables access to a different state stored inside the same file. This enables the output
of, e.g., the results of a full mfRG flow into a single file, where each '':math:`\Lambda` layer index'' corresponds to a
different value of the regulator. Alternatively, this structure can be used to store the results of all the iterations
done during a solution of the parquet equations separately. Of course, one can equally well use just a single layer to
store the end result of a computation, such as a converged solution of the parquet equations or the result of a PT2 computation.

.. doxygenfunction:: write_state_to_hdf(const H5std_string FILE_NAME, double Lambda_i, const int Lambda_size, const State<Q,diff>& state_in, const bool verbose=true, const bool is_converged=false)

.. doxygenfunction:: add_state_to_hdf(const H5std_string FILE_NAME, int Lambda_it, const State<Q,diff>& state_in, const bool is_converged=false, const bool verbose=true)

.. doxygenfunction:: read_state_from_hdf(const H5std_string& filename, int Lambda_it)