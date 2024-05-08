Parquet equations
=================

.. doxygenfunction:: run_parquet(const fRG_config &config, const std::vector<double> &U_NRG_list, const int version, const bool overwrite_old_results)

.. doxygenfunction:: parquet_solver(const std::string filename, State<Q>& state_in, const double Lambda, const int version, const double accuracy=1e-6, const int Nmax=6, const bool overwrite_old_results=true, const double mixing_ratio=1.0, const bool use_last_state_anyway=false)

.. doxygenfunction:: parquet_iteration(State<Q>& state_out, const State<Q>& state_in, const double Lambda, const int it_Lambda, const int version)

.. doxygenfunction:: compute_BSE(Vertex<Q,false>& Gamma_BSE, Vertex<Q,false>& Gamma_BSE_L, Vertex<Q,false>& Gamma_BSE_R, const State<Q>& state_in, const double Lambda)

.. doxygenfunction:: compute_SDE(SelfEnergy<Q>& Sigma_SDE, const State<Q>& state_in, const double Lambda, int version)