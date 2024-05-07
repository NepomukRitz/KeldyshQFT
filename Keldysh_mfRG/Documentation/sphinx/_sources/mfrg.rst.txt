mfRG
====

.. doxygenfunction:: n_loop_flow(const std::string& outputFileName, const fRG_config &config)

.. doxygenfunction:: n_loop_flow(const std::string& inputFileName, const fRG_config& config, unsigned int it_start)

.. doxygenfunction:: ode_solver(Y& result, const Y& state_ini, const System& rhs, const ODE_solver_config& config=ODE_solver_config(), const bool verbose=true)

.. doxygenfunction:: rhs_n_loop_flow(const State<Q>& Psi, const double Lambda, const int nloops_max, const vec<size_t> opt, const fRG_config& config, mfRG_stats& stats)

.. doxygenfunction:: selfEnergyOneLoopFlow(SelfEnergy<Q>& dPsiSelfEnergy, const Vertex<Q,false>& PsiVertex, const Propagator<Q>& S)

.. doxygenfunction:: vertexOneLoopFlow(Vertex<Q,true>& dPsiVertex, const Vertex<Q,false>& PsiVertex, const Bubble_Object& dPi, const fRG_config& config)

.. doxygenfunction:: calculate_dGammaL(const GeneralVertex<Q, symmetric_r_irred,true>& dPsiVertex, const Vertex<Q,false>& PsiVertex, const Bubble_Object& Pi, const fRG_config& config)

.. doxygenfunction:: calculate_dGammaR(const GeneralVertex<Q, symmetric_r_irred,true>& dPsiVertex, const Vertex<Q,false>& PsiVertex, const Bubble_Object& Pi, const fRG_config& config)

.. doxygenfunction:: calculate_dGammaC_right_insertion(const Vertex<Q,false>& PsiVertex, const GeneralVertex<Q, non_symmetric_diffleft,true>& nonsymVertex, const Bubble_Object& Pi, const fRG_config& config)

.. doxygenfunction:: calculate_dGammaC_left_insertion(const GeneralVertex<Q, non_symmetric_diffright,true>& nonsymVertex, const Vertex<Q,false>& PsiVertex, const Bubble_Object& Pi, const fRG_config& config)

.. doxygenfunction:: selfEnergyFlowCorrections(SelfEnergy<Q>& dPsiSelfEnergy, const GeneralVertex<Q,symmetric_r_irred,true>& dGammaC_tbar, const State<Q>& Psi, const Propagator<Q>& G)

.. doxygenfunction:: postRKstep_stuff(Y& y_run, System& rhs, double x_run, const vec<double>& x_vals, int iteration, const std::string& filename, const ODE_solver_config& config, bool verbose)

