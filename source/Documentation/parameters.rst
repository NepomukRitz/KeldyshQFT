==========
Parameters
==========

Pre-processor flags
===================

- ``ADAPTIVE_GRID``
    | If defined, use optimization routine to find the best scale factor :math:`A` of the frequency grid;
    | If undefined, just rescale the grid.
    | **Warning**: Can be expensive and unreliable in the presence of numerical artifacts.

- ``ANALYTIC_TAILS``
    | 0 for false; 1 for true.
    | If true, the analytic expression for the bare bubble is used to treat the high-frequency asymptotics during bubble computations in the finite-:math:`T` MF.

- ``BARE_SE_FEEDBACK``
    | If defined, only bare selfenergy is used. Only makes sense if ``STATIC_FEEDBACK`` is defined. Useful for benchmarks with previous Keldysh fRG schemes.

- ``CONTOUR_BASIS``
    | 0 for false, 1 for true.
    | If true, no Keldysh rotation is performed and the contour basis is used instead to parametrize the Keldysh components of all correlation functions. Useful for comparisons with results that use this convention. Not as well tested and thus not recommended for production runs.

- ``DEBUG_SYMMETRIES``
    | 0 for false; 1 for true.
    | Performs computations without use of symmetries if true. Useful for debugging purposes.

- ``GRID``
    | Controls which frequency grid is to be used.
    | 0 for the non-linear grid, 1 for the hybrid grid, 2 for the polar grid.
    | **Recommendation**: 0.

- ``KATANIN``
    | If defined, the Katanin extension is used during fRG computations.

- ``KELDYSH_FORMALISM``
    | Determines whether calculations shall be done in the Keldysh or Matsubara formalism.
    | 0 for Matsubara formalism (MF); 1 for Keldysh formalism (KF).

- ``MAX_DIAG_CLASS``
    | Defines the diagrammatic classes that will be considered.
    | 1 for only :math:`K_1`, 2 for :math:`K_1` and :math:`K_2` and 3 for the full dependencies.
    | Useful for debugging purposes and for computations in second-order perturbation theory, or if ``STATIC_FEEDBACK`` is defined, when only :math:`K_1` is required.

- ``NDEBUG``
    | If defined, assert functions are switched off. Recommended setting for production runs.

- ``PARTICLE_HOLE_SYMM``
    | 0 for false; 1 for true.
    | If true, particle-hole symmetry is assumed.

- ``PT2_FLOW``
    | If defined, only compute the flow equations up to :math:`O(U^2)`. Only makes sense for pure :math:`K_1` calculations. Useful as a consistency check together with independent PT2 calculations.

- ``REG``
    | Specifies the mfRG flow regulator to be used.
    | :math:`2`: :math:`\Delta`-flow, :math:`3`: :math:`\omega`-flow, :math:`4`: :math:`U`-flow, :math:`5`: :math:`T`-flow.

- ``REPARAMETRIZE_FLOWGRID``
    | If defined, the flow parameter is reparametrized. Only recommended for the :math:`\Delta`-flow.

- ``SBE_DECOMPOSITION``
    | 0 for false; 1 for true.
    | If true, the SBE decomposition is used to parametrize the vertex and the flow equations. Only implemented in the MF!

- ``SELF_ENERGY_FLOW_CORRECTIONS``
    | 0 for false; 1 for true.
    | If true, corrections to the flow equations for the vertex from the self-energy, starting at :math:`\ell=3`, are included.

- ``STATIC_FEEDBACK``
    | If defined, use static :math:`K_1` inter-channel feedback. Only makes sense for pure :math:`K_1` calculations.

- ``SWITCH_SUM_N_INTEGRAL``
    | 0 for false; 1 for true.
    | If true, the sum over internal Keldysh indices is done before the frequency integration.
    | **Recommended setting**: 1.

- ``USE_ANDERSON_ACCELERATION``
    | 0 for false; 1 for true.
    | If true, Anderson acceleration is used to converge parquet iterations and self-energy iterations in mfRG faster.

- ``USE_MPI``
    | If defined, MPI is used for parallelization across multiple nodes.

- ``USE_SBEb_MFRG_EQS``
    | Determines which version of the SBE approximation shall be used.
    | 0 for SBEa, 1 for SBEb.
    | **Warning**: Only implemented in the MF!

- ``VECTORIZED_INTEGRATION``
    | 0 for false; 1 for true.
    | If true, integrals are performed with vector-valued integrands. For Keldysh, vectorization over Keldysh indices. For Matsubara at finite :math:`T`, vectorization over the Matsubara sum.

- ``ZERO_TEMP``
    | 0 for false; 1 for true.
    | If true, temperature :math:`T = 0` is assumed.


Global parameters
=================


- ``converged_tol``
    Tolerance for loop convergence in mfRG.

- ``COUNT``
    Used to set the number of frequency points in the MF. For details, see the definitions in the file ``frequency_parameters.hpp``.

- ``Delta_factor_K1``
    Scale factor for the frequency grid of the :math:`K_1` vertex class.

- ``Delta_factor_SE``
    Scale factor for the frequency grid of the self-energy.

- ``Delta_factor_K2_w``
    Scale factor for the frequency grid of the bosonic frequency of the :math:`K_2` and :math:`K_{2'}` vertex classes.

- ``Delta_factor_K2_v``
    Scale factor for the frequency grid of the fermionic frequency of the :math:`K_2` and :math:`K_{2'}` vertex classes.

- ``Delta_factor_K3_w``
    Scale factor for the frequency grid of the bosonic frequency of the :math:`K_3` vertex class.

- ``Delta_factor_K3_v``
    Scale factor for the frequency grid of the fermionic frequencies of the :math:`K_3` vertex class.

- ``EQUILIBRIUM``
    If ``true``, use equilibrium FDTs for propagators.

- ``glb_mu``
    Chemical potential -- w.l.o.g. **ALWAYS** set to 0.0 for the AM!

- ``integrator_tol``
    Integrator tolerance.

- ``inter_tol``
    Tolerance for closeness to grid points when interpolating.

- ``INTERPOLATION``
    | Interpolation method to me used.
    | ``linear``: linear interpolation on the frequency grid.
    | ``linear_on_aux``: linear interpolation on the grid for the auxiliary frequency :math:`\Omega`.
    | ``cubic``: Interpolation with cubic splines (**warning**: expensive!).


- ``Lambda_ini``
    Initial value of the regulator :math:`\Lambda` for an mfRG flow.

- ``Lambda_fin``
    Final value of the regulator :math:`\Lambda` for an mfRG flow.

- ``Lambda_scale``
    Scale of the log substitution, relevant in the hybridization flow.

- ``dLambda_initial``
    Initial step size for ODE solvers with adaptive step size control.

- ``nBOS``
    Number of bosonic frequency points for the :math:`K_1` vertex class.

- ``nFER``
    Number of fermionic frequency points for the self-energy.

- ``nBOS2``
    Number of bosonic frequency points for the :math:`K_2` and :math:`K_{2'}` vertex classes.

- ``nFER2``
    Number of fermionic frequency points for the :math:`K_2` and :math:`K_{2'}` vertex classes.

- ``nBOS3``
    Number of bosonic frequency points for the :math:`K_3` vertex class.

- ``nFER3``
    Number of fermionic frequency points for the :math:`K_3` vertex class.

- ``U_NRG``
    Vector with the values of :math:`U` in units of :math:`\Delta` that an mfRG flow should cover. Serve as checkpoints for the flow. Useful for benchmarking purposes if data from other methods at precise parameter points are available.

- ``VERBOSE``
    | If ``true``, detailed information about all computational steps is written into the log file.
    | Recommended setting for production runs: ``false``

- ``nmax_Selfenergy_iterations``
    | Maximal number of self-energy iterations to be done during an mfRG flow for :math:`\ell \geq 3`.
    | Default value: 10

- ``tol_selfenergy_correction_abs``
    | Absolute tolerance for self-energy iterations in mfRG.
    | Default value: :math:`10^{-9}`

- ``tol_selfenergy_correction_rel``
    | Relative tolerance for self-energy iterations in mfRG.
    | Default value: :math:`10^{-5}`