Perturbation theory
===================

Hartree-Fock
++++++++++++

.. doxygenclass:: Hartree_Solver
    :members:

Second-order perturbation theory (PT2)
++++++++++++++++++++++++++++++++++++++

.. doxygenfunction:: sopt_state(State<Q>& Psi, const bool diff = false)

.. doxygenfunction:: sopt_state_impl(State<Q>& Psi, const Bubble_Object& Pi, const State<Q>& bareState)

.. doxygenfunction:: vertexInSOPT(Vertex<Q, false>& PsiVertex, const State<Q>& bareState, const Bubble_Object& Pi)

.. doxygenfunction:: selfEnergyInSOPT(SelfEnergy<Q>& PsiSelfEnergy, const State<Q>& bareState, const Bubble_Object& Pi)
