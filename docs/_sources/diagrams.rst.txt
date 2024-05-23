==================
Computing diagrams
==================

All diagrammatic computations for self-energies and vertices in the setting discussed here ultimately require the
computation of diagrams of one of two types, which we call ''bubbles'' and ''loops''. The former uses an instance of
the ``Bubble`` class to connect two four-point vertices with a propagator pair, yielding another four-point object.
The latter contracts two external fermionic legs of a four-point vertex with a propagator, yielding a two-point object.

Bubbles
=======

.. doxygenfunction:: bubble_function(vertexType_result& dgamma, const vertexType_left& vertex1, const vertexType_right& vertex2, const Bubble_Object& Pi, const char channel, const fRG_config& config, const std::array<bool,3> tobecomputed = {true,true,true})

.. doxygenfunction:: bubble_function(vertexType_result& dgamma, const vertexType_left& vertex1, const vertexType_right& vertex2, const Propagator<Q>& G, const Propagator<Q>& S, const char channel, const bool diff, const fRG_config& config, const std::array<bool,3> tobecomputed = {true,true,true})

.. doxygenclass:: BubbleFunctionCalculator
    :members:

.. doxygenclass:: Integrand
    :members:

Loops
=====

.. doxygenfunction:: loop(SelfEnergy<Q>& self, const vertType& fullvertex, const Propagator<Q>& prop)

.. doxygenclass:: LoopCalculator
    :members:

.. doxygenclass:: IntegrandSE
    :members:
