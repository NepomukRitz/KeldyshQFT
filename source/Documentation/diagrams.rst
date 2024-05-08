Computing diagrams
==================

Bubbles
+++++++

.. doxygenfunction:: bubble_function(vertexType_result& dgamma, const vertexType_left& vertex1, const vertexType_right& vertex2, const Bubble_Object& Pi, const char channel, const fRG_config& config, const std::array<bool,3> tobecomputed = {true,true,true})

.. doxygenfunction:: bubble_function(vertexType_result& dgamma, const vertexType_left& vertex1, const vertexType_right& vertex2, const Propagator<Q>& G, const Propagator<Q>& S, const char channel, const bool diff, const fRG_config& config, const std::array<bool,3> tobecomputed = {true,true,true})

.. doxygenclass:: BubbleFunctionCalculator
    :members:

.. doxygenclass:: Integrand
    :members:

Loops
+++++

.. doxygenfunction:: loop(SelfEnergy<Q>& self, const vertType& fullvertex, const Propagator<Q>& prop)

.. doxygenclass:: LoopCalculator
    :members:

.. doxygenclass:: IntegrandSE
    :members:
