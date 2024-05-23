=====================
Correlation Functions
=====================

This is the documentation of the classes that have been implemented to represent the correlation functions and related objects
that are computed in all numerical calculations. Those are the different vertex classes and the class used to
represent the self-energy, which are used together to represent the ''state'' of a calculation. In addition, there are
two helper classes, which don't hold any additional non-trivial information, but which are used for convenience.
First, there is a class to represent propagators, which combines the trivial non-interacting part determined by the
physical system and the non-trivial self-energy. Second, there is another class that combines one pair of propagators
in such a way that it is used to contract bubble-type diagrams.

Vertices
========

In total, the code contains the four classes ``irreducible``, ``rvert``, ``fullvert``, and ``GeneralVertex`` to store different types of four-point vertices.

.. doxygenclass:: irreducible
   :members:

.. doxygenclass:: rvert
   :members:

.. doxygenclass:: fullvert
   :members:

.. doxygenclass:: GeneralVertex
   :members:


Self-Energy
===========

.. doxygenclass:: SelfEnergy
   :members:

State
=====

Instances of the ``State`` class are the high-level objects that are mainly used by the high-level algorithms
discussed in :ref:`algorithms`. The ``State`` class combines a ``GeneralVertex`` and a ``SelfEnergy``,
which together contain all non-trivial information that one might wish to compute.

.. doxygenclass:: State
   :members:


Propagator
==========

The ``Propagator`` class is special in the sense that it stores almost no data itself. Instead, it references
instances of the ``SelfEnergy`` class, which hold the non-trivial data from the computations. The ``Propagator``
class is then used to combine the analytical form of the bare propagator :math:`G_0` with the self-energy via the
Dyson equation, :math:`G = 1 / [(G_0)^{-1} - \Sigma]`.

.. doxygenclass:: Propagator
   :members:

Bubble
======

Lastly, the ``Bubble`` class combines two propagators to yield a bubble in one of the three two-particle channels
:math:`a`, :math:`p` and :math:`t`.
For evaluating differentiated bubbles in the context of mfRG computations, one of the
propagators can be chosen to be the single-scale propagator :math:`S` or the fully differentiated one :math:`\dot{G}`.
In that case, the bubble already takes care of the product rule, giving (symbolically)
:math:`\dot{\Pi}^S = G\cdot S + S \cdot G` or :math:`\dot{\Pi} = G\cdot\dot{G} + \dot{G}\cdot G`.
Otherwise, it just yields :math:`\Pi = G\cdot G`.

.. doxygenclass:: Bubble
   :members:
