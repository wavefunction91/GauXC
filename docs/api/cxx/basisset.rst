Basis set and shell types
=========================

These APIs define the core data structures for Gaussian-type orbital (GTO) basis sets.
The :cpp:struct:`GauXC::BasisSet` is an ordered collection of shells that make up a full basis for a molecule or atom set.
A :cpp:class:`GauXC::Shell` represents a contracted Gaussian shell with its exponents, coefficients, origin, and angular momentum.

.. doxygenstruct:: GauXC::BasisSet
   :members:

.. doxygenclass:: GauXC::Shell
   :members: