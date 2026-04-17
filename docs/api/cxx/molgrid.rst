Molecular grid types
====================

These APIs provide the quadrature grid infrastructure for numerical integration.
The :cpp:class:`GauXC::Grid` class encapsulates an atom-centered spherical quadrature with batching support.
:cpp:class:`GauXC::MolGrid` aggregates per-element grids based on the atomic number to form a complete molecular integration grid.
The :cpp:class:`GauXC::MolGridFactory` creates molecular grids from molecules and basis sets.

For applying standard grid sizes and pruning schemes, the enums :cpp:enum:`GauXC::AtomicGridSizeDefault`, :cpp:enum:`GauXC::RadialQuad`, and :cpp:enum:`GauXC::PruningScheme` are provided.
Finally, the :cpp:class:`GauXC::MolecularWeights` class provides support for computing molecular integration weights, with configuration via :cpp:struct:`GauXC::MolecularWeightsSettings` and construction through :cpp:class:`GauXC::MolecularWeightsFactory`.

.. doxygenclass:: GauXC::MolGrid
   :members:

.. doxygenclass:: GauXC::Grid
   :members:

.. doxygenclass:: GauXC::MolGridFactory
   :members:

.. doxygenenum:: GauXC::AtomicGridSizeDefault

.. doxygenenum:: GauXC::RadialQuad

.. doxygenenum:: GauXC::PruningScheme

.. doxygenclass:: GauXC::MolecularWeights
   :members:

.. doxygenclass:: GauXC::MolecularWeightsFactory
   :members:

.. doxygenstruct:: GauXC::MolecularWeightsSettings
   :members: