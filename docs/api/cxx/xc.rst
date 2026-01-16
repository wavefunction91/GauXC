Exchange-correlation integrator
===============================

For performing exchange-correlation integration, GauXC provides the :cpp:class:`GauXC::XCIntegrator` interface.
The :cpp:class:`GauXC::XCIntegratorFactory` creates integrators based on user-defined settings encapsulated in :cpp:struct:`GauXC::IntegratorSettingsXC`, :cpp:struct:`GauXC::IntegratorSettingsKS`, and :cpp:struct:`GauXC::IntegratorSettingsEXC_GRAD`.

.. doxygenclass:: GauXC::XCIntegrator
   :members:

.. doxygenclass:: GauXC::XCIntegratorFactory
   :members:

.. doxygenstruct:: GauXC::IntegratorSettingsXC
   :members:

.. doxygenstruct:: GauXC::IntegratorSettingsKS
   :members:

.. doxygenstruct:: GauXC::IntegratorSettingsEXC_GRAD
   :members:

.. doxygenstruct:: GauXC::IntegratorSettingsEXX
   :members:

.. doxygenstruct:: GauXC::IntegratorSettingsSNLinK
   :members: