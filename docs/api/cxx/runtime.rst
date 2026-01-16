Runtime environment
===================

These APIs manage execution contexts and work distribution for XC integration.
The :cpp:class:`GauXC::RuntimeEnvironment` encapsulates MPI communication and provides
process-local information, while :cpp:class:`GauXC::DeviceRuntimeEnvironment` extends
it with GPU memory management for device-accelerated computations.

.. doxygenclass:: GauXC::RuntimeEnvironment
   :members:

.. doxygenclass:: GauXC::DeviceRuntimeEnvironment
   :members: