Load balancing
==============

The :cpp:class:`GauXC::LoadBalancer` distributes quadrature tasks across processes and provides access to local work units.
To create load balancers, the :cpp:class:`GauXC::LoadBalancerFactory` is provided, which can generate different load balancing strategies based on user-defined settings.
For specifying the execution context, the :cpp:enum:`GauXC::ExecutionSpace` enum is available.

.. doxygenclass:: GauXC::LoadBalancer
   :members:

.. doxygenclass:: GauXC::LoadBalancerFactory
   :members:

.. doxygenenum:: GauXC::ExecutionSpace