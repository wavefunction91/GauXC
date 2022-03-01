.. gau2grid documentation master file, created by
   sphinx-quickstart on Sat Sep  1 17:41:02 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

========
gau2grid
========

*gau2grid is a python-generated C library for vectorized computation of grid to gaussian collocation matrices*

The core of gau2grid is generating the collocation matrices between a real
space grid and a gaussian basis set expanded to a given angular momenta.
Where a simple gaussian can be represented with the cartesian form as:

.. math::

    \phi({\bf r}) = x^l y^m z^n e^{-\alpha r^2}

where for a given angular momenta :math:`\ell`, a gaussian basis has all
possible combinations of :math:`l, m, n` that satisfy :math:`l + m + n =
\ell`. These gaussians can also take a `spherical harmonic <https://en.wikipedia.org/wiki/Spherical_harmonics>`_ form of:

.. math::

    \phi({\bf r}) = Y_\ell^m (\hat{\bf r}) e^{-\alpha r^2}

where :math:`m` ranges from :math:`+\ell` to :math:`-\ell`. The spherical
form offers a more compact representation at higher angular momenta, but is
more difficult to work with when examining cartesian derivates.

In quantum chemistry, an individual basis is often represented as a sum of
several gaussian with different exponents and coefficients together:

.. math::

    \phi({\bf r}) = Y_\ell^m (\hat{\bf r}) \sum_i c_i e^{-\alpha_i r^2}

Collocation matrices between a single basis set and multiple grid points can
then be represented as follows:

.. math::

    \phi_{m p} = Y_\ell^m (\widehat{{\bf r}_p -{\bf r}_{\rm center}}) \sum_i c_i e^{-\alpha_i ({\bf r}_{\rm center} - {\bf r}_p) ^2}

where the basis is evaluated at every point :math:`p` for every
component of the basis i.e. basis function :math:`m`. The
:math:`\phi_{m p}` matrices are the primary focus on the gau2grid
library.


Index
-----

**Getting Started**

* :doc:`py_install`
* :doc:`c_install`
* :doc:`order`

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Getting Started

   py_install
   c_install
   order

**Python Interface**

* :doc:`py_example`
* :doc:`py_api`

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Python API

   py_api
   py_example

**C Interface**

* :doc:`c_example`
* :doc:`c_api`

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: C API

   c_api
   c_example
