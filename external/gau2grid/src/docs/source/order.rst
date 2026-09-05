.. _gpo_order:

Gaussian Component Orders
=========================

The order of the individual components can vary between use cases. gau2grid can
produce any resulting order that a user requires. The ``C`` version of the code
must be compiled to a given order. The currently supported orders are as
follows.


Cartesian Order
---------------

gau2grid currently supports both the ``cca`` and ``molden`` orders. The number
of components per angular momentum can be computed as:

.. math::

    N_{\rm cartesian} = (\ell + 1) (\ell + 2) / 2

Row Order
+++++++++

The ``cca`` order iterates over the upper triangular hyper diagonal and has the
following pattern:

 - ``S`` (:math:`\ell = 0`): `1`
 - ``P`` (:math:`\ell = 1`): ``X``, ``Y``, ``Z``
 - ``D`` (:math:`\ell = 2`): ``XX``, ``XY``, ``XZ``, ``YY``, ``YZ``, ``ZZ``
 - ``F`` (:math:`\ell = 3`): ``XXX``, ``XXY``, ``XXZ``, ``XYY``, ``XYZ``, ``XZZ``, ``YYY``, ``YYZ``, ``YZZ``, ``ZZZ``

Molden Order
++++++++++++

The ``molden`` order is primarily found in a Molden format and only has a
determined values for :math:`0 \leq \ell < 4`.

 - ``S`` (:math:`\ell = 0`): `1`
 - ``P`` (:math:`\ell = 1`): ``X``, ``Y``, ``Z``
 - ``D`` (:math:`\ell = 2`): ``XX``, ``YY``, ``ZZ``, ``XY``, ``XZ``, ``YZ``
 - ``F`` (:math:`\ell = 3`): ``XXX``, ``YYY``, ``ZZZ``, ``XYY``, ``XXY``, ``XXZ``, ``XZZ``, ``YZZ``, ``YYZ``, ``XYZ``



Spherical Order
---------------

gau2grid currently supports both the ``CCA`` and ``gaussian`` orders. The numnber
of components per angular momentum can be computed as:

.. math::

     N_{\rm spherical} = 2\ell + 1


CCA Order
+++++++++

An industry standard order known as the Common Component Architecture:

 - ``S`` (:math:`\ell = 0`): :math:`Y_0^0`
 - ``P`` (:math:`\ell = 1`): :math:`Y_1^{-1}`, :math:`Y_1^{0}`, :math:`Y_1^{+1}`,
 - ``D`` (:math:`\ell = 2`): :math:`Y_2^{-2}`, :math:`Y_2^{-1}`, :math:`Y_2^{0}`, :math:`Y_2^{+1}`, :math:`Y_2^{+2}`

Gaussian Order
++++++++++++++

The ``gaussian`` order as used by the Gaussian program:

 - ``S`` (:math:`\ell = 0`): :math:`Y_0^0`
 - ``P`` (:math:`\ell = 1`): :math:`Y_1^{0}`, :math:`Y_1^{+1}`, :math:`Y_1^{-1}`,
 - ``D`` (:math:`\ell = 2`): :math:`Y_2^{0}`, :math:`Y_2^{+1}`, :math:`Y_2^{-1}`, :math:`Y_2^{+2}`, :math:`Y_2^{-2}`
