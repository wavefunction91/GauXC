Collocation Example
===================

Single Collocation
------------------

A collocation grid between a single basis and a Cartesian grid can be computed
with the :func:`~gau2grid.collocation` function. For example, we will use a grid
starting at the origin along the ``z`` axis:

.. code-block:: python

    >>> import gau2grid
    >>> import numpy as np
    >>> xyz = np.zeros((3, 5))
    >>> xyz[2] = np.arange(5)

We can then create a gaussian with only a single coefficient and exponent of 1
centered on the origin:

.. code-block:: python

    >>> L = 0
    >>> coef = [1]
    >>> exp = [1]
    >>> center = [0, 0, 0]

The collocation grid can then be computed as:

.. code-block:: python

    >>> ret = gau2grid.collocation(xyz, L, coef, exp, center)
    >>> ret["PHI"]
    [[  1.00000e+00   3.67879e-01   1.83156e-02   1.23409e-04   1.12535e-07]]

The ``p`` gaussian can be also be computed. Note that since our grid points are
along the ``z`` axis, the ``x`` and ``y`` components are orthogonal and thus
zero. 

.. code-block:: python

    >>> L = 1
    >>> ret = gau2grid.collocation(xyz, L, coef, exp, center, spherical=False, grad=1)
    >>> ret["PHI"]
    [[  0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]  # P_x
     [  0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]  # P_y
     [  0.00000e+00   3.67879e-01   3.66312e-02   3.70229e-04   4.50140e-07]] # P_z

As the previous execution used ``grad=1``, the ``X``, ``Y``, and ``Z``
cartesian gradients are also available and can be accessed as:

.. code-block:: python

    >>> ret["PHI_Z"]
    [[  0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]
     [  0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]
     [  1.00000e+00  -3.67879e-01  -1.28209e-01  -2.09797e-03  -3.48859e-06]]

Basis Collocation
-----------------

Often it is beneficial to compute the collocation matrix between several basis
functions and a set of grid points at once the
:func:`~gau2grid.collocation_basis` helper function provides this
functionality. To begin, a set of basis sets can be constructed with the
following form:

.. code-block:: python

    >>> basis = [{
        'center': [0., 0., 0.],
        'exp': [38, 6, 1],
        'coef': [0.4, 0.6, 0.7],
        'am': 0
    }, {
        'center': [0., 0., 0.],
        'exp': [0.3],
        'coef': [0.3],
        'am': 1
    }]

Execution of this basis results in a collocation matrix where basis results are
vertically stacked on top of each other:

.. code-block:: python

    >>> ret = gau2grid.collocation_basis(xyz, basis, spherical=False)
    >>> ret["PHI"]
    [[  1.70000e+00   2.59003e-01   1.28209e-02   8.63869e-05   7.87746e-08]  # S
     [  0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]  # P_x
     [  0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]  # P_y
     [  0.00000e+00   2.22245e-01   1.80717e-01   6.04850e-02   9.87570e-03]] # P_z
