API Reference
=============

Helper Functions
++++++++++++++++

A collection of function ment to provide information and the gau2grid library.

.. c:function:: int gg_max_L();

    Returns the maximum compiled angular momentum

.. c:function:: int gg_ncomponents(const int L, const int spherical)

    Returns the number of components for a given angular momentum.

    :param L: The angular momentum of the basis function.
    :param spherical: Boolean that returns spherical (1) or cartesian (0) basis representations.

The following enums are also specified:

 - ``GG_SPHERICAL_CCA`` - CCA spherical output.
 - ``GG_SPHERICAL_GAUSSIAN`` - Gaussian spherical output.
 - ``GG_CARTESIAN_CCA`` - CCA cartesian output.
 - ``GG_CARTESIAN_MOLDEN`` - Molden cartesian output.

Transpose Functions
+++++++++++++++++++

Transposes matrices if input or output order is incorrect.

.. c:function:: void gg_naive_transpose(unsigned long n, unsigned long m, const double* PRAGMA_RESTRICT input, double* PRAGMA_RESTRICT output)

    Transposes a matrix using a simple for loop.

    :param n: The number of rows in the input matrix.
    :param m: The number of rows in the output matrix.
    :param input: The ``(n x m)`` input matrix.
    :param output: The ``(m x n)`` output matrix.


.. c:function:: void gg_fast_transpose(unsigned long n, unsigned long m, const double* PRAGMA_RESTRICT input, double* PRAGMA_RESTRICT output)

    Transposes a matrix using a small on-cache temporary array. Is usually faster than :c:func:`~gg_naive_transpose`.

    :param n: The number of rows in the input matrix.
    :param m: The number of rows in the output matrix.
    :param input: The ``(n x m)`` input matrix.
    :param output: The ``(m x n)`` output matrix.

Orbital Functions
+++++++++++++++++

Computes orbitals on a grid.


.. c:function:: void gg_orbitals(int L, const double* PRAGMA_RESTRICT C, const unsigned long norbitals, const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT orbital_out)

    Computes orbital a section on a grid. This function performs the following
    contraction inplace.

    .. math::

        C_{im} \phi_{m p} \rightarrow ret_{i p}

    This is often more efficient than generating :math:`\phi_{m p}` and then
    contracting with the orbitals C as there is greater cache locality.

    :param L: The angular momentum of the basis function.
    :param C: A ``(norbitals, ncomponents)`` matrix of orbital coefficients.
    :param norbitals: The number of orbs to compute.
    :param npoints: The number of grid points to compute.
    :param xyz: A ``(npoints, 3)`` or (npoints, n) array of the xyz coordinates.
    :param xyz_stride: The stride of the xyz input array. 1 for ``xx..., yy..., zz...`` style input, 3 for ``xyz, xyz, xyz, ...`` style input.
    :param nprim: The number of primitives (exponents and coefficients) in the basis set
    :param coeffs: A ``(nprim, )`` array of coefficients (:math:`c`).
    :param exponents: A ``(nprim, )`` array of exponents (:math:`\alpha`).
    :param center: A ``(3, )`` array of x, y, z coordinate of the basis center.
    :param order: Enum that specifies the output order.
    :param orbital_out: ``(norbitals, npoints)`` array of orbitals on the grid.

Collocation Functions
+++++++++++++++++++++

Creates collocation matrices between a gaussian function and a set of grid points.


.. c:function:: void gg_collocation(int L, const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out)

    Computes the collocation array:

    .. math::

        \phi_{m p} = Y_\ell^m \sum_i c_i e^{-\alpha_i |\phi_{\rm center} - p| ^2}

    :param L: The angular momentum of the basis function.
    :param npoints: The number of grid points to compute.
    :param xyz: A ``(npoints, 3)`` or (npoints, n) array of the xyz coordinates.
    :param xyz_stride: The stride of the xyz input array. 1 for ``xx..., yy..., zz...`` style input, 3 for ``xyz, xyz, xyz, ...`` style input.
    :param nprim: The number of primitives (exponents and coefficients) in the basis set
    :param coeffs: A ``(nprim, )`` array of coefficients (:math:`c`).
    :param exponents: A ``(nprim, )`` array of exponents (:math:`\alpha`).
    :param center: A ``(3, )`` array of x, y, z coordinate of the basis center.
    :param order: Enum that specifies the output order.
    :param phi_out: ``(ncomponents, npoints)`` collocation array.

.. c:function:: void gg_collocation_deriv1(int L, const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out)

    Computes the collocation array and the corresponding first cartesian derivatives:

    .. math::

        \phi_{m p} = Y_\ell^m \sum_i c_i e^{-\alpha_i |\phi_{\rm center} - p| ^2}

    :param L: The angular momentum of the basis function.
    :param npoints: The number of grid points to compute.
    :param xyz: A ``(npoints, 3)`` or (npoints, n) array of the xyz coordinates.
    :param xyz_stride: The stride of the xyz input array. 1 for ``xx..., yy..., zz...`` style input, 3 for ``xyz, xyz, xyz, ...`` style input.
    :param nprim: The number of primitives (exponents and coefficients) in the basis set
    :param coeffs: A ``(nprim, )`` array of coefficients (:math:`c`).
    :param exponents: A ``(nprim, )`` array of exponents (:math:`\alpha`).
    :param center: A ``(3, )`` array of x, y, z coordinate of the basis center.
    :param order: Enum that specifies the output order.
    :param phi_out: ``(ncomponents, npoints)`` collocation array.
    :param phi_x_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``x``.
    :param phi_y_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``y``.
    :param phi_z_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``z``.


.. c:function:: void gg_collocation_deriv2(int L, const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out, double* PRAGMA_RESTRICT phi_xx_out, double* PRAGMA_RESTRICT phi_xy_out, double* PRAGMA_RESTRICT phi_xz_out, double* PRAGMA_RESTRICT phi_yy_out, double* PRAGMA_RESTRICT phi_yz_out, double* PRAGMA_RESTRICT phi_zz_out)

    Computes the collocation array and the corresponding first and second cartesian derivatives:

    .. math::

        \phi_{m p} = Y_\ell^m \sum_i c_i e^{-\alpha_i |\phi_{\rm center} - p| ^2}

    :param L: The angular momentum of the basis function.
    :param npoints: The number of grid points to compute.
    :param xyz: A ``(npoints, 3)`` or (npoints, n) array of the xyz coordinates.
    :param xyz_stride: The stride of the xyz input array. 1 for ``xx..., yy..., zz...`` style input, 3 for ``xyz, xyz, xyz, ...`` style input.
    :param nprim: The number of primitives (exponents and coefficients) in the basis set
    :param coeffs: A ``(nprim, )`` array of coefficients (:math:`c`).
    :param exponents: A ``(nprim, )`` array of exponents (:math:`\alpha`).
    :param center: A ``(3, )`` array of x, y, z coordinate of the basis center.
    :param order: Enum that specifies the output order.
    :param phi_out: ``(ncomponents, npoints)`` collocation array.
    :param phi_x_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``x``.
    :param phi_y_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``y``.
    :param phi_z_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``z``.
    :param phi_xx_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``xx``.
    :param phi_xy_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``xy``.
    :param phi_xz_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``xz``.
    :param phi_yy_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``yy``.
    :param phi_yz_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``yz``.
    :param phi_zz_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``zz``.

.. c:function:: void gg_collocation_deriv3(int L, const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out, double* PRAGMA_RESTRICT phi_xx_out, double* PRAGMA_RESTRICT phi_xy_out, double* PRAGMA_RESTRICT phi_xz_out, double* PRAGMA_RESTRICT phi_yy_out, double* PRAGMA_RESTRICT phi_yz_out, double* PRAGMA_RESTRICT phi_zz_out, double* PRAGMA_RESTRICT phi_xxx_out, double* PRAGMA_RESTRICT phi_xxy_out, double* PRAGMA_RESTRICT phi_xxz_out, double* PRAGMA_RESTRICT phi_xyy_out, double* PRAGMA_RESTRICT phi_xyz_out, double* PRAGMA_RESTRICT phi_xzz_out, double* PRAGMA_RESTRICT phi_yyy_out, double* PRAGMA_RESTRICT phi_yyz_out, double* PRAGMA_RESTRICT phi_yzz_out, double* PRAGMA_RESTRICT phi_zzz_out)

    Computes the collocation array and the corresponding first, second, and third cartesian derivatives:

    .. math::

        \phi_{m p} = Y_\ell^m \sum_i c_i e^{-\alpha_i |\phi_{\rm center} - p| ^2}

    :param L: The angular momentum of the basis function.
    :param npoints: The number of grid points to compute.
    :param xyz: A ``(npoints, 3)`` or (npoints, n) array of the xyz coordinates.
    :param xyz_stride: The stride of the xyz input array. 1 for ``xx..., yy..., zz...`` style input, 3 for ``xyz, xyz, xyz, ...`` style input.
    :param nprim: The number of primitives (exponents and coefficients) in the basis set
    :param coeffs: A ``(nprim, )`` array of coefficients (:math:`c`).
    :param exponents: A ``(nprim, )`` array of exponents (:math:`\alpha`).
    :param center: A ``(3, )`` array of x, y, z coordinate of the basis center.
    :param order: Enum that specifies the output order.
    :param phi_out: ``(ncomponents, npoints)`` collocation array.
    :param phi_x_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``x``.
    :param phi_y_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``y``.
    :param phi_z_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``z``.
    :param phi_xx_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``xx``.
    :param phi_xy_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``xy``.
    :param phi_xz_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``xz``.
    :param phi_yy_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``yy``.
    :param phi_yz_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``yz``.
    :param phi_zz_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``zz``.
    :param phi_xxx_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``xxx``.
    :param phi_xxy_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``xxy``.
    :param phi_xxz_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``xxz``.
    :param phi_xyy_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``xyy``.
    :param phi_xyz_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``xyz``.
    :param phi_xzz_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``xzz``.
    :param phi_yyy_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``yyy``.
    :param phi_yyz_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``yyz``.
    :param phi_yzz_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``yzz``.
    :param phi_zzz_out: ``(ncomponents, npoints)`` collocation derivative with respect to ``zzz``.
