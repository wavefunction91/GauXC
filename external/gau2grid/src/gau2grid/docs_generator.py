"""
Contains several docstrings as there are several duplicate functions
"""

__doc_header = r"""

    .. math::

        \phi_{m p} = Y_\ell^m \sum_i c_i e^{ -\alpha_{i} | \phi_{\rm center} - p | ^2}

    Where for a given angular momentum :math:`\ell`, components :math:`m` range from :math:`+\ell` to :math:`-\ell`
    for each grid point :math:`p`.

"""

__basis_str = """
    basis : list of dicts
        Each dict should contain the following keys (L, coeffs, exponents, center).

        L : int
            The angular momentum of the gaussian
        coeffs : array_like
            The coefficients of the gaussian
        exponents : array_like
            The exponents of the gaussian
        center : array_like
            The cartesian center of the gaussian"""

__doc_notes = ""
# __doc_notes = """
#     Notes
#     -----
#     For cartesian output the "row" order is used:
#     L_0 = .
#     L_1 = X, Y, Z
#     L_2 = XX, XY, XZ, YY, YZ, ZZ.
#     ...

#     For spherical harmonics a 0-based ordering is used:
#     L_0 = R_00
#     L_1 = R_10, R_11c, R_11s
#     L_2 = R_20, R_21c, R_21s, R_22c, R_22s
#     ...
# """


def build_collocation_docs(insert=""):

    doc_header = "    Computes the collocation matrix for a given gaussian basis of the form:"
    doc_header += __doc_header

    param_data = """

    Parameters
    ----------
    xyz : array_like
        The ``(3, N)`` cartesian points to compute the grid on
    L : int
        The angular momentum of the gaussian
    coeffs : array_like
        The coefficients of the gaussian
    exponents : array_like
        The exponents of the gaussian
    center : array_like
        The cartesian center of the gaussian
    grad : int, optional (default: 0)
        Can return cartesian gradient and Hessian per point if requested.
    spherical : bool, optional (default: True)
        Transform the resulting cartesian gaussian to spherical
    out : dict, optional
        A dictionary of output NumPy arrays to write the data to.
    %s
    Returns
    -------
    dict of array_like
        Returns a dictionary containing the requested arrays (``PHI``, ``PHI_X``, ``PHI_XX``, etc).
        Where each matrix is of shape ``(ngaussian_basis x npoints)``
    """

    ret = doc_header
    if insert == "":
        ret += "\n"
    else:
        ret += "    " + insert

    ret += param_data % __doc_notes
    return ret


def build_orbital_docs(insert=""):

    doc_header = "    Computes a array of a given orbital on a grid for a given gaussian basis of the form:"
    doc_header += __doc_header

    param_data = """

    Parameters
    ----------
    orbitals : array_like
        The ``(norb, nval)`` section of orbitals.
    xyz : array_like
        The ``(3, N)`` cartesian points to compute the grid on
    L : int
        The angular momentum of the gaussian
    coeffs : array_like
        The coefficients of the gaussian
    exponents : array_like
        The exponents of the gaussian
    center : array_like
        The cartesian center of the gaussian
    spherical : bool, optional (default: True)
        Transform the resulting cartesian gaussian to spherical
    out : dict, optional
        A dictionary of output NumPy arrays to write the data to.
    %s

    Returns
    -------
    array_like
        Returns a ``(norb, N)`` array of the orbitals on a grid.
    """

    ret = doc_header
    if insert == "":
        ret += "\n"
    else:
        ret += "    " + insert

    ret += param_data % __doc_notes
    return ret


def build_collocation_basis_docs(insert=""):

    doc_header = "    Computes the collocation matrix for a given gaussian basis of the form:"
    doc_header += __doc_header

    param_data = """

    xyz : array_like
        The ``(3, N)`` cartesian points to compute the grid on
    %s
    grad : int, default=0
        Can return cartesian gradient and Hessian per point if requested.
    spherical : bool, default=True
        Transform the resulting cartesian gaussian to spherical
    out : dict, optional
        A dictionary of output NumPy arrays to write the data to.
    %s

    Returns
    -------
    dict of array_like
        Returns a dictionary containing the requested arrays (``PHI``, ``PHI_X``, ``PHI_XX``, etc).
        Where each matrix is of shape (ngaussian_basis x npoints)
    """

    ret = doc_header
    if insert == "":
        ret += "\n"
    else:
        ret += "    " + insert

    ret += param_data % (__basis_str, __doc_notes)
    return ret


def build_orbital_basis_docs(insert=""):

    doc_header = "    Computes a array of a given orbital on a grid for a given gaussian basis of the form:"
    doc_header += "    " + __doc_header

    param_data = """

    orbital : array_line
        A ``(norb, nao)`` orbital array aligned to the orbitals basis
    xyz : array_like
        The (3, N) cartesian points to compute the grid on
    %s
    spherical : bool, default=True
        Transform the resulting cartesian gaussian to spherical
    out : dict, optional
        A dictionary of output NumPy arrays to write the data to.
    %s
    Returns
    -------
    array_like
        Returns a ``(norb, N)`` array of the orbitals on a grid.
    """

    ret = doc_header
    # if insert == "":
    #     ret += "\n"
    # else:
    #     ret += "    " + insert

    ret += param_data % (__basis_str, __doc_notes)
    return ret
