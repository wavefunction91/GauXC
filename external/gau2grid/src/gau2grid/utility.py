"""
Provides utility functions for the gau2grid program
"""

import numpy as np


def get_deriv_indices(grad):
    """
    Returns the indices of the derivatives involved in the grid derivatives

    Examples
    --------
    >>> get_deriv_indices(1)
    ["x", "y", "z"]
    """
    if grad == 0:
        return []
    elif grad == 1:
        return ["x", "y", "z"]
    elif grad == 2:
        return ["x", "y", "z", "xx", "xy", "xz", "yy", "yz", "zz"]
    elif grad == 3:
        return ["x", "y", "z", "xx", "xy", "xz", "yy", "yz", "zz",
                "xxx", "xxy", "xxz", "xyy", "xyz", "xzz", "yyy", "yyz", "yzz", "zzz" ]
    else:
        raise ValueError("Only grid derivatives up to grad=3 are supported.")


def get_output_keys(grad):
    """
    Returns the output keys required for a given derivative

    Examples
    --------
    >>> get_output_keys(1)
    ["PHI", "PHI_X", "PHI_Y", "PHI_Z"]
    """

    phi = ["PHI"]

    if grad == 0:
        return phi

    deriv_keys = ["PHI_" + x.upper() for x in get_deriv_indices(grad)]
    return phi + deriv_keys


def validate_coll_output(grad, shape, out):
    """
    Validates the a collocation output, constructs a new
    output array if necessary
    """
    keys_needed = get_output_keys(grad)
    if out is None:
        out = {k: np.zeros(shape) for k in keys_needed}
    else:
        if not isinstance(out, dict):
            raise TypeError("Output parameter must be a dictionary.")
        missing = set(keys_needed) - set(out)
        if len(missing):
            raise KeyError("Missing output keys '%s'" % str(missing))

        for key in keys_needed:
            out[key] = np.asarray(out[key])
            if out[key].shape != shape:
                raise ValueError(
                    "Shape of each output array must be (ntotal, npoints). Shape of key '%s' is incorrect." % key)
    return out


def nspherical(L):
    """
    Computes the number of spherical functions for a given angular momentum.

    Parameters
    ----------
    L : int
        The input angular momentum

    Returns
    -------
    nspherical : int
        The number of spherical functions
    """

    return L * 2 + 1


def ncartesian(L):
    """
    Computes the number of cartesian functions for a given angular momentum.

    Parameters
    ----------
    L : int
        The input angular momentum

    Returns
    -------
    ncartesian : int
        The number of cartesian functions
    """

    return int((L + 1) * (L + 2) / 2)


def _parse_basis(basis, spherical):
    # Check the basis
    parsed_basis = []
    for num, func in enumerate(basis):
        # TODO more checks

        # Either list or dict form
        if isinstance(func, (list, tuple)):
            if len(func) != 4:
                raise ValueError("Basis should have 4 components (L, coeffs, exponents, center).")
            parsed_basis.append(func)

        elif isinstance(func, dict):
            missing = {"am", "coef", "exp", "center"} - set(func)
            if len(missing):
                raise KeyError("Missing function keys '%s'" % str(missing))

            tmp = [func["am"], func["coef"], func["exp"], func["center"]]
            parsed_basis.append(tmp)
        else:
            raise TypeError("Basis type not recognized!")

    # The total number of output parameters
    if spherical:
        nfunc = [nspherical(func[0]) for func in parsed_basis]
        ntotal = sum(nfunc)
    else:
        nfunc = [ncartesian(func[0]) for func in parsed_basis]
        ntotal = sum(nfunc)

    return parsed_basis, nfunc, ntotal


def wrap_basis_collocation(coll_function, xyz, basis, grad, spherical, out, cartesian_order, spherical_order):
    """
    Wraps collocation computers to apply to entire basis sets.

    Expects the basis to take the form of:
        [L, coeffs, exponents, center]
    """

    # A few checkers
    if grad > 3:
        raise IndexError("Can only compute up to 3rd derivatives of the grid (grad=3).")

    parsed_basis, nfunc, ntotal = _parse_basis(basis, spherical)
    npoints = xyz.shape[1]

    # Handle output
    out = validate_coll_output(grad, (ntotal, npoints), out)

    # Loop over functions in the basis set
    start = 0
    for n, func in enumerate(parsed_basis):
        # Build slice
        nvals = nfunc[n]
        sl = slice(start, start + nvals)
        start += nvals

        # Build temporary output views
        tmp_out = {k: v[sl] for k, v in out.items()}

        coll_function(
            xyz,
            *func,
            grad=grad,
            spherical=spherical,
            out=tmp_out,
            cartesian_order=cartesian_order,
            spherical_order=spherical_order)

    return out


def wrap_basis_orbital(orbital_function, orbs, xyz, basis, spherical, out, cartesian_order, spherical_order):
    """
    Wraps orbital computers to apply to entire basis sets.

    Expects the basis to take the form of:
        [L, coeffs, exponents, center]
    """

    parsed_basis, nfunc, ntotal = _parse_basis(basis, spherical)
    npoints = xyz.shape[1]
    norbs = orbs.shape[0]

    # Handle output
    if out is not None:
        out = {"PHI": out}
    out = validate_coll_output(0, (norbs, npoints), out)["PHI"]

    # Loop over functions in the basis set
    start = 0
    for n, func in enumerate(parsed_basis):
        # Build slice
        nvals = nfunc[n]
        sl = slice(start, start + nvals)
        start += nvals

        # Build temporary output views
        tmp_orbs = np.array(orbs[:, sl])

        orbital_function(
            tmp_orbs,
            xyz,
            *func,
            spherical=spherical,
            out=out,
            cartesian_order=cartesian_order,
            spherical_order=spherical_order)

    return out
