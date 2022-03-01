"""
A Python wrapper for the compiled GG functions.
"""

import ctypes
import ctypes.util
import os

import numpy as np

from . import docs_generator
from . import utility

# Attempt to load the compiled C code
__lib_found = False
__libgg_path = None
cgg = None

# First check the local folder
try:
    abs_path = os.path.dirname(os.path.abspath(__file__))
    cgg = np.ctypeslib.load_library("gg", abs_path)
    __libgg_path = os.path.join(abs_path, cgg._name)
    __lib_found = True
except OSError:
    try:
        cgg = np.ctypeslib.load_library("libgg", abs_path)
        __libgg_path = os.path.join(abs_path, cgg._name)
        __lib_found = True
    except OSError:
        cgg = None

__order_enum = {
    "spherical": {
        "cca": 300,
        "gaussian": 301,
    },
    "cartesian": {
        "cca": 400,
        "molden": 401,
    }
}


def _build_collocation_ctype(nout, orbital=False):
    """
    Builds the ctypes signatures for the libgg C lib
    """
    ret = [
        # L, npoints
        ctypes.c_int,
        ctypes.c_ulong,

        # XYZ, stride
        np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags=("C", "A")),
        ctypes.c_ulong,

        # Gaussian, nprim, coef, exp, center
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags=("C", "A")),  # coef
        np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags=("C", "A")),  # exp
        np.ctypeslib.ndpointer(dtype=np.double, shape=(3, ), ndim=1, flags=("C", "A")),  # center

        # Spherical
        ctypes.c_int,
    ]
    if orbital:
        ret.insert(1, ctypes.c_ulong)  # norbs
        ret.insert(1, np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags=("C", "A")))  # orbs

    # Pushback output
    for n in range(nout):
        ret.append(np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags=("W", "C", "A")))

    return tuple(ret)


# Bind the C object
if cgg is not None:

    # Helpers
    cgg.gg_ncomponents.argtypes = (ctypes.c_int, ctypes.c_int)
    cgg.gg_ncomponents.restype = ctypes.c_int

    # Transposes
    cgg.gg_naive_transpose.restype = None
    cgg.gg_naive_transpose.argtypes = (ctypes.c_ulong, ctypes.c_ulong, np.ctypeslib.ndpointer(),
                                       np.ctypeslib.ndpointer())

    cgg.gg_fast_transpose.restype = None
    cgg.gg_fast_transpose.argtypes = (ctypes.c_ulong, ctypes.c_ulong, np.ctypeslib.ndpointer(),
                                      np.ctypeslib.ndpointer())

    # Collocation generators
    cgg.gg_orbitals.restype = None
    cgg.gg_orbitals.argtypes = _build_collocation_ctype(1, orbital=True)

    cgg.gg_collocation.restype = None
    cgg.gg_collocation.argtypes = _build_collocation_ctype(1)

    cgg.gg_collocation.restype = None
    cgg.gg_collocation_deriv1.argtypes = _build_collocation_ctype(4)

    cgg.gg_collocation.restype = None
    cgg.gg_collocation_deriv2.argtypes = _build_collocation_ctype(10)

    cgg.gg_collocation.restype = None
    cgg.gg_collocation_deriv3.argtypes = _build_collocation_ctype(20)


def c_compiled():
    """
    Checks if the c code has been compiled or not.
    """
    return __lib_found


def _validate_c_import():
    """
    Throws an error if libgg is not found.
    """
    if c_compiled() is False:
        raise ImportError("Compiled libgg not found. Please compile gau2grid before calling these\n"
                          "  functions. Alternatively, use the NumPy-based collocation functions found at\n"
                          "  gau2grid.ref.collocation or gau2grid.ref.collocation_basis. It should\n"
                          "  be noted that these functions are dramatically slower (4-20x).\n")


def cgg_path():
    """
    Returns the path to the found libgg.so/dylib/dll object.
    """
    _validate_c_import()
    return __libgg_path


def get_cgg_shared_object():
    """
    Returns the compiled C shared object.
    """
    _validate_c_import()

    return cgg


def max_L():
    """
    Return the maximum compiled angular momentum.
    """

    return cgg.gg_max_L()


def ncomponents(L, spherical=True):
    """
    Computes the number of components for spherical and cartesian gaussians of a given L

    Parameters
    ----------
    L : int
        The angular momentum of the basis function
    spherical : bool, optional
        Spherical (True) or Cartesian (False) number of components

    Returns
    -------
    int
        The number of components in the gaussian.
    """

    return cgg.gg_ncomponents(L, spherical)


def _wrapper_checks(L, xyz, spherical, spherical_order, cartesian_order):
    if L > cgg.gg_max_L():
        raise ValueError("LibGG was only compiled to AM=%d, requested AM=%d." % (cgg.max_L(), L))

    # Check XYZ
    if xyz.shape[0] != 3:
        raise ValueError("XYZ array must be of shape (3, N), found %s" % str(xyz.shape))


# Validate the input
    try:
        if spherical:
            order_name = spherical_order
            order_enum = __order_enum["spherical"][spherical_order]
        else:
            order_name = cartesian_order
            order_enum = __order_enum["cartesian"][cartesian_order]
    except KeyError:
        raise KeyError("Order Spherical=%s:%s not understood." % (spherical, order_name))

    return order_enum


def collocation_basis(xyz, basis, grad=0, spherical=True, out=None, cartesian_order="cca", spherical_order="cca"):

    return utility.wrap_basis_collocation(collocation,
                                          xyz,
                                          basis,
                                          grad,
                                          spherical=spherical,
                                          out=out,
                                          cartesian_order=cartesian_order,
                                          spherical_order=spherical_order)


# Write common docs
collocation_basis.__doc__ = docs_generator.build_collocation_basis_docs(
    "This function uses a optimized C library as a backend.")


def orbital_basis(orbs, xyz, basis, spherical=True, out=None, cartesian_order="cca", spherical_order="cca"):

    return utility.wrap_basis_orbital(orbital,
                                      orbs,
                                      xyz,
                                      basis,
                                      spherical=spherical,
                                      out=out,
                                      cartesian_order=cartesian_order,
                                      spherical_order=spherical_order)


orbital_basis.__doc__ = docs_generator.build_orbital_basis_docs(
    "This function uses a optimized C library as a backend.")


def collocation(xyz,
                L,
                coeffs,
                exponents,
                center,
                grad=0,
                spherical=True,
                out=None,
                cartesian_order="cca",
                spherical_order="cca"):

    # Validates we loaded correctly
    _validate_c_import()

    order_enum = _wrapper_checks(L, xyz, spherical, spherical_order, cartesian_order)

    # Check gaussian
    coeffs = np.asarray(coeffs, dtype=np.double)
    exponents = np.asarray(exponents, dtype=np.double)
    center = np.asarray(center, dtype=np.double)
    if coeffs.shape[0] != exponents.shape[0]:
        raise ValueError("Coefficients (N=%d) and exponents (N=%d) must have the same shape." %
                         (coeffs.shape[0], exponents.shape[0]))

    # Find the output shape
    npoints = xyz.shape[1]
    if spherical:
        nvals = utility.nspherical(L)
    else:
        nvals = utility.ncartesian(L)

    # Build the outputs
    out = utility.validate_coll_output(grad, (nvals, npoints), out)

    # Select the correct function
    if grad == 0:
        cgg.gg_collocation(L, npoints, xyz.ravel(), 1, coeffs.shape[0], coeffs, exponents, center, order_enum,
                           out["PHI"])
    elif grad == 1:
        cgg.gg_collocation_deriv1(L, npoints, xyz.ravel(), 1, coeffs.shape[0], coeffs, exponents, center, order_enum,
                                  out["PHI"], out["PHI_X"], out["PHI_Y"], out["PHI_Z"])
    elif grad == 2:
        cgg.gg_collocation_deriv2(L, npoints, xyz.ravel(), 1, coeffs.shape[0], coeffs, exponents, center, order_enum,
                                  out["PHI"], out["PHI_X"], out["PHI_Y"], out["PHI_Z"], out["PHI_XX"], out["PHI_XY"],
                                  out["PHI_XZ"], out["PHI_YY"], out["PHI_YZ"], out["PHI_ZZ"])
    elif grad == 3:
        cgg.gg_collocation_deriv3(L, npoints, xyz.ravel(), 1, coeffs.shape[0], coeffs, exponents, center, order_enum,
                                  out["PHI"], out["PHI_X"], out["PHI_Y"], out["PHI_Z"], out["PHI_XX"], out["PHI_XY"],
                                  out["PHI_XZ"], out["PHI_YY"], out["PHI_YZ"], out["PHI_ZZ"], out["PHI_XXX"],
                                  out["PHI_XXY"], out["PHI_XXZ"], out["PHI_XYY"], out["PHI_XYZ"], out["PHI_XZZ"],
                                  out["PHI_YYY"], out["PHI_YYZ"], out["PHI_YZZ"], out["PHI_ZZZ"])
    else:
        raise ValueError("Only up to grad=3 is supported.")

    return out


collocation.__doc__ = docs_generator.build_collocation_docs("This function uses a optimized C library as a backend.")


def orbital(orbs,
            xyz,
            L,
            coeffs,
            exponents,
            center,
            spherical=True,
            out=None,
            cartesian_order="cca",
            spherical_order="cca"):

    # Validates we loaded correctly
    _validate_c_import()

    order_enum = _wrapper_checks(L, xyz, spherical, spherical_order, cartesian_order)

    # Check gaussian
    orbs = np.asarray(orbs, dtype=np.double)
    coeffs = np.asarray(coeffs, dtype=np.double)
    exponents = np.asarray(exponents, dtype=np.double)
    center = np.asarray(center, dtype=np.double)
    if coeffs.shape[0] != exponents.shape[0]:
        raise ValueError("Coefficients (N=%d) and exponents (N=%d) must have the same shape." %
                         (coeffs.shape[0], exponents.shape[0]))

    # Find the output shape
    npoints = xyz.shape[1]
    if spherical:
        nvals = utility.nspherical(L)
    else:
        nvals = utility.ncartesian(L)

    if nvals != orbs.shape[1]:
        raise ValueError("Orbital block, must be equal to the shell size.")

    # Build the outputs
    if out is not None:
        out = {"PHI": out}
    out = utility.validate_coll_output(0, (orbs.shape[0], npoints), out)["PHI"]

    # Select the correct function
    cgg.gg_orbitals(L, orbs, orbs.shape[0], npoints, xyz.ravel(), 1, coeffs.shape[0], coeffs, exponents, center,
                    order_enum, out)

    return out


orbital.__doc__ = docs_generator.build_orbital_docs("This function uses a optimized C library as a backend.")
