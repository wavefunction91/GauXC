"""
Compare the generated C code against the NumPy reference code.
"""

import os
import sys
import time

import numpy as np
import pytest

import gau2grid as gg

np.set_printoptions(linewidth=120, suppress=True)

# Import locals
from . import ref_basis
from . import test_helper as th

# Tweakers
npoints = int(1.e3)
npoints2 = int(npoints / 2)

# Global points
np.random.seed(0)
xyzw = np.random.rand(3, npoints)

# LR points
xyzw[:, npoints2:] += 5 * np.random.rand(3, npoints2)

# Make sure the C-side has been compiled
if "GAU2GRID_FORCE_C_TEST" in os.environ:
    skip_c_test = False
else:
    skip_c_test = gg.c_compiled() is False
check_compile = pytest.mark.skipif(skip_c_test, reason="Could not find the C compiled SO for gau2grid")

test_basis_keys = list(ref_basis.test_basis.keys())
# test_basis_keys = ["single-1s", "single-1p", "single-1d"]

test_orders = [
    ("cartesian", "cca"),
    ("cartesian", "molden"),
    ("spherical", "cca"),
    ("spherical", "gaussian")
] # yapf: disable

@check_compile
@pytest.mark.parametrize("basis_name", test_basis_keys)
@pytest.mark.parametrize("spherical, order_name", test_orders)
def test_generator_collocation(basis_name, spherical, order_name):
#
    kwargs = {"grad": 2, "spherical": "spherical" == spherical}

    if kwargs["spherical"]:
        kwargs["spherical_order"] = order_name
    else:
        kwargs["cartesian_order"] = order_name


    basis = ref_basis.test_basis[basis_name]

    max_L = max(x["am"] for x in basis)
    if (order_name == "molden") and (max_L > 4):
        pytest.skip("Molden only goes to L=4.")

    t = time.time()
    gen_results = gg.collocation_basis(xyzw, basis, **kwargs)
    gg_time = time.time() - t

    t = time.time()
    ref_results = gg.ref.collocation_basis(xyzw, basis, **kwargs)
    ref_time = time.time() - t

    # Print time with py.test -s flags
    print("")
    print("%s-%s time REF: %8.4f GG: %8.4f" % (basis_name, spherical, ref_time, gg_time))

    # print(ref_results["PHI"])
    # print(gen_results["PHI"])
    th.compare_collocation_results(gen_results, ref_results)

@check_compile
@pytest.mark.parametrize("xyz_shape", [3, 4, 5])
def test_generator_collocation_transposed(xyz_shape):

    cgg = gg.get_cgg_shared_object()

    # Collocation data
    npoints = 2000
    L = 2
    nelem = 2 * L + 1
    order_enum = 300
    coeffs = np.array([2., 1])
    exponents = np.array([1., 2])
    center = np.array([0., 0, 0])

    # Generate random points
    data = np.random.rand(xyz_shape, npoints)
    xyz = data[:3].copy()
    xyz_t = data.transpose().copy()

    out = np.zeros((nelem, npoints))
    out_t = np.zeros((nelem, npoints))

    cgg.gg_collocation(L, npoints, xyz.ravel(), 1, coeffs.shape[0], coeffs, exponents, center, order_enum,
                       out)
    cgg.gg_collocation(L, npoints, xyz_t.ravel(), xyz_shape, coeffs.shape[0], coeffs, exponents, center, order_enum,
                       out_t)
    th.compare_collocation_results({"PHI": out}, {"PHI": out_t})


@check_compile
@pytest.mark.parametrize("basis_name", test_basis_keys)
@pytest.mark.parametrize("spherical, order_name", test_orders)
def test_generator_orbital(basis_name, spherical, order_name):

    kwargs = {"grad": 0, "spherical": "spherical" == spherical}

    if kwargs["spherical"]:
        kwargs["spherical_order"] = order_name
    else:
        kwargs["cartesian_order"] = order_name

    basis = ref_basis.test_basis[basis_name]

    max_L = max(x["am"] for x in basis)
    if (order_name == "molden") and (max_L > 4):
        pytest.skip("Molden only goes to L=4.")

    t = time.time()
    phi = gg.collocation_basis(xyzw, basis, **kwargs)["PHI"]
    orbs = np.random.rand(5, phi.shape[0])
    benchmark = np.dot(orbs, phi)
    ref_time = time.time() - t

    t = time.time()
    del kwargs["grad"]
    ref_results = gg.orbital_basis(orbs, xyzw, basis, **kwargs)
    gg_time = time.time() - t

    # Print time with py.test -s flags
    print("")
    print("%s-%s time REF: %8.4f GG: %8.4f" % (basis_name, spherical, ref_time, gg_time))
    # print(benchmark)
    # print(ref_results)

    th.compare_collocation_results({"ORBITALS": benchmark}, {"ORBITALS": ref_results})

@check_compile
@pytest.mark.parametrize("xyz_shape", [3, 4, 5])
def test_generator_orbital_transposed(xyz_shape):

    cgg = gg.get_cgg_shared_object()

    # Collocation data
    npoints = 2000
    L = 2
    nelem = 2 * L + 1
    order_enum = 300
    coeffs = np.array([2., 1])
    exponents = np.array([1., 2])
    center = np.array([0., 0, 0])

    C = np.random.rand(2, nelem)

    # Generate random points
    data = np.random.rand(xyz_shape, npoints)
    xyz = data[:3].copy()
    xyz_t = data.transpose().copy()

    out = np.zeros((nelem, npoints))
    out_t = np.zeros((nelem, npoints))

    cgg.gg_orbitals(L, C, C.shape[0], npoints, xyz.ravel(), 1, coeffs.shape[0], coeffs, exponents, center, order_enum,
                       out)
    cgg.gg_orbitals(L, C, C.shape[0], npoints, xyz_t.ravel(), xyz_shape, coeffs.shape[0], coeffs, exponents, center, order_enum,
                       out_t)
    th.compare_collocation_results({"ORB": out}, {"ORB": out_t})


@check_compile
@pytest.mark.parametrize("grad", [0, 1, 2, 3])
def test_generator_derivs(grad):
    kwargs = {"spherical_order": "cca", "cartesian_order": "cca", "spherical": False, "grad": grad}

    basis = ref_basis.test_basis["cc-pVDZ"]

    gen_results = gg.collocation_basis(xyzw, basis, **kwargs)
    ref_results = gg.ref.collocation_basis(xyzw, basis, **kwargs)

    th.compare_collocation_results(gen_results, ref_results)


@check_compile
@pytest.mark.parametrize("grad", [0, 1, 2, 3])
def test_generator_derivs_spherical(grad):
    kwargs = {"spherical_order": "cca", "cartesian_order": "cca", "spherical": True, "grad": grad}

    basis = ref_basis.test_basis["cc-pVDZ"]

    gen_results = gg.collocation_basis(xyzw, basis, **kwargs)
    ref_results = gg.ref.collocation_basis(xyzw, basis, **kwargs)

    th.compare_collocation_results(gen_results, ref_results)


@check_compile
def test_libgg_path():
    assert "gg" in gg.cgg_path()


@check_compile
@pytest.mark.parametrize("am,spherical,result", [
    (0, True, 1),
    (0, False, 1),
    (1, True, 3),
    (1, False, 3),
    (2, True, 5),
    (2, False, 6),
    (3, True, 7),
    (3, False, 10),
])
def test_ncomponents(am, spherical, result):
    assert gg.ncomponents(am, spherical) == result


@check_compile
@pytest.mark.parametrize("spherical", [True, False])
@pytest.mark.parametrize("am", [0, 1, 2, 3, 4])
def test_generator_orbitals_am(spherical, am):
    kwargs = {"spherical_order": "cca", "cartesian_order": "cca", "spherical": spherical, "grad": 0}

    # Build a single orbital
    coeffs = [0.44135347600549724, 0.6934968471367846, 0.6641842253258472, 0.0001]
    exponents = [38.36, 5.77, 1.24, 1.e-2]
    center = [0., 0., 0.]
    L = am

    ret = gg.collocation(xyzw, L, coeffs, exponents, center, **kwargs)["PHI"]
    orb = np.random.rand(3, ret.shape[0])
    bench = np.dot(orb, ret)

    del kwargs["grad"]
    ret = gg.orbital(orb, xyzw, L, coeffs, exponents, center, **kwargs)

    # Compare the results
    th.compare_collocation_results({"ORBITALS": bench}, {"ORBITALS": ret})
