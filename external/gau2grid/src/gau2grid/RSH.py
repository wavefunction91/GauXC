"""
Cartesian to regular solid harmonics conversion code.
"""

import decimal
import os
import pickle
import platform

import numpy as np

from . import order
from . import utility

_MAX_AM = 17
_DECIMAL_PREC = 60
_saved_rsh_coefs = {}
_saved_factorials = {}


def _factorial(n):
    decimal.getcontext().prec = _DECIMAL_PREC
    if n in _saved_factorials:
        return _saved_factorials[n]

    if n == 0:
        return decimal.Decimal("1.0")
    else:
        return n * _factorial(n - 1)


class RSH_Memoize(object):
    """
    Simple memoize class for RSH_coefs which is quite expensive
    """

    def __init__(self, func):
        self.func = func
        self.mem = {}

    def __call__(self, AM, **kwargs):

        # Bypass Memoize for testing
        if kwargs.get("force_call", False):
            return self.func(AM)

        if AM not in self.mem:
            self.mem[AM] = self.func(AM)

        return self.mem[AM]


@RSH_Memoize
def _cart_to_RSH_coeffs_gen(l):
    """
    Generates a coefficients [ coef, x power, y power, z power ] for each component of
    a regular solid harmonic (in terms of raw Cartesians) with angular momentum l.

    See eq. 23 of ACS, F. C. Pickard, H. F. Schaefer and B. R. Brooks, JCP, 140, 184101 (2014)

    Returns coeffs with order 0, +1, -1, +2, -2, ...
    """

    # Arbitrary precision math with 50 decimal places
    decimal.getcontext().prec = _DECIMAL_PREC

    terms = []
    for m in range(l + 1):
        thisterm = {}
        p1 = ((_factorial(l - m)) / (_factorial(l + m))).sqrt() * ((_factorial(m)) / (2**l))
        if m:
            p1 *= decimal.Decimal("2.0").sqrt()

        # Loop over cartesian components
        for lz in range(l + 1):
            for ly in range(l - lz + 1):

                lx = l - ly - lz
                xyz = lx, ly, lz
                j = int((lx + ly - m) / 2)
                if (lx + ly - m) % 2 == 1 or j < 0:
                    continue

                # P2
                p2 = decimal.Decimal(0.0)
                for i in range(int((l - m) / 2) + 1):
                    if i >= j:
                        p2 += (-1)**i * _factorial(2 * l - 2 * i) / (_factorial(l - i) * _factorial(i - j) *
                                                                     _factorial(l - m - 2 * i))

                # P3
                p3 = decimal.Decimal(0.0)
                for k in range(j + 1):
                    if (j >= k) and (lx >= 2 * k) and (m + 2 * k >= lx):
                        p3 += (-1)**k / (_factorial(j - k) * _factorial(k) * _factorial(lx - 2 * k) *
                                         _factorial(m - lx + 2 * k))

                p = p1 * p2 * p3

                # Add in part if not already present
                if xyz not in thisterm:
                    thisterm[xyz] = [decimal.Decimal(0.0), decimal.Decimal(0.0)]

                # Add the two components
                if (m - lx) % 2:
                    # imaginary
                    sign = decimal.Decimal(-1.0)**decimal.Decimal((m - lx - 1) / 2.0)
                    thisterm[xyz][1] += sign * p
                else:
                    # real
                    sign = decimal.Decimal(-1.0)**decimal.Decimal((m - lx) / 2.0)
                    thisterm[xyz][0] += sign * p

        tmp_R = []
        tmp_I = []
        for k, v in thisterm.items():
            if abs(v[0]) > 0:
                tmp_R.append((k, v[0]))
            if abs(v[1]) > 0:
                tmp_I.append((k, v[1]))

        if m == 0:
            # name_R = "R_%d%d" % (l, m)
            terms.append(tmp_R)
        else:
            # name_R = "R_%d%dc" % (l, m)
            # name_I = "R_%d%ds" % (l, m)
            terms.append(tmp_R)
            terms.append(tmp_I)
            # terms[name_R] = tmp_R
            # terms[name_I] = tmp_I

        # for k, v in terms.items():
        #     print(k, v)

    return terms


def cart_to_RSH_coeffs(L, order="gaussian", force_call=False):
    """
    Allows coefficients either to be generated or pulled from disk

    Allowed orders:
        "gaussian":
            R_0, R^+_1, R^-_1, ..., R^+_l, R^-_l
        "CCA":
            R^-_(l), R^-_(l-1), ..., R_0, ..., R^+_(l-1), R^+_l
    """

    # Gen the coefficients (may be memoized)
    data = _cart_to_RSH_coeffs_gen(L, force_call=force_call)

    if order.lower() == "gaussian":
        return data
    elif order.lower() == "cca":
        ret = []

        # Add negative
        for l in range(L):
            ret.append(data[2 + l * 2])

        # Reverse so we get (-L, 0) not (0, L)
        ret.reverse()

        # Add in zero
        ret.append(data[0])

        # Add positive
        for l in range(L):
            ret.append(data[1 + l * 2])

        return ret

    else:
        raise KeyError("Order '%s' not understood" % order)


def cart_to_spherical_transform(data, L, cartesian_order, spherical_order):
    """
    Transforms a cartesian x points matrix into a spherical x points matrix.
    """

    cartesian_order = {x[1:]: x[0] for x in order.cartesian_order_factory(L, cartesian_order)}
    RSH_coefs = cart_to_RSH_coeffs(L, order=spherical_order)

    nspherical = len(RSH_coefs)
    ret = np.zeros((nspherical, data.shape[1]))

    idx = 0
    for spherical in RSH_coefs:
        for cart_index, scale in spherical:
            ret[idx] += float(scale) * data[cartesian_order[cart_index]]
        idx += 1

    return ret


def transformation_c_generator(cg, L, cartesian_order, spherical_order, function_name="", prefix=None, align=32):
    """
    Builds a conversion from cartesian to spherical coordinates in C
    """

    if function_name == "":
        if prefix:
            function_name = "gg_%s_cart_to_spherical_L%d" % (prefix, L)
        else:
            function_name = "gg_cart_to_spherical_L%d" % L

    cartesian_order = {x[1:]: x[0] for x in order.cartesian_order_factory(L, cartesian_order)}
    RSH_coefs = cart_to_RSH_coeffs(L, order=spherical_order)

    signature = "void %s(const unsigned long size, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT spherical, const unsigned long nspherical)" % function_name

    # Start function
    cg.start_c_block(signature)
    cg.write("ASSUME_ALIGNED(cart, %d)" % align)

    cg.write("// R_%d0 Transform" % L)
    _c_spherical_trans(cg, 0, RSH_coefs, cartesian_order)
    cg.blankline()

    for l in range(L):
        cg.write("// R_%d%dc Transform" % (L, l + 1))
        sidx = 2 * l + 1
        _c_spherical_trans(cg, sidx, RSH_coefs, cartesian_order)

        sidx = 2 * l + 2
        cg.write("// R_%d%ds Transform" % (L, l + 1))
        _c_spherical_trans(cg, sidx, RSH_coefs, cartesian_order)
        cg.blankline()

    # End function
    cg.close_c_block()
    return signature


def _c_spherical_trans(cg, sidx, RSH_coefs, cartesian_order):
    # cg.write("#pragma clang loop vectorize(assume_safety)")
    cg.start_c_block("for (unsigned long i = 0; i < size; i++)")

    # Figure out where we are summing to
    if sidx == 0:
        lhs = "spherical[i]"
    elif sidx == 1:
        lhs = "spherical[nspherical + i]"
    else:
        lhs = "spherical[%d * nspherical + i]" % sidx

    op = " ="
    for cart_index, scale in RSH_coefs[sidx]:

        # Figure out car idx
        idx = cartesian_order[cart_index]
        if idx == 0:
            rhs = "cart[i]"
        elif idx == 1:
            rhs = "cart[ncart + i]"
        else:
            rhs = "cart[%d * ncart + i]" % idx

        # Scales
        if scale != 1.0:
            cg.write("%s %s % .16f * %s" % (lhs, op, scale, rhs))
        else:
            cg.write("%s %s %s" % (lhs, op, rhs))
        op = "+="
    cg.blankline()

    cg.close_c_block()


def transformation_c_generator_sum(cg, L, cartesian_order, spherical_order, function_name="", prefix=None, align=32):
    """
    Builds a conversion from cartesian to spherical coordinates in C
    """

    if function_name == "":
        if prefix:
            function_name = "gg_%s_cart_to_spherical_sum_L%d" % (prefix, L)
        else:
            function_name = "gg_cart_to_spherical_sum_L%d" % L

    cartesian_order = {x[1:]: x[0] for x in order.cartesian_order_factory(L, cartesian_order)}
    RSH_coefs = cart_to_RSH_coeffs(L, order=spherical_order)

    signature = "void %s(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT output, const unsigned long nspherical)" % function_name

    # Start function
    cg.start_c_block(signature)
    cg.write("ASSUME_ALIGNED(cart, %d)" % align)

    cg.write("// temps")
    cg.write("double tmp")

    cg.write("// R_%d0 Transform" % L)
    _c_spherical_trans_vector_sum(cg, 0, RSH_coefs, cartesian_order)
    cg.blankline()

    for l in range(L):
        cg.write("// R_%d%dc Transform" % (L, l + 1))
        sidx = 2 * l + 1
        _c_spherical_trans_vector_sum(cg, sidx, RSH_coefs, cartesian_order)

        sidx = 2 * l + 2
        cg.write("// R_%d%ds Transform" % (L, l + 1))
        _c_spherical_trans_vector_sum(cg, sidx, RSH_coefs, cartesian_order)
        cg.blankline()

    # End function
    cg.close_c_block()
    return signature


def _c_spherical_trans_vector_sum(cg, sidx, RSH_coefs, cartesian_order):
    # cg.write("#pragma clang loop vectorize(assume_safety)")
    cg.start_c_block("for (unsigned long i = 0; i < size; i++)")

    lhs = "tmp"

    op = " ="
    for cart_index, scale in RSH_coefs[sidx]:

        # Figure out car idx
        idx = cartesian_order[cart_index]
        if idx == 0:
            rhs = "cart[i]"
        elif idx == 1:
            rhs = "cart[ncart + i]"
        else:
            rhs = "cart[%d * ncart + i]" % idx

        # Scales
        if scale != 1.0:
            cg.write("%s %s % .16f * %s" % (lhs, op, scale, rhs))
        else:
            cg.write("%s %s %s" % (lhs, op, rhs))
        op = "+="

    cg.write("output[i] += tmp * vector[%s]" % sidx)
    cg.blankline()

    cg.close_c_block()
