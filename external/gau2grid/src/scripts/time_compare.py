import numpy as np
import time

np.random.seed(0)

import gau2grid as gg

### Options

#npoints = int(5)
npoints = int(1.e5)

L = 2
nprim = 1

spherical = False
spherical = True

do_transpose = False
#do_transpose = True

### Test

xyz = np.random.rand(3, npoints)

grad_inds = ["PHI_X", "PHI_Y", "PHI_Z"]
hess_inds = ["PHI_XX", "PHI_XY", "PHI_XZ", "PHI_YY", "PHI_YZ", "PHI_ZZ"]


def compare(test, ref, grad):
    """
    Compares two results
    """
    print("%-6s %s" % ("PHI", np.allclose(test["PHI"], ref["PHI"])))
    if grad > 0:
        print('--')
        for key in grad_inds:
            print("%-6s %s" % (key, np.allclose(test[key], ref[key])))
    if grad > 1:
        print('--')
        for key in hess_inds:
            print("%-6s %s" % (key, np.allclose(test[key], ref[key])))


def transpose_dict(inp, out):
    return

#pygg.fast_transpose(inp[k], out[k])


def build_out(nvals, npoints, grad):
    """
    Builds output zeros to prevent cost effecting timings
    """
    inds = ["PHI"]
    if grad > 0:
        inds += grad_inds
    if grad > 1:
        inds += hess_inds

    return {k: np.zeros((nvals, npoints)) for k in inds}


ncart = int((L + 1) * (L + 2) / 2)
nsph = L * 2 + 1

nvals = ncart
if spherical:
    nvals = nsph

coefs = np.arange(nprim, dtype=np.double) + 1
exps = np.arange(nprim, dtype=np.double) + 2
#center = np.array([5, 5, 5])
center = np.array([0, 0, 0], dtype=np.double)

### Points

# Call pyGG
gg_out = build_out(nvals, npoints, 0)
tran_out = build_out(npoints, nvals, 0)
t = time.time()
if do_transpose:
    transpose_dict(gg_out, tran_out)

gg.collocation(xyz, L, coefs, exps, center, grad=0, spherical=spherical, out=gg_out)
#gg_out["PHI"] = gg_out["PHI"].copy().reshape(npoints, nvals).T
ctime = (time.time() - t)

# Call NP GG
t = time.time()
np_out = gg.ref.collocation(xyz, L, coefs, exps, center, grad=0, spherical=spherical, cartesian_order="row")
pytime = (time.time() - t)

# print(c_out.shape)
# print(np_out["PHI"].shape)
print("PHI")
compare(gg_out, np_out, 0)
#print(np_out["PHI"])
#print(gg_out["PHI"])

print("C time  %12.6f" % ctime)
print("Py time %12.6f" % pytime)
print("Ratio   %12.6f" % (pytime / ctime))

### Derivatives

print("\nDerivative")
# Call pyGG
gg_out = build_out(nvals, npoints, 1)
tran_out = build_out(npoints, nvals, 1)
t = time.time()
if do_transpose:
    transpose_dict(gg_out, tran_out)
gg.collocation(xyz, L, coefs, exps, center, grad=1, spherical=spherical, out=gg_out)
ctime = (time.time() - t)

# Call NP GG
t = time.time()
np_out = gg.ref.collocation(xyz, L, coefs, exps, center, grad=1, spherical=spherical, cartesian_order="row")
pytime = (time.time() - t)

compare(gg_out, np_out, 1)

print("C time  %12.6f" % ctime)
print("Py time %12.6f" % pytime)
print("Ratio   %12.6f" % (pytime / ctime))

### Hessian

print("\nHessian")
gg_out = build_out(nvals, npoints, 2)
tran_out = build_out(npoints, nvals, 2)
t = time.time()
gg.collocation(xyz, L, coefs, exps, center, grad=2, spherical=spherical, out=gg_out)
if do_transpose:
    transpose_dict(gg_out, tran_out)
ctime = (time.time() - t)

# Call NP GG
t = time.time()
np_out = gg.ref.collocation(xyz, L, coefs, exps, center, grad=2, spherical=spherical, cartesian_order="row")
pytime = (time.time() - t)
#print(np_out["PHI_X"])
#print(np_out["PHI_Y"])
#print(np_out["PHI_Z"])

compare(gg_out, np_out, 2)

print("C time  %12.6f" % ctime)
print("Py time %12.6f" % pytime)
print("Ratio   %12.6f" % (pytime / ctime))
