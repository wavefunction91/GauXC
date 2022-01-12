"""
Contains several testing helper function
"""

import glob
import os

import numpy as np


def compare_collocation_results(test, ref):
    if set(test) != set(ref):
        raise KeyError("Test and Ref results dicts do not match")

    for k in ref.keys():
        match = np.allclose(test[k], ref[k], atol=1.e-14, rtol=1.e-10)
        if not match:
            tnorm = np.linalg.norm(test[k])
            rnorm = np.linalg.norm(ref[k])
            raise ValueError("Test (norm=%10.9f) and Ref (norm=%10.9f) results do not match for %s" %
                             (tnorm, rnorm, k))


def find_pygau2grid():
    """
    Finds a compiled pygau2grid code and attempts to run it
    """
    base_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    # Can expand this later
    found = False
    search_paths = ["objdir"]
    for path in search_paths:
        folder = os.path.join(base_folder, path)
        find = glob.glob(os.path.join(folder, "pygau2grid") + "*.so")
        if len(find) == 1:
            found = os.path.dirname(find[0])
            break
        elif len(find) > 1:
            raise ImportError("Found multiple pygau2grid's. How is that possible?")

    return found
