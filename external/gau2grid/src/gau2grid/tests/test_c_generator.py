"""
Compare the generated C code against the NumPy reference code.
"""

import shutil
import tempfile

import numpy as np
import pytest

import gau2grid as gg

np.set_printoptions(linewidth=120, suppress=True)

# Import locals

# Simply test that the code runs for now
c_gen_tests = []
for AM in range(4):
    for grad in range(3):
        c_gen_tests.append((AM, grad))


@pytest.mark.parametrize("AM,grad", c_gen_tests)
def test_c_collocation_codgen(AM, grad):
    cg = gg.codegen.CodeGen(cgen=True)
    gg.c_gen.shell_c_generator(cg, AM, grad=grad)


@pytest.mark.parametrize("AM", list(range(4)))
def test_c_spherical_trans_codgen(AM):
    cg = gg.codegen.CodeGen(cgen=True)
    gg.RSH.transformation_c_generator(cg, AM, "row", "gaussian")


def test_library_gen():

    temp_dir = tempfile.mkdtemp()
    gg.c_gen.generate_c_gau2grid(4, path=temp_dir)
    shutil.rmtree(temp_dir)


def test_pybind11_gen():
    cg = gg.codegen.CodeGen(cgen=True)
    gg.c_util_generator.pybind11_func(cg, "something", 0, "somthingelse", 6)
    gg.c_util_generator.pybind11_func(cg, "something", 1, "somthingelse", 6)
    gg.c_util_generator.pybind11_func(cg, "something", 2, "somthingelse", 6)
    gg.c_util_generator.pybind11_transpose(cg, "t1", "t2")
