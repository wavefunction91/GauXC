"""
Additional functionality not directly related to gau2grid.
"""

import os

# Testing
def test():
    """
    Runs a smoke test suite through pytest.
    """

    try:
        import pytest
    except ImportError:
        raise RuntimeError('Testing module `pytest` is not installed. Run `conda install pytest`')

    abs_test_dir = os.path.sep.join([os.path.abspath(os.path.dirname(__file__)), "tests"])
    retcode = pytest.main(['-rws', '-v', '--capture=sys', abs_test_dir])
    return retcode
