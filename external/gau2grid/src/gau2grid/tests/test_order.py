"""
Tests the cartesian ordering code.
"""

import pytest

import gau2grid as gg

_benchmark_data = {
    "molden": {
        0: [""],
        1: ["X", "Y", "Z"],
        2: ["XX", "YY", "ZZ", "XY", "XZ", "YZ"],
        3: ["XXX", "YYY", "ZZZ", "XYY", "XXY", "XXZ", "XZZ", "YZZ", "YYZ", "XYZ"],
        4: [
            "XXXX", "YYYY", "ZZZZ", "XXXY", "XXXZ", "XYYY", "YYYZ", "XZZZ", "YZZZ", "XXYY", "XXZZ", "YYZZ", "XXYZ",
            "XYYZ", "XYZZ"
        ],
    },
    "row": {
        0: [""],
        1: ["X", "Y", "Z"],
        2: ["XX", "XY", "XZ", "YY", "YZ", "ZZ"],
        3: ["XXX", "XXY", "XXZ", "XYY", "XYZ", "XZZ", "YYY", "YYZ", "YZZ", "ZZZ"],
        4: [
            "XXXX", "XXXY", "XXXZ", "XXYY", "XXYZ", "XXZZ", "XYYY", "XYYZ", "XYZZ", "XZZZ", "YYYY", "YYYZ", "YYZZ",
            "YZZZ", "ZZZZ"
        ],
    },
    "libint": {
        0: [""],
        1: ["X", "Y", "Z"],
        2: ["XX", "YY", "ZZ", "XY", "YZ", "XZ"],
        3: ["XXX", "YYY", "ZZZ", "XXY", "XXZ", "XYY", "YYZ", "XZZ", "YZZ", "XYZ"],
        4: [
            "XXXX", "YYYY", "ZZZZ", "XXXY", "XXXZ", "XYYY", "YYYZ", "XZZZ", "YZZZ", "XXYY", "XXZZ", "YYZZ", "XXYZ",
            "XYYZ", "XYZZ"
        ],
    }
}


@pytest.mark.parametrize("order", ["molden", "row"])
@pytest.mark.parametrize("L", [0, 1, 2, 3, 4])
def test_cartesian_order(order, L):

    data = _benchmark_data[order][L]

    order_list = []
    for idx, l, m, n in gg.order.cartesian_order_factory(L, order=order):
        order = "X" * l + "Y" * m + "Z" * n
        assert order == data[idx]
        order_list.append(order)

    # Check all values are unique
    assert len(order_list) == len(set(order_list))

    # Check all lengths are correct
    assert all(len(x) == L for x in order_list)
