#!/usr/bin/env python
"""
Configure top level, common  processes for testing.

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

# flake8: noqa W503 W504

import pytest

reltol = 0.0001


def myapprox(test_arg):
    """Appoximate value of testarg to specified relative tolerance.

    Args:
        test_arg: argument to approximate

    Returns:
        approximation of argument
    """
    return pytest.approx(test_arg, rel=reltol)


def assert_field(field_arr, field_truth_arr):
    """Test that output field matches truth data.

    Args:
        field_arr: field array of type pint.Quantity
    """
    assert (
           field_arr.to_base_units().magnitude
           ==
           pytest.approx(field_truth_arr.to_base_units().magnitude, rel=reltol)
           )
