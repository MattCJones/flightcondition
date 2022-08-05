#!/usr/bin/env python
"""
Test normal shock computations.

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

from common import assert_field
from flightcondition import NormalShock, dimless

# Atmospheric ground truth data
M = [1.0, 3.0] * dimless
M2_truth_arr = [1.0, 0.47519096] * dimless
p2_by_p1_truth_arr = [1.0, 10.3333333] * dimless
rho2_by_rho1_truth_arr = [1.0, 3.85714285] * dimless
T2_by_T1_truth_arr = [1.0, 2.67901234] * dimless
p02_by_p01_truth_arr = [1.0, 0.32834388] * dimless


def test_M2():
    """Test Mach number. """
    assert_field(NormalShock.M2(M), M2_truth_arr)


def test_p2_by_p1():
    """Test static pressure. """
    assert_field(NormalShock.p2_by_p1(M), p2_by_p1_truth_arr)


def test_rho2_by_rho1():
    """Test density. """
    assert_field(NormalShock.rho2_by_rho1(M), rho2_by_rho1_truth_arr)


def test_T2_by_T1():
    """Test temperature. """
    assert_field(NormalShock.T2_by_T1(M), T2_by_T1_truth_arr)


def test_p02_by_p01():
    """Test stagnation pressure. """
    assert_field(NormalShock.p02_by_p01(M), p02_by_p01_truth_arr)
