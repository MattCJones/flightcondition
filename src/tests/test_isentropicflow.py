#!/usr/bin/env python
"""
Test isentropic flow relations.

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

import numpy as np

from common import assert_field
from flightcondition import IsentropicFlow, dimless
from flightcondition.constants import PhysicalConstants as Phys

# Atmospheric ground truth data
M = [0.5, 1.0, 3.0] * dimless
p0_by_p_truth_arr = [1/0.84301917, 1/0.52828178, 1/0.02722368] * dimless
T0_by_T_truth_arr = [1/0.95238095, 1/0.83333333, 1/0.35714285] * dimless
rho0_by_rho_truth_arr = [1/0.88517013, 1/0.63393814, 1/0.07622631] * dimless

gamma_air = Phys.gamma_air
R_air = Phys.R_air
a0_by_a_truth_arr = np.sqrt(T0_by_T_truth_arr)

A_star_by_A_truth_arr = [1/1.33984375, 1.0, 1/4.23456790] * dimless


def test_p0_by_p():
    """Test pressure ratio. """
    assert_field(IsentropicFlow.p0_by_p(M), p0_by_p_truth_arr)
    assert_field(1/IsentropicFlow.p_by_p0(M), p0_by_p_truth_arr)


def test_T0_by_T():
    """Test temperature ratio. """
    assert_field(IsentropicFlow.T0_by_T(M), T0_by_T_truth_arr)
    assert_field(1/IsentropicFlow.T_by_T0(M), T0_by_T_truth_arr)


def test_a0_by_a():
    """Test sound speed ratio. """
    assert_field(IsentropicFlow.a0_by_a(M), a0_by_a_truth_arr)
    assert_field(1/IsentropicFlow.a_by_a0(M), a0_by_a_truth_arr)


def test_rho0_by_rho():
    """Test density ratio. """
    assert_field(IsentropicFlow.rho0_by_rho(M), rho0_by_rho_truth_arr)
    assert_field(1/IsentropicFlow.rho_by_rho0(M), rho0_by_rho_truth_arr)


def test_A_star_by_A():
    """Test area ratio. """
    assert_field(IsentropicFlow.A_star_by_A(M), A_star_by_A_truth_arr)
    assert_field(1/IsentropicFlow.A_by_A_star(M), A_star_by_A_truth_arr)


def test_M_from_p0_by_p():
    """Test Mach number from pressure ratio. """
    assert_field(IsentropicFlow.M_from_p0_by_p(p0_by_p_truth_arr), M)
    assert_field(IsentropicFlow.M_from_p_by_p0(1/p0_by_p_truth_arr), M)


def test_M_from_T0_by_T():
    """Test Mach number from temperature ratio. """
    assert_field(IsentropicFlow.M_from_T0_by_T(T0_by_T_truth_arr), M)
    assert_field(IsentropicFlow.M_from_T_by_T0(1/T0_by_T_truth_arr), M)


def test_M_from_a0_by_a():
    """Test Mach number from sound speed. """
    assert_field(IsentropicFlow.M_from_a0_by_a(a0_by_a_truth_arr), M)
    assert_field(IsentropicFlow.M_from_a_by_a0(1/a0_by_a_truth_arr), M)


def test_M_from_rho0_by_rho():
    """Test Mach number from density. """
    assert_field(IsentropicFlow.M_from_rho0_by_rho(rho0_by_rho_truth_arr), M)
    assert_field(IsentropicFlow.M_from_rho_by_rho0(1/rho0_by_rho_truth_arr), M)
