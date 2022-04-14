#!/usr/bin/env python
"""
Test flight condition functionality.

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

# flake8: noqa E203

import pytest

from numpy import array

from flightcondition import Atmosphere, FlightCondition, unit, dimless
from common import assert_field, myapprox

h_geom_arr = [0, 30e3] * unit('ft')


def test_TAS():
    """Test true airspeed calculations. """

    fc = FlightCondition(h_geom_arr, TAS=300*unit('knots'))

    TAS_truth = array([300, 300]) * unit('knots')
    assert_field(fc.TAS, TAS_truth)

    CAS_truth = array([300, 187.7518]) * unit('knots')
    assert_field(fc.CAS, CAS_truth)

    EAS_truth = array([300, 183.6448]) * unit('knots')
    assert_field(fc.EAS, EAS_truth)

    mach_truth = array([0.4535, 0.5090]) * dimless
    assert_field(fc.mach, mach_truth)


def test_CAS():
    """Test calibrated airspeed calculations. """

    fc = FlightCondition(h_geom_arr, CAS=300*unit('knots'))

    TAS_truth = array([300, 465.6309]) * unit('knots')
    assert_field(fc.TAS, TAS_truth)

    CAS_truth = array([300, 300]) * unit('knots')
    assert_field(fc.CAS, CAS_truth)

    EAS_truth = array([300, 285.0357]) * unit('knots')
    assert_field(fc.EAS, EAS_truth)

    mach_truth = array([0.4535, 0.7900]) * dimless
    assert_field(fc.mach, mach_truth)


def test_EAS():
    """Test equivalent airspeed calculations. """

    fc = FlightCondition(h_geom_arr, EAS=300*unit('knots'))

    TAS_truth = array([300, 490.0764]) * unit('knots')
    assert_field(fc.TAS, TAS_truth)

    CAS_truth = array([300, 317.3602]) * unit('knots')
    assert_field(fc.CAS, CAS_truth)

    EAS_truth = array([300, 300]) * unit('knots')
    assert_field(fc.EAS, EAS_truth)

    mach_truth = array([0.4535, 0.8314]) * dimless
    assert_field(fc.mach, mach_truth)


def test_mach():
    """Test Mach number calculations. """

    fc = FlightCondition(h_geom_arr, mach=0.88*dimless)

    TAS_truth = array([582.1012, 518.7004]) * unit('knots')
    assert_field(fc.TAS, TAS_truth)

    CAS_truth = array([582.1012, 337.977]) * unit('knots')
    assert_field(fc.CAS, CAS_truth)

    EAS_truth = array([582.1012, 317.5222]) * unit('knots')
    assert_field(fc.EAS, EAS_truth)

    mach_truth = array([0.88, 0.88]) * dimless
    assert_field(fc.mach, mach_truth)


def test_reynolds_number():
    """Test Reynolds number calculations. """

    ell = 5.34 * unit('ft')
    h_geom = 44.5 * unit('km')
    M_inf = 0.93 * dimless
    fc = FlightCondition(h_geom, mach=M_inf)

    Re_test = fc.reynolds_number(ell).magnitude
    Re_truth = 62278

    assert Re_test == myapprox(Re_truth)


def test_input_altitude_bounds():
    """Test that input altitude is properly bounded. Both FlightCondition
    and Atmosphere are covered in test since embedded Atmosphere object
    raises error.
    """

    M_inf = 0.44 * dimless
    atm = Atmosphere(0*unit('km'))

    h_below_min = atm._h_min*1.01
    with pytest.raises(ValueError) as e_info:
        FlightCondition(h_below_min, mach=M_inf)

    h_above_max = atm._h_max*1.01
    with pytest.raises(ValueError) as e_info:
        FlightCondition(h_above_max, mach=M_inf)


def test_mach_bounds():
    """Test that input is properly bounded. """

    h_geom = 13.37 * unit('km')
    fc = FlightCondition(h_geom, mach=0.42*dimless)

    M_below_min = fc._mach_min - (0.00001*dimless)
    with pytest.raises(ValueError) as e_info:
        FlightCondition(h_geom, mach=M_below_min)

    M_above_max = fc._mach_max*1.01
    with pytest.raises(ValueError) as e_info:
        FlightCondition(h_geom, mach=M_above_max)
