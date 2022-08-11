#!/usr/bin/env python
"""
Test flight condition functionality.

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

# flake8: noqa E203

import warnings
import re

from shlex import split
from subprocess import run

import pytest

from numpy import array

from flightcondition import Atmosphere, FlightCondition, unit, dimless
from common import assert_field, myapprox

h_geom_arr = [0, 30e3] * unit('ft')


def test_TAS():
    """Test true airspeed calculations. """

    fc = FlightCondition(h_geom_arr, TAS=300*unit('knots'))

    TAS_truth = array([300, 300]) * unit('knots')
    assert_field(fc.vel.TAS, TAS_truth)

    CAS_truth = array([300, 187.7518]) * unit('knots')
    assert_field(fc.vel.CAS, CAS_truth)

    EAS_truth = array([300, 183.6448]) * unit('knots')
    assert_field(fc.vel.EAS, EAS_truth)

    M_truth = array([0.4535, 0.5090]) * dimless
    assert_field(fc.vel.M, M_truth)


def test_CAS():
    """Test calibrated airspeed calculations. """

    fc = FlightCondition(h_geom_arr, CAS=300*unit('knots'))

    TAS_truth = array([300, 465.6309]) * unit('knots')
    assert_field(fc.vel.TAS, TAS_truth)

    CAS_truth = array([300, 300]) * unit('knots')
    assert_field(fc.vel.CAS, CAS_truth)

    EAS_truth = array([300, 285.0357]) * unit('knots')
    assert_field(fc.vel.EAS, EAS_truth)

    M_truth = array([0.4535, 0.7900]) * dimless
    assert_field(fc.vel.M, M_truth)


def test_EAS():
    """Test equivalent airspeed calculations. """

    fc = FlightCondition(h_geom_arr, EAS=300*unit('knots'))

    TAS_truth = array([300, 490.0764]) * unit('knots')
    assert_field(fc.vel.TAS, TAS_truth)

    CAS_truth = array([300, 317.3602]) * unit('knots')
    assert_field(fc.vel.CAS, CAS_truth)

    EAS_truth = array([300, 300]) * unit('knots')
    assert_field(fc.vel.EAS, EAS_truth)

    M_truth = array([0.4535, 0.8314]) * dimless
    assert_field(fc.vel.M, M_truth)


def test_mach():
    """Test Mach number calculations. """

    fc = FlightCondition(h_geom_arr, M=0.88*dimless)

    TAS_truth = array([582.1012, 518.7004]) * unit('knots')
    assert_field(fc.vel.TAS, TAS_truth)

    CAS_truth = array([582.1012, 337.977]) * unit('knots')
    assert_field(fc.vel.CAS, CAS_truth)

    EAS_truth = array([582.1012, 317.5222]) * unit('knots')
    assert_field(fc.vel.EAS, EAS_truth)

    M_truth = array([0.88, 0.88]) * dimless
    assert_field(fc.vel.M, M_truth)


def test_other_vel_properties():
    """Test other velocity property calculations. """

    fc = FlightCondition(h_geom_arr, TAS=300*unit('knots'), unit_system='US')

    q_c_truth = array([320.6898, 121.7655]) * unit('lbf/ft^2')
    assert_field(fc.vel.q_c, q_c_truth)

    p0_truth = array([2436.9064, 751.4328]) * unit('lbf/ft^2')
    assert_field(fc.vel.p0, p0_truth)

    T0_truth = array([540.0069, 433.1758]) * unit('degR')
    assert_field(fc.vel.T0, T0_truth.to('degK'))

    Tr_lamr_truth = array([536.8064, 429.9753]) * unit('degR')
    assert_field(fc.vel.Tr_lamr, Tr_lamr_truth.to('degK'), reltol=0.001)

    Tr_turb_truth = array([537.6599, 430.8287]) * unit('degR')
    assert_field(fc.vel.Tr_turb, Tr_turb_truth.to('degK'), reltol=0.001)

    Re_by_L_truth = array([2.6837e+5, 1.2096e+5]) * unit('1/in')
    assert_field(fc.vel.Re_by_L, Re_by_L_truth)


def test_reynolds_number():
    """Test Reynolds number calculations. """

    L = 5.34 * unit('ft')
    h_geom = 44.5 * unit('km')
    M_ = 0.93 * dimless
    fc = FlightCondition(h_geom, M=M_, L=L)

    Re_test = fc.len.Re.magnitude
    Re_truth = 62278

    assert Re_test == myapprox(Re_truth)


def test_access_byname():
    """Test that quantities are properly accessible by name. """

    L = 5.34 * unit('ft')
    h_geom = 44.5 * unit('km')
    M_ = 0.93 * dimless
    fc = FlightCondition(h_geom, M=M_, L=L, unit_system='US')

    # Check that sub-objects .byname works properly
    assert fc.atm.p == fc.atm.byname.pressure
    assert fc.atm.T == fc.atm.byname.temperature
    assert fc.atm.rho == fc.atm.byname.density
    assert fc.atm.nu == fc.atm.byname.kinematic_viscosity

    assert fc.vel.M == fc.vel.byname.mach_number
    assert fc.vel.TAS == fc.vel.byname.true_airspeed
    assert fc.vel.CAS == fc.vel.byname.calibrated_airspeed
    assert fc.vel.EAS == fc.vel.byname.equivalent_airspeed

    assert fc.len.Re == fc.len.byname.reynolds_number

    # # Check that base object .byname works properly
    # assert fc.byname.pressure == fc.atm.byname.pressure
    # assert fc.byname.temperature == fc.atm.byname.temperature
    # assert fc.byname.density == fc.atm.byname.density
    # assert fc.byname.kinematic_viscosity == fc.atm.byname.kinematic_viscosity

    # assert fc.byname.mach_number == fc.vel.byname.mach_number
    # assert fc.byname.true_airspeed == fc.vel.byname.true_airspeed
    # assert fc.byname.calibrated_airspeed == fc.vel.byname.calibrated_airspeed
    # assert fc.byname.equivalent_airspeed == fc.vel.byname.equivalent_airspeed

    # assert fc.byname.reynolds_number == fc.len.byname.reynolds_number


def test_input_altitude_bounds():
    """Test that input altitude is properly bounded. Both FlightCondition
    and Atmosphere are covered in test since embedded Atmosphere object
    raises error.
    """

    M_ = 0.44 * dimless
    atm = Atmosphere(0*unit('km'))

    h_below_min = atm._h_min*1.01
    with pytest.raises(ValueError) as e_info:
        FlightCondition(h_below_min, M=M_)

    h_above_max = atm._h_max*1.01
    with pytest.raises(ValueError) as e_info:
        FlightCondition(h_above_max, M=M_)


def test_mach_bounds():
    """Test that input is properly bounded. """

    h_geom = 13.37 * unit('km')
    fc = FlightCondition(h_geom, M=0.42*dimless)

    M_below_min = fc._mach_min - (0.00001*dimless)
    with warnings.catch_warnings():  # catch warning from sqrt
        warnings.simplefilter("ignore")
        with pytest.raises(ValueError) as e_info:
            FlightCondition(h_geom, M=M_below_min)

    M_above_max = fc._mach_max*1.01
    with pytest.raises(ValueError) as e_info:
        FlightCondition(h_geom, M=M_above_max)


def test_command_line_interface():
    """Test that command line interface is running properly. """
    cmd_str = "flightcondition --alt 23 kft --EAS 233 kt --len 4 ft"
    out = run(split(cmd_str), capture_output=True)
    out_str = out.stdout.decode()
    out_regex = r"""[=]+
\s+Flight Condition.*
[=]+
[-]+\s+Altitude Quantities\s+[-]+
.*h\s+= 23 \w+
.*
[-]+\s+Airspeed Quantities\s+[-]+
.*EAS\s+= 233 kt
.*
[-]+\s+Length Quantities\s+[-]+
.*L\s+= 4 ft
.*"""
    re_out = re.search(out_regex, out_str, re.DOTALL)
    # Helpful print statements for debugging failure:
    print("re_out:\n", re_out)
    print("out_str:\n", out_str)
    print("out_regex:\n", out_regex)
    assert re_out
