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
from flightcondition.boundarylayer import BoundaryLayer
from flightcondition.nondimensional import NonDimensional
from common import assert_field, myapprox

h_geom_arr = [0, 30e3] * unit('ft')


def test_TAS():
    """Test true airspeed calculations. """

    fc = FlightCondition(h_geom_arr, TAS=300*unit('knots'))

    TAS_truth = array([300, 300]) * unit('knots')
    assert_field(fc.byvel.TAS, TAS_truth)

    CAS_truth = array([300, 187.7518]) * unit('knots')
    assert_field(fc.byvel.CAS, CAS_truth)

    EAS_truth = array([300, 183.6448]) * unit('knots')
    assert_field(fc.byvel.EAS, EAS_truth)

    M_truth = array([0.4535, 0.5090]) * dimless
    assert_field(fc.byvel.M, M_truth)


def test_CAS():
    """Test calibrated airspeed calculations. """

    fc = FlightCondition(h_geom_arr, CAS=300*unit('knots'))

    TAS_truth = array([300, 465.6309]) * unit('knots')
    assert_field(fc.byvel.TAS, TAS_truth)

    CAS_truth = array([300, 300]) * unit('knots')
    assert_field(fc.byvel.CAS, CAS_truth)

    EAS_truth = array([300, 285.0357]) * unit('knots')
    assert_field(fc.byvel.EAS, EAS_truth)

    M_truth = array([0.4535, 0.7900]) * dimless
    assert_field(fc.byvel.M, M_truth)


def test_EAS():
    """Test equivalent airspeed calculations. """

    fc = FlightCondition(h_geom_arr, EAS=300*unit('knots'))

    TAS_truth = array([300, 490.0764]) * unit('knots')
    assert_field(fc.byvel.TAS, TAS_truth)

    CAS_truth = array([300, 317.3602]) * unit('knots')
    assert_field(fc.byvel.CAS, CAS_truth)

    EAS_truth = array([300, 300]) * unit('knots')
    assert_field(fc.byvel.EAS, EAS_truth)

    M_truth = array([0.4535, 0.8314]) * dimless
    assert_field(fc.byvel.M, M_truth)


def test_mach():
    """Test Mach number calculations. """

    fc = FlightCondition(h_geom_arr, M=0.88*dimless)

    TAS_truth = array([582.1012, 518.7004]) * unit('knots')
    assert_field(fc.byvel.TAS, TAS_truth)

    CAS_truth = array([582.1012, 337.977]) * unit('knots')
    assert_field(fc.byvel.CAS, CAS_truth)

    EAS_truth = array([582.1012, 317.5222]) * unit('knots')
    assert_field(fc.byvel.EAS, EAS_truth)

    M_truth = array([0.88, 0.88]) * dimless
    assert_field(fc.byvel.M, M_truth)


def test_other_vel_properties():
    """Test other velocity property calculations. """

    fc = FlightCondition(h_geom_arr, TAS=300*unit('knots'), units='US')

    q_c_truth = array([320.6898, 121.7655]) * unit('lbf/ft^2')
    assert_field(fc.byvel.q_c, q_c_truth)

    p0_truth = array([2436.9064, 751.4328]) * unit('lbf/ft^2')
    assert_field(fc.byvel.p0, p0_truth)

    T0_truth = array([540.0069, 433.1758]) * unit('degR')
    assert_field(fc.byvel.T0, T0_truth.to('degK'))

    Tr_lamr_truth = array([536.8064, 429.9753]) * unit('degR')
    assert_field(fc.byvel.Tr_lamr, Tr_lamr_truth.to('degK'), reltol=0.001)

    Tr_turb_truth = array([537.6599, 430.8287]) * unit('degR')
    assert_field(fc.byvel.Tr_turb, Tr_turb_truth.to('degK'), reltol=0.001)

    Re_by_L_truth = array([2.6837e+5, 1.2096e+5]) * unit('1/in')
    assert_field(fc.byvel.Re_by_L, Re_by_L_truth)


def test_other_len_properties():
    """Test other length property calculations. """

    fc = FlightCondition(h_geom_arr, TAS=300*unit('knots'), L=1*unit.ft,
                         units='US')

    Re_truth = array([3.2204e+6, 1.4516e+6]) * dimless
    assert_field(fc.bylen.Re, Re_truth)

    h_BL_lamr_truth = array([0.0348, 0.0518]) * unit('in')
    assert_field(fc.bylen.h_BL_lamr.to('in'), h_BL_lamr_truth, reltol=0.04)

    h_BL_turb_truth = array([0.2245, 0.2633]) * unit('in')
    assert_field(fc.bylen.h_BL_turb.to('in'), h_BL_turb_truth, reltol=0.001)

    Cf_lamr_truth = array([0.00074, 0.0011]) * dimless
    assert_field(fc.bylen.Cf_lamr, Cf_lamr_truth, reltol=0.005)

    Cf_turb_truth = array([0.0036, 0.0042]) * dimless
    assert_field(fc.bylen.Cf_turb, Cf_turb_truth, reltol=0.05)

    # Test all skin friction methods
    x = fc.bylen.L
    Re_x = NonDimensional.reynolds_number(U=fc.byvel.TAS, L=x, nu=fc.byalt.nu)

    Cf_turb = BoundaryLayer.flat_plate_skin_friction_coeff_turb(
        Re_x=Re_x, source='calibrated_powerlaw')
    assert_field(fc.bylen.Cf_turb, Cf_turb_truth, reltol=0.05)

    Cf_turb = BoundaryLayer.flat_plate_skin_friction_coeff_turb(
        Re_x=Re_x, source='schlichting')
    assert_field(fc.bylen.Cf_turb, Cf_turb_truth, reltol=0.05)

    Cf_turb = BoundaryLayer.flat_plate_skin_friction_coeff_turb(
        Re_x=Re_x, source='prandtl1927')
    assert_field(fc.bylen.Cf_turb, Cf_turb_truth, reltol=0.05)

    Cf_turb = BoundaryLayer.flat_plate_skin_friction_coeff_turb(
        Re_x=Re_x, source='prandtl_schlichting1932')
    assert_field(fc.bylen.Cf_turb, Cf_turb_truth, reltol=0.05)

    Cf_turb = BoundaryLayer.flat_plate_skin_friction_coeff_turb(
        Re_x=Re_x, source='schultz_grunov1940')
    assert_field(fc.bylen.Cf_turb, Cf_turb_truth, reltol=0.05)

    Cf_turb = BoundaryLayer.flat_plate_skin_friction_coeff_turb(
        Re_x=Re_x, source='ittc1957')
    assert_field(fc.bylen.Cf_turb, Cf_turb_truth, reltol=0.05)

    Cf_turb = BoundaryLayer.flat_plate_skin_friction_coeff_turb(
        Re_x=Re_x, source='granville1977')
    assert_field(fc.bylen.Cf_turb, Cf_turb_truth, reltol=0.05)

    Cf_turb = BoundaryLayer.flat_plate_skin_friction_coeff_turb(
        Re_x=Re_x, source='roskam1987')
    assert_field(fc.bylen.Cf_turb, Cf_turb_truth, reltol=0.05)

    Cf_turb = BoundaryLayer.flat_plate_skin_friction_coeff_turb(
        Re_x=Re_x, source='white2006')
    assert_field(fc.bylen.Cf_turb, Cf_turb_truth, reltol=0.05)

    # Test yplus computations
    wall_distance_from_yplus1 = fc.bylen.wall_distance_from_yplus(1)
    wall_distance_from_yplus1.ito('m')
    wall_distance_from_yplus1_truth = [2.4e-6, 5.0e-6] * unit.m
    assert_field(wall_distance_from_yplus1, wall_distance_from_yplus1_truth,
                 reltol=0.10)
    assert_field(wall_distance_from_yplus1, fc.bylen.h_yplus1)

    wall_distance_from_yplus30 = fc.bylen.wall_distance_from_yplus(30)
    wall_distance_from_yplus30.ito('m')
    wall_distance_from_yplus30_truth = [7.2e-5, 1.5e-4] * unit.m
    assert_field(wall_distance_from_yplus30, wall_distance_from_yplus30_truth,
                 reltol=0.10)

    wall_distance_from_yplus100 = fc.bylen.wall_distance_from_yplus(100)
    wall_distance_from_yplus100.ito('m')
    wall_distance_from_yplus100_truth = [2.4e-4, 5.0e-4] * unit.m
    assert_field(wall_distance_from_yplus100, wall_distance_from_yplus100_truth,
                 reltol=0.10)


def test_reynolds_number():
    """Test Reynolds number calculations. """

    L = 5.34 * unit('ft')
    h_geom = 44.5 * unit('km')
    M_ = 0.93 * dimless
    fc = FlightCondition(h_geom, M=M_, L=L)

    Re_test = fc.bylen.Re.magnitude
    Re_truth = 62278

    assert Re_test == myapprox(Re_truth)


def test_reynolds_number_velocity():
    """Test Reynolds number calculations. """

    L = 5.34 * unit('ft')
    h_geom = 44.5 * unit('km')
    M_ = 0.93 * dimless
    fc = FlightCondition(h_geom, M=M_, L=L)

    # Set velocity based on Reynolds number
    Re_arr = [1e4, 1e5, 1e6]
    fc.bylen.Re = Re_arr
    M_truth = [0.14933, 1.4933, 14.933] * dimless
    assert_field(fc.byvel.M, M_truth)

    # Check setting Re at initlization
    fc2 = FlightCondition(h_geom, L=L, Re=Re_arr)
    assert_field(fc.byvel.TAS, fc2.byvel.TAS)


def test_holding_velocity_constant():
    """Test that proper velocity quantity is held constant when altitude is
    dynamically changed. """
    fc = FlightCondition(h=[4, 7, 44]*unit.kft, M=0.5)

    # Dynamically change velocity --> Mach number should stay constant
    M = fc.M
    fc.h = [2, 13, 24] * unit.kft
    assert_field(M, fc.M)

    # Dynamically change velocity --> true airspeed should stay constant
    fc.TAS = [100, 200, 300] * unit('m/s')
    TAS = fc.TAS
    fc.h = [1, 9, 33] * unit.kft
    assert_field(TAS, fc.TAS)

    # Dynamically change velocity --> calibrated airspeed should stay constant
    fc.CAS = [110, 220, 330] * unit('m/s')
    CAS = fc.CAS
    fc.h = [11, 19, 25] * unit.kft
    assert_field(CAS, fc.CAS)

    # Dynamically change velocity --> calibrated airspeed should stay constant
    fc.EAS = [111, 222, 333] * unit('m/s')
    EAS = fc.EAS
    fc.h = [7, 15, 21] * unit.kft
    assert_field(EAS, fc.EAS)


def test_access_byname():
    """Test that quantities are properly accessible by name. """

    L = 5.34 * unit('ft')
    h_geom = 44.5 * unit('km')
    M_ = 0.93 * dimless
    fc = FlightCondition(h_geom, M=M_, L=L, units='US')

    # Check that sub-objects .byname works properly
    assert fc.byalt.p == fc.byalt.byname.pressure
    assert fc.byalt.T == fc.byalt.byname.temperature
    assert fc.byalt.rho == fc.byalt.byname.density
    assert fc.byalt.nu == fc.byalt.byname.kinematic_viscosity

    assert fc.byvel.M == fc.byvel.byname.mach_number
    assert fc.byvel.TAS == fc.byvel.byname.true_airspeed
    assert fc.byvel.CAS == fc.byvel.byname.calibrated_airspeed
    assert fc.byvel.EAS == fc.byvel.byname.equivalent_airspeed

    assert fc.bylen.Re == fc.bylen.byname.reynolds_number

    # Check that base object .byname works properly
    assert fc.byname.pressure == fc.byalt.byname.pressure
    assert fc.byname.temperature == fc.byalt.byname.temperature
    assert fc.byname.density == fc.byalt.byname.density
    assert fc.byname.kinematic_viscosity == fc.byalt.byname.kinematic_viscosity

    assert fc.byname.mach_number == fc.byvel.byname.mach_number
    assert fc.byname.true_airspeed == fc.byvel.byname.true_airspeed
    assert fc.byname.calibrated_airspeed == fc.byvel.byname.calibrated_airspeed
    assert fc.byname.equivalent_airspeed == fc.byvel.byname.equivalent_airspeed

    assert fc.byname.reynolds_number == fc.bylen.byname.reynolds_number


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

def test_dynamic_array_warnings_and_errors():
    """Test that proper warnings and errors are thrown when dynamically
    changing property arrays. """

    mismatch_substr = "Non-singular arrays must be equal in size"

    # Setting altitude with of incorrect size should throw a warning
    with warnings.catch_warnings(record=True) as wrn:
        fc = FlightCondition()
        fc.M = [0.2, 0.5, 0.7]
        fc.h = [2, 3] * unit.km;
        assert mismatch_substr in str(wrn[-1].message)

    # Setting Mach number with of incorrect size should throw an error
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            fc = FlightCondition()
            fc.h = [2, 3] * unit.km;
            fc.M = [0.2, 0.5, 0.7]
        except AttributeError as err:
            assert (mismatch_substr in str(err))

    # Setting TAS with of incorrect size should throw an error
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            fc = FlightCondition()
            fc.h = [2, 3] * unit.km;
            fc.TAS = [120, 300, 540] * unit('m/s')
        except AttributeError as err:
            assert (mismatch_substr in str(err))

    # Setting CAS with of incorrect size should throw an error
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            fc = FlightCondition()
            fc.h = [2, 3] * unit.km;
            fc.CAS = [120, 300, 540] * unit('m/s')
        except AttributeError as err:
            assert (mismatch_substr in str(err))

    # Setting EAS with of incorrect size should throw an error
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            fc = FlightCondition()
            fc.h = [2, 3] * unit.km;
            fc.EAS = [120, 300, 540] * unit('m/s')
        except AttributeError as err:
            assert (mismatch_substr in str(err))

    # Setting Reynolds number of incorrect size should throw an error
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            fc = FlightCondition()
            fc.M = [0.2, 0.5, 0.7]
            fc.Re = [1e6, 2e5];
        except AttributeError as err:
            assert (mismatch_substr in str(err))

def test_command_line_interface():
    """Test that command line interface is running properly. """
    cmd_str = "flightcondition --h 23 kft --EAS 233 kt --L 4 ft"
    out = run(split(cmd_str), capture_output=True)
    out_str = out.stdout.decode()
    out_regex = r"""[=]+
\s+Flight Condition.*
[=]+
[-]+\s+Altitude Quantities\s+[-]+
.*h\s+= 23 \w+
.*
[-]+\s+Velocity Quantities\s+[-]+
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
