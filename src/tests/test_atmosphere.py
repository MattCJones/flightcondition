#!/usr/bin/env python
"""
Test atmospheric capabilities.

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

# flake8: noqa E203

import pytest

from flightcondition import Atmosphere, unit
from common import assert_field

# Atmospheric ground truth data
h_geom_truth_arr = [
    0, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000,
    50000, 55000, 60000, 65000, 70000, 75000, 80000] * unit('m')
h_geop_truth_arr = [
    0.0           , 4996.07027357 , 9984.29343877 , 14964.68796877,
    19937.27227877, 24902.06472628, 29859.08361133, 34808.34717666,
    39749.87360801, 44683.68103426, 49609.78752775, 54528.2111044 ,
    59438.969724  , 64342.08129041, 69237.56365177, 74125.4346007 ,
    79005.71187457] * unit('m')
T_inf_truth_arr = [
    288.15      , 255.67554322, 223.25209265, 216.65      , 216.65      ,
    221.55206473, 226.50908361, 236.51337209, 250.3496461 , 264.1643069 ,
    270.65      , 260.77100891, 247.02088477, 233.29217239, 219.58482178,
    208.3991308 , 198.63857625] * unit('K')
p_inf_truth_arr = [
    1.01325000e+05, 5.40482622e+04, 2.64998731e+04, 1.21117861e+04,
    5.52929078e+03, 2.54921293e+03, 1.19702628e+03, 5.74591263e+02,
    2.87142182e+02, 1.49100428e+02, 7.97788547e+01, 4.25248048e+01,
    2.19584937e+01, 1.09296242e+01, 5.22085015e+00, 2.38812369e+00,
    1.05246447e+00] * unit('Pa')
rho_inf_truth_arr = [
    1.22500002e+00, 7.36428613e-01, 4.13510330e-01, 1.94754547e-01,
    8.89096382e-02, 4.00837567e-02, 1.84101009e-02, 8.46333291e-03,
    3.99565628e-03, 1.96626868e-03, 1.02687569e-03, 5.68095211e-04,
    3.09675594e-04, 1.63208648e-04, 8.28279701e-05, 3.99207802e-05,
    1.84578859e-05] * unit('kg/m^3')
a_inf_truth_arr = [
    340.29398803, 320.54540686, 299.53166026, 295.06949351, 295.06949351,
    298.38903875, 301.70866004, 308.29949587, 317.18924664, 325.82322113,
    329.798731  , 323.72379141, 315.0734446 , 306.19285211, 297.06133141,
    289.39626128, 282.53793156] * unit('m/s')
nu_inf_truth_arr = [
    1.46071857e-05, 2.21100607e-05, 3.52509330e-05, 7.29951161e-05,
    1.59894147e-04, 3.61349481e-04, 8.01340459e-04, 1.80625312e-03,
    4.00667357e-03, 8.49961937e-03, 1.65908919e-02, 2.91173140e-02,
    5.11412252e-02, 9.26179078e-02, 1.73576137e-01, 3.44655513e-01,
    7.15580116e-01] * unit('m^2/s')


@pytest.fixture
def atm():
    """Fixture to compute atmospheric properties just once before tests.
    :returns: Atmosphere object

    """
    return Atmosphere(h_geom_truth_arr)

def test_h_geop(atm):
    """Test geopential altitude calculations. """
    assert_field(atm.H, h_geop_truth_arr)

def test_T_inf(atm):
    """Test temperature calculations. """
    assert_field(atm.T, T_inf_truth_arr)

def test_p_inf(atm):
    """Test pressure calculations. """
    assert_field(atm.p, p_inf_truth_arr)

def test_rho_inf(atm):
    """Test density calculations. """
    assert_field(atm.rho, rho_inf_truth_arr)

def test_a_inf(atm):
    """Test sound speed calculations. """
    assert_field(atm.a, a_inf_truth_arr)

def test_nu_inf(atm):
    """Test kinematic viscosity calculations. """
    assert_field(atm.nu, nu_inf_truth_arr)
