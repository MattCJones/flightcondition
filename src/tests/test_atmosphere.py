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

from numpy import array

from flightcondition import Atmosphere, unit
from common import assert_field

# Standard Atmospheric ground truth data
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

# NRL MSIS ground truth data
reltol_msis = 0.001
# Version 2.0:
datetime_msis2 = "2014-03-25T00:00"
f107_msis2 = 150
f107a_msis2 = 150
ap_msis2 = 7
h_msis2 = array([500.0, 1000.0]) * unit('km')
lon_msis2 = 15.0
lat_msis2 = 65.0
O_msis2 = array([2.363E+07, 1.776E+04]) * unit('1/cm^3')
N2_msis2 = array([5.129E+05, 1.738E+00]) * unit('1/cm^3')
O2_msis2 = array([8.620E+03, 4.869E-03]) * unit('1/cm^3')
#rho_msis2 = array([6.752E-16, 3.075E-18]) * unit('1/cm^3')
T_msis2 = array([1052.8, 1052.9]) * unit('degK')
He_msis2 = array([2.102E+06, 3.480E+05]) * unit('1/cm^3')
Ar_msis2 = array([1.919E+01, 3.039E-07]) * unit('1/cm^3')
H_msis2 = array([6.582E+04, 4.198E+04]) * unit('1/cm^3')
N_msis2 = array([3.529E+05, 6.497E+02]) * unit('1/cm^3')


@pytest.fixture
def atm():
    """Fixture to compute atmospheric properties just once before tests.

    Returns:
        Atmosphere object

    """
    return Atmosphere(h_geom_truth_arr)

def test_h_geop(atm):
    """Test geopotential altitude calculations. """
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

@pytest.fixture
def atm_msis2():
    """Fixture to compute MSIS 2.0 atmospheric properties just once before
    tests.

    Returns:
        Atmosphere object

    """
    return Atmosphere(h_msis2, datetime=datetime_msis2, lon=lon_msis2, 
				      lat=lat_msis2, f107=f107_msis2, f107a=f107a_msis2, ap=7,
                      model="msis2.0")

def test_O_msis2(atm_msis2):
    """Test MSIS 2.0 O. """
    assert_field(atm_msis2.species.O.to('1/cm^3'), O_msis2, reltol=reltol_msis)

def test_N2_msis2(atm_msis2):
    """Test MSIS 2.0 N2. """
    assert_field(atm_msis2.species.N2.to('1/cm^3'), N2_msis2,
                 reltol=reltol_msis)

def test_O2_msis2(atm_msis2):
    """Test MSIS 2.0 O2. """
    assert_field(atm_msis2.species.O2.to('1/cm^3'), O2_msis2,
                 reltol=reltol_msis)

def test_T_msis2(atm_msis2):
    """Test MSIS 2.0 T. """
    assert_field(atm_msis2.T.to('degK'), T_msis2, reltol=reltol_msis)

def test_He_msis2(atm_msis2):
    """Test MSIS 2.0 He. """
    assert_field(atm_msis2.species.He.to('1/cm^3'), He_msis2,
                 reltol=reltol_msis)

def test_Ar_msis2(atm_msis2):
    """Test MSIS 2.0 Ar. """
    assert_field(atm_msis2.species.Ar.to('1/cm^3'), Ar_msis2,
                 reltol=reltol_msis)

def test_H_msis2(atm_msis2):
    """Test MSIS 2.0 H. """
    assert_field(atm_msis2.species.H.to('1/cm^3'), H_msis2, reltol=reltol_msis)

def test_N_msis2(atm_msis2):
    """Test MSIS 2.0 N. """
    assert_field(atm_msis2.species.N.to('1/cm^3'), N_msis2, reltol=reltol_msis)
