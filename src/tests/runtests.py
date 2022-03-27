#!/usr/bin/env python
"""
Unit tests to verify code functionality.

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

import unittest
import sys

from numpy import array

from aeroutils import Atmosphere, FlightCondition, unit, dimless

# Atmospheric ground truth data
h_geom_truth_arr = [
        0,  5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000,
        50000, 55000, 60000, 65000, 70000, 75000, 80000] * unit('m')
h_geop_truth_arr = [
        0.0,            4996.07027357,  9984.29343877,  14964.68796877,
        19937.27227877, 24902.06472628, 29859.08361133, 34808.34717666,
        39749.87360801, 44683.68103426, 49609.78752775, 54528.2111044,
        59438.969724,   64342.08129041, 69237.56365177, 74125.4346007,
        79005.71187457] * unit('m')
T_inf_truth_arr = [
        288.15,       255.67554322, 223.25209265, 216.65,       216.65,
        221.55206473, 226.50908361, 236.51337209, 250.3496461,  264.1643069,
        270.65,       260.77100891, 247.02088477, 233.29217239, 219.58482178,
        208.3991308,  198.63857625] * unit('K')
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
        329.798731,   323.72379141, 315.0734446,  306.19285211, 297.06133141,
        289.39626128, 282.53793156] * unit('m/s')
nu_inf_truth_arr = [
        1.46071857e-05, 2.21100607e-05, 3.52509330e-05, 7.29951161e-05,
        1.59894147e-04, 3.61349481e-04, 8.01340459e-04, 1.80625312e-03,
        4.00667357e-03, 8.49961937e-03, 1.65908919e-02, 2.91173140e-02,
        5.11412252e-02, 9.26179078e-02, 1.73576137e-01, 3.44655513e-01,
        7.15580116e-01] * unit('m^2/s')

maxpercdiff = 0.01


def percdiff(arg1, arg2):
    """Percentage difference of arg1 relative to arg2 """
    return 100*(arg1 - arg2)/arg2


def ut_print(*args, **kwargs):
    """Print that takes into account verbosity level of test suite. """
    if ('-v' in sys.argv) or ('--verbose' in sys.argv):
        print(*args, **kwargs)


@unittest.skipIf(False, "Skipping for debug")
class TestStandardAtm(unittest.TestCase):
    """Run through test sets and verify output. """

    def setUp(self):
        """Set up fields and testing function. """

        ut_print("\nComputing standard atmosphere fields")
        atm = Atmosphere(h_geom_truth_arr)
        self.h_geop_arr = atm.H
        self.T_inf_arr = atm.T
        self.p_inf_arr = atm.p
        self.rho_inf_arr = atm.rho
        self.a_inf_arr = atm.a
        self.nu_inf_arr = atm.nu

        def test_field(field_arr, field_truth_arr, field_name, unit_str):
            """Test that output field matches truth data. """
            ut_print("\nTesting {field_name}")
            for h_geom_truth, field, field_truth in zip(
                    h_geom_truth_arr.to('km').magnitude,
                    field_arr.to(unit_str).magnitude,
                    field_truth_arr.to(unit_str).magnitude):
                pd = percdiff(field, field_truth)
                self.assertAlmostEqual(
                    0, pd, delta=maxpercdiff,
                    msg=f"Failed at z={h_geom_truth:.2g} km")
        self.test_field = test_field

    def test_h_geop(self):
        self.test_field(self.h_geop_arr[1:], h_geop_truth_arr[1:], 'h_geop',
                        'km')

    def test_T_inf(self):
        self.test_field(self.T_inf_arr, T_inf_truth_arr, 'T_inf', 'degC')

    def test_p_inf(self):
        self.test_field(self.p_inf_arr, p_inf_truth_arr, 'p_inf', 'kPa')

    def test_rho_inf(self):
        self.test_field(self.rho_inf_arr, rho_inf_truth_arr, 'rho_inf',
                        unit_str='kg/m^3')

    def test_a_inf(self):
        self.test_field(self.a_inf_arr, a_inf_truth_arr, 'a_inf', 'm/s')

    def test_nu_inf(self):
        self.test_field(self.nu_inf_arr, nu_inf_truth_arr, 'nu_inf', 'm^2/s')


@unittest.skipIf(False, "Skipping for debug")
class TestFlightCondition(unittest.TestCase):
    """Run through FlightCondition test sets and verify output. """

    def setUp(self):
        """Set up fields and testing function. """
        ut_print("\nSetting up FlightCondition test variables")
        self.h_geom_arr = [0, 30e3] * unit('ft')

        def test_field(field_arr, field_truth_arr, field_name, unit_str):
            """Test that output field matches truth data. """
            ut_print("\nTesting {field_name}")
            for h_geom, field, field_truth in zip(
                    self.h_geom_arr.to('kft').magnitude,
                    field_arr.to(unit_str).magnitude,
                    field_truth_arr.to(unit_str).magnitude):
                pd = percdiff(field, field_truth)
                self.assertAlmostEqual(0, pd, delta=maxpercdiff,
                                       msg=f"Failed at h={h_geom:.2g} kft")
        self.test_field = test_field

    def test_TAS(self):
        fc = FlightCondition(self.h_geom_arr, TAS=300*unit('knots'))
        TAS_truth = array([300, 300]) * unit('knots')
        self.test_field(fc.TAS, TAS_truth, 'KTAS', 'knots')
        CAS_truth = array([300, 187.7518]) * unit('knots')
        self.test_field(fc.CAS, CAS_truth, 'KCAS', 'knots')
        EAS_truth = array([300, 183.6448]) * unit('knots')
        self.test_field(fc.EAS, EAS_truth, 'KEAS', 'knots')
        mach_truth = array([0.4535, 0.5090]) * dimless
        self.test_field(fc.mach, mach_truth, 'mach', 'dimensionless')

    def test_CAS(self):
        fc = FlightCondition(self.h_geom_arr, CAS=300*unit('knots'))
        TAS_truth = array([300, 465.6309]) * unit('knots')
        self.test_field(fc.TAS, TAS_truth, 'KTAS', 'knots')
        CAS_truth = array([300, 300]) * unit('knots')
        self.test_field(fc.CAS, CAS_truth, 'KCAS', 'knots')
        EAS_truth = array([300, 285.0357]) * unit('knots')
        self.test_field(fc.EAS, EAS_truth, 'KEAS', 'knots')
        mach_truth = array([0.4535, 0.7900]) * dimless
        self.test_field(fc.mach, mach_truth, 'mach', 'dimensionless')

    def test_EAS(self):
        fc = FlightCondition(self.h_geom_arr, EAS=300*unit('knots'))
        TAS_truth = array([300, 490.0764]) * unit('knots')
        self.test_field(fc.TAS, TAS_truth, 'KTAS', 'knots')
        CAS_truth = array([300, 317.3602]) * unit('knots')
        self.test_field(fc.CAS, CAS_truth, 'KCAS', 'knots')
        EAS_truth = array([300, 300]) * unit('knots')
        self.test_field(fc.EAS, EAS_truth, 'KEAS', 'knots')
        mach_truth = array([0.4535, 0.8314]) * dimless
        self.test_field(fc.mach, mach_truth, 'mach', 'dimensionless')

    def test_mach(self):
        fc = FlightCondition(self.h_geom_arr, mach=0.88*dimless)
        TAS_truth = array([582.1012, 518.7004]) * unit('knots')
        self.test_field(fc.TAS, TAS_truth, 'KTAS', 'knots')
        CAS_truth = array([582.1012, 337.977]) * unit('knots')
        self.test_field(fc.CAS, CAS_truth, 'KCAS', 'knots')
        EAS_truth = array([582.1012, 317.5222]) * unit('knots')
        self.test_field(fc.EAS, EAS_truth, 'KEAS', 'knots')
        mach_truth = array([0.88, 0.88]) * dimless
        self.test_field(fc.mach, mach_truth, 'mach', 'dimensionless')

    def test_reynolds_number(self):
        ell = 5.34 * unit('ft')
        h_geom = 44.5 * unit('km')
        M_inf = 0.93 * dimless
        fc = FlightCondition(h_geom, mach=M_inf)

        Re_test = fc.reynolds_number(ell).magnitude
        Re_truth = 62278
        self.assertAlmostEqual(Re_test, Re_truth,
                               delta=1,
                               msg="Reynolds number failed.")

        ell_unit = 'in'
        Re_by_ell_test = fc.reynolds_number_per_unit_length(ell_unit).magnitude
        Re_by_ell_truth = 971.88
        self.assertAlmostEqual(Re_by_ell_test/1e3, Re_by_ell_truth/1e3,
                               delta=0.01,
                               msg="Reynolds number per length failed.")

    def test_input_altitude_bounds(self):
        """Test that input altitude is properly bounded. Both FlightCondition
        and Atmosphere are covered in test since embedded Atmosphere object
        raises error.
        """
        M_inf = 0.44 * dimless
        atm = Atmosphere(0*unit('km'))

        h_below_min = atm._h_min*1.01
        with self.assertRaises(ValueError):
            FlightCondition(h_below_min, mach=M_inf)

        h_above_max = atm._h_max*1.01
        with self.assertRaises(ValueError):
            FlightCondition(h_above_max, mach=M_inf)

    def test_mach_bounds(self):
        """Test that input is properly bounded. """
        h_geom = 13.37 * unit('km')
        fc = FlightCondition(h_geom, mach=0.42*dimless)

        M_below_min = fc._mach_min - (0.00001*dimless)
        with self.assertRaises(ValueError):
            FlightCondition(h_geom, mach=M_below_min)

        M_above_max = fc._mach_max*1.01
        with self.assertRaises(ValueError):
            FlightCondition(h_geom, mach=M_above_max)


if __name__ == '__main__':
    unittest.main()
