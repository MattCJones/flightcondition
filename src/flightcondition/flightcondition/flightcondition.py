#!/usr/bin/env python
"""Easily convert between Mach number, true airspeed (TAS), calibrated airspeed
(CAS), and equivalent airspeed (EAS) for given altitude(s).  Additional flight
condition data and atmospheric data is computed.

Dependencies: numpy, pint

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

from numpy import atleast_1d, ones, sqrt, shape

from flightcondition.atmosphere import Atmosphere
from flightcondition.atmosphere.atmosphere import _atleast_1d  # DEPRECATED
from flightcondition.constants import PhysicalConstants as Phys
from flightcondition.units import unit, dimless
from flightcondition.units.units import check_dimensioned,\
    check_area_dimensioned, check_length_dimensioned, to_base_units_wrapper


class FlightCondition:

    """Easily convert between Mach number, true airspeed (TAS), calibrated
    airspeed (CAS), and equivalent airspeed (EAS) for given altitude(s).
    Additional flight condition data and atmospheric data is computed.

    All inputs must be dimensional unit quantities.

    Usage:

        from flightcondition import FlightCondition, unit, dimless

        # Compute flight conditions for a scalar or array of altitudes
        altitudes = [0, 10e3, 33.5e3] * unit('ft')
        fc = FlightCondition(altitudes, EAS=300*unit('knots'))

        # Print flight condition data:
        print(f"{fc}")

        # Print extended output in Imperial units:
        print(f"\n{fc.tostring(full_output=True, imperial_units=True)}")

        # Access flight speed formats individually
        M_inf, U_inf, U_CAS, U_EAS = fc.mach, fc.TAS, fc.CAS, fc.EAS

        # Access atmospheric data alone (see Atmosphere class for options)
        atm = fc.atm  # access Atmosphere object 'atm'
        p, T, rho, nu, a = atm.p, atm.T, atm.rho, atm.nu, atm.a

        # Input true/calibrated/equivalent airspeed or Mach number
        fc_TAS = FlightCondition(altitudes, TAS=300*unit('knots'))
        fc_CAS = FlightCondition(altitudes, CAS=300*unit('knots'))
        fc_EAS = FlightCondition(altitudes, EAS=300*unit('knots'))
        fc_mach = FlightCondition(altitudes, mach=0.4535*dimless)

        # Specify desired units on input and output
        altitudes_in_km = [0, 3.048, 10.2108] * unit('km')
        fc_other_units = FlightCondition(altitudes, EAS=154.33*unit('m/s'))
        U_TAS = fc_other_units.TAS
        print(f"\nThe true airspeed in m/s is {U_TAS.to('m/s'):.5g}")
        print(f"The true airspeed in km/s is {U_TAS.to('km/s'):.5g}")

        # Compute additional derived quantities (see class for all options)
        print(f"\nThe dynamic pressure in psi is {fc.q_inf.to('psi'):.5g}")
        ell = 60 * unit('in')  # arbitrary length scale of interest
        print(f"The Reynolds number is {fc.reynolds_number(ell):.5g}")
        print(f"The Reynolds number per-unit-length [1/in] is "
            f"{fc.reynolds_number_per_unit_length('in'):.5g}")
    """

    def __init__(self, altitude, mach=None, TAS=None, CAS=None, EAS=None,):
        """Constructor based on altitude and input speed in terms of Mach
        number, TAS, CAS, or EAS.  At least one input speed is required.  All
        inputs must be dimensional unit quantities.

        :altitude: geometric altitude
        :mach: Mach number
        :TAS: true airspeed
        :CAS: calibrated airspeed
        :EAS: equivalent airspeed

        """

        # Automatically process altitude through Atmosphere class
        self.atm = Atmosphere(altitude)
        self.altitude = self.atm.h

        h0 = 0 * unit('kft')
        self._atm0 = Atmosphere(h0)

        p_inf = self.atm.p
        p_inf_h0 = self._atm0.p
        self.delta = p_inf/p_inf_h0

        if mach is not None:
            self.mach = self.checkandsize(mach)
            self.TAS = self.TAS_from_M_inf(self.mach)
            self.EAS = self.EAS_from_TAS(self.TAS, self.mach)
            self.q_c = self.q_c_from_M_inf(self.mach)
            self.CAS = self.CAS_from_q_c(self.q_c)
        elif TAS is not None:
            self.TAS = self.checkandsize(TAS)
            self.mach = self.M_inf_from_TAS(TAS)
            self.EAS = self.EAS_from_TAS(self.TAS, self.mach)
            self.q_c = self.q_c_from_M_inf(self.mach)
            self.CAS = self.CAS_from_q_c(self.q_c)
        elif CAS is not None:
            self.CAS = self.checkandsize(CAS)
            self.q_c = self.q_c_from_CAS(self.CAS)
            self.mach = self.M_inf_from_q_c(self.q_c)
            self.TAS = self.TAS_from_M_inf(self.mach)
            self.EAS = self.EAS_from_TAS(self.TAS, self.mach)
        elif EAS is not None:
            self.EAS = self.checkandsize(EAS)
            self.mach = self.M_inf_from_EAS(self.EAS)
            self.TAS = self.TAS_from_M_inf(self.mach)
            self.q_c = self.q_c_from_M_inf(self.mach)
            self.CAS = self.CAS_from_q_c(self.q_c)
        else:
            raise TypeError("Input mach, TAS, CAS, or EAS")
        self.q_inf = self.q_inf_from_TAS(self.TAS)

        # Check that computations are valid
        M = _atleast_1d(self.mach)
        self._mach_min = 0 * dimless
        self._mach_max = 1 * dimless
        if (M < self._mach_min).any() or (self._mach_max < M).any():
            raise ValueError(
                f"Mach number is out of bounds "
                f"({self._mach_min:.5g} < mach < {self._mach_max:.5g})"
                )

    def __repr__(self):
        """Output string representation of class object. Default for __str__
        :returns: string output

        """
        return self.tostring(full_output=False, imperial_units=False)

    @staticmethod
    def isentropic_mach(p_0, p):
        """Isentropic flow equation for Mach number

        :p_0: stagnation pressure
        :p: static pressure
        :returns: Mach number

        """
        y = Phys.gamma_air
        M = sqrt((2/(y-1))*((p_0/p + 1)**((y-1)/y) - 1))
        return M

    @staticmethod
    def isentropic_stagnation_pressure(M, p):
        """Isentropic flow equation for stagnation pressure

        :M: Mach number
        :p: static pressure
        :returns: stagnation pressure

        """
        y = Phys.gamma_air
        p_0 = p*(-1 + (1 + ((y-1)/2)*M**2)**(y/(y-1)))
        return p_0

    @staticmethod
    def mach_number(U, a):
        """Compute Mach number

        :U: airspeed
        :a: sound speed
        :returns: Mach number

        """
        return U/a

    @staticmethod
    def mach_airspeed(M, a):
        """Compute airspeed from Mach number

        :M: Mach number
        :a: sound speed
        :returns: airspeed

        """
        return M*a

    def checkandsize(self, inpvar):
        """Check that input is correctly typed, then size input array. If
        scalar, leave as scalar.

        :inpvar: input array (or scalar)
        :returns: sized array (or scalar)

        """
        check_dimensioned(inpvar)
        if shape(self.altitude):  # if altitude is an array
            if shape(inpvar):  # if inpvar is an array
                if inpvar.size > self.altitude.size:
                    raise TypeError("Input airspeed array size must be less "
                                    "than or equal to the altitude array "
                                    "size.")

            sizedarr = ones(shape(self.altitude))*inpvar
        else:
            sizedarr = inpvar

        return sizedarr

    def tostring(self, full_output=False, imperial_units=False):
        """String representation of data structure.

        :full_output: set to True for full output
        :imperial_units: set to True for Imperial units
        :returns: string representation

        """
        atm_str = self.atm.tostring(full_output, imperial_units)
        atm_long_str = self.atm.tostring(full_output, imperial_units)

        M_str = f"mach = {self.mach:8.5g~P}"
        if imperial_units:
            TAS_str = f"TAS = {self.TAS.to('ft/s'):8.5g~P}"
            CAS_str = f"CAS = {self.CAS.to('ft/s'):8.5g~P}"
            EAS_str = f"EAS = {self.EAS.to('ft/s'):8.5g~P}"
            q_str = f"q_inf = {self.q_inf.to('lbf/ft^2'):8.5g~P}"
        else:  # SI units
            TAS_str = f"TAS = {self.TAS.to('m/s'):8.5g~P}"
            CAS_str = f"CAS = {self.CAS.to('m/s'):8.5g~P}"
            EAS_str = f"EAS = {self.EAS.to('m/s'):8.5g~P}"
            q_str = f"q_inf = {self.q_inf.to('Pa'):8.5g~P}"

        if full_output:
            repr_str = (f"{TAS_str}\n{CAS_str}\n{EAS_str}\n{M_str}\n{q_str}\n"
                        f"{atm_long_str}")
        else:
            repr_str = f"{TAS_str}\n{CAS_str}\n{EAS_str}\n{atm_str}"
        return repr_str

    @to_base_units_wrapper
    def M_inf_from_TAS(self, TAS):
        """Compute Mach number from true airspeed.

        :TAS: true airspeed
        :returns: Mach number

        """
        a_inf = self.atm.a
        mach = __class__.mach_number(U=TAS, a=a_inf)
        return mach

    @to_base_units_wrapper
    def TAS_from_M_inf(self, mach):
        """Compute true airspeed from Mach number.

        :mach: Mach number
        :returns: true airspeed

        """
        a_inf = self.atm.a
        TAS = __class__.mach_airspeed(M=mach, a=a_inf)
        return TAS

    @to_base_units_wrapper
    def M_inf_from_EAS(self, EAS):
        """Computer Mach number from equivalent airspeed.

        :EAS: equivalent airspeed
        :returns: Mach number

        """
        a_inf_h0 = self._atm0.a
        mach = EAS/(a_inf_h0*sqrt(self.delta))
        return mach

    @to_base_units_wrapper
    def EAS_from_M_inf(self, mach):
        """Computer Mach number from equivalent airspeed.

        :mach: Mach number
        :returns: equivalent airspeed

        """
        a_inf_h0 = self._atm0.a
        EAS = mach*(a_inf_h0*sqrt(self.delta))
        return EAS

    @to_base_units_wrapper
    def q_inf_from_TAS(self, TAS):
        """Compute dynamic pressure from true airspeed.

        :TAS: true airspeed
        :returns: dynamic pressure

        """
        rho_inf = self.atm.rho
        q_inf = 0.5*rho_inf*TAS**2
        return q_inf

    @to_base_units_wrapper
    def q_c_from_CAS(self, CAS):
        """Compute impact pressure from calibrated airspeed (accounting for
           compressibility).

        :CAS: calibrated airspeed
        :returns: impact pressure

        """
        a_h0 = self._atm0.a
        p_h0 = self._atm0.p
        # Account for compressibility with the isentropic flow equation
        M_ = __class__.mach_number(U=CAS, a=a_h0)
        q_c = __class__.isentropic_stagnation_pressure(M=M_, p=p_h0)
        return q_c

    @to_base_units_wrapper
    def CAS_from_q_c(self, q_c):
        """Compute calibrated airspeed from impact pressure (accounting for
           compressibility).

        :CAS: calibrated airspeed
        :returns: impact pressure

        """
        a_h0 = self._atm0.a
        p_h0 = self._atm0.p
        # Account for compressibility with the isentropic flow equation
        M_ = __class__.isentropic_mach(p_0=q_c, p=p_h0)
        CAS = __class__.mach_airspeed(M_, a_h0)
        return CAS

    @to_base_units_wrapper
    def M_inf_from_q_c(self, q_c):
        """Compute Mach number from impact pressure.

        :q_c: impact pressure
        :returns: Mach number

        """
        p_inf = self.atm.p
        # Isentropic flow equation
        mach = __class__.isentropic_mach(p_0=q_c, p=p_inf)

        return mach

    @to_base_units_wrapper
    def q_c_from_M_inf(self, mach):
        """Compute impact pressure from Mach number.

        :mach: Mach number
        :returns: impact pressure

        """
        p_inf = self.atm.p
        # Solve for impact pressure from isentropic flow equation:
        q_c = __class__.isentropic_stagnation_pressure(M=mach, p=p_inf)
        return q_c

    @to_base_units_wrapper
    def EAS_from_TAS(self, TAS, mach):
        """Convert airspeed to true calibrated airspeed.

        :TAS: true airspeed
        :mach: Mach number
        :returns: equivalent airspeed

        """
        a_inf_h0 = self._atm0.a
        EAS = mach*(a_inf_h0*sqrt(self.delta))
        return EAS

    @to_base_units_wrapper
    def redim_force(self, S_ref):
        """Factor to redimensionalize force from force coefficient e.g.,

            .. math::

                L & = C_L q_\\infty S_ref
                  & = C_L \\text{redim_force}

        :S_ref: reference area
        :returns: redimensionalization factor

        """
        check_area_dimensioned(S_ref)
        return self.q_inf * S_ref

    @to_base_units_wrapper
    def redim_moment(self, S_ref, ell):
        """Factor to redimensionalize moment from moment coefficient e.g.,

            .. math::

                M & = C_M q_\\infty S_ref ell
                  & = C_M \\text{redim_moment}

        :S_ref: reference area
        :ell: reference length
        :returns: redimensionalization factor

        """
        check_area_dimensioned(S_ref)
        check_length_dimensioned(ell)
        return self.q_inf * S_ref * ell

    @to_base_units_wrapper
    def reynolds_number(self, ell):
        """Compute Reynolds number for the flight condition.

        :ell: length scale
        :returns: Reynolds number

        """
        nu = self.atm.nu
        TAS = self.TAS
        Re_ell = TAS*ell/nu
        return Re_ell

    def reynolds_number_per_unit_length(self, length_unit='in'):
        """Compute Reynolds number divided by length unit.

        Reynolds number is,
            Re = U_inf * ell / nu

        For some fluid dynamics solvers the user must input the Reynolds number
        in terms of Reynolds number per-unit-length,
            Re_by_ell = Re_ell / ell_in_grid_units
                      = U_inf * (ell/ell_in_grid_units) / nu

        Assume the length scale is 5 feet.  When grid units are the same as
        standard length units, such as feet vs. feet, the term
        (ell_in_grid_units/ell) is unity, and does not change the Reynolds
        number magnitude.  For example,
            ell/ell_in_grid_units = (5 ft)/(5 ft) = 1 ft/ft

        When grid units differ from standard length units, such as inches vs.
        feet, the term (ell_in_grid_units/ell) becomes,
            ell/ell_in_grid_units = (5 ft)/(60 in) = 1/12 ft/in
        which is the conversion factor between feet and inches.

        So, the Re_by_ell term becomes
            Re_by_ell = Re_ell * (1/12 1/in)

        In this case, it is helpful to compute Reynolds number per-unit-length
        in inches (length_unit='in').


        :length_unit: desired length unit as string, ('in', 'mm', 'cm')
        :returns: Reynolds number

        """
        nu_inf = self.atm.nu
        TAS = self.TAS
        Re_by_length_unit = TAS/nu_inf
        Re_by_length_unit.ito(f"1/{length_unit}")
        return Re_by_length_unit
