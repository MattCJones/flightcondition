#!/usr/bin/env python
"""Easily convert between Mach number, true airspeed (TAS), calibrated airspeed
(CAS), and equivalent airspeed (EAS) for given altitude(s).  Additional flight
condition data and atmospheric data is computed.

Dependencies: numpy, pint, scipy

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

from numpy import ones, sqrt, shape

from ..atmosphere import Atmosphere
from ..constants import Physical as Phys
from ..units import unit, check_dimensioned, to_base_units_wrapper, \
        check_length_dimensioned, to_base_units_wrapper


class FlightCondition:

    """Easily convert between Mach number, true airspeed (TAS), calibrated
    airspeed (CAS), and equivalent airspeed (EAS) for given altitude(s).
    Additional flight condition data and atmospheric data is computed.

    All inputs must be dimensional unit quantities.

    Usage:

        from aeroutils import FlightCondition, unit, dimless

        # Compute flight conditions for a scalar or array of altitudes
        altitudes = [0, 10e3, 33.5e3] * unit('ft')
        fc = FlightCondition(altitudes, EAS=300*unit('knots'))
        print(f"Flight condition data including mach, TAS, CAS, EAS"
              f"+ atmospheric properties:\n{fc}")
        print(f"\nEven more data:\n{fc.tostring()}")

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
        fc_alt_units = FlightCondition(altitudes, EAS=154.33*unit('m/s'))
        U_TAS = fc_alt_units.TAS
        print(f"\nThe true airspeed in m/s is {U_TAS.to('m/s'):.5g}")
        print(f"The true airspeed in km/s is {U_TAS.to('km/s'):.5g}")

        # Compute additional derived quantities (see class for all options)
        print(f"\nThe dynamic pressure in psi is {fc.q_inf.to('psi'):.5g}")
        ell = 60 * unit('in')  # arbitrary length scale of interest
        print(f"The Reynolds number is {fc.reynolds_number(ell):.5g}")
        print(f"The Reynolds number per-unit-length [1/in] is "
            f"{fc.reynolds_number_per_unit_length('in'):.5g}")

    """

    def __init__(self, h_geom, mach=None, TAS=None, CAS=None, EAS=None,):
        """Constructor based on altitude and input speed in terms of Mach
        number, TAS, CAS, or EAS.  At least one input speed is required.  All
        inputs must be dimensional unit quantities.

        :h_geom: geometric altitude
        :mach: Mach number
        :TAS: true airspeed
        :CAS: calibrated airspeed
        :EAS: equivalent airspeed

        """

        h_geom = Atmosphere._process_input_altitude(h_geom)
        h_geom = h_geom[0] if h_geom.size == 1 else h_geom
        self.h_geom = h_geom

        h0 = 0 * unit('kft')
        self._atm0 = Atmosphere(h0)
        self.atm = Atmosphere(self.h_geom)

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

    def checkandsize(self, inpvar):
        """Check that input is correctly typed, then size input array. If
        scalar, leave as scalar.

        :inpvar: input array (or scalar)
        :returns: sized array (or scalar)

        """
        check_dimensioned(inpvar)
        if shape(self.h_geom):  # if h_geom is an array
            if shape(inpvar):  # if inpvar is an array
                if inpvar.size > self.h_geom.size:
                    raise TypeError("Input airspeed array size must be less "
                                    "than or equal to the altitude array "
                                    "size.")

            sizedarr = ones(shape(self.h_geom))*inpvar
        else:
            sizedarr = inpvar

        return sizedarr

    def tostring(self, short_repr=False):
        """String representation of data structure.

        :short_repr: set to True for limited output
        :returns: string representation

        """
        atm_str = self.atm.tostring(short_repr=True, imperial_units=True)
        atm_long_str = self.atm.tostring(short_repr=False, imperial_units=True)
        TAS_str = f"TAS = {self.TAS.to('knots'):8.5g~P}"
        CAS_str = f"CAS = {self.CAS.to('knots'):8.5g~P}"
        EAS_str = f"EAS = {self.EAS.to('knots'):8.5g~P}"
        M_str = f"mach = {self.mach:8.5g~P}"
        q_str = f"q_inf = {self.q_inf.to('psi'):8.5g~P}"
        if short_repr:
            repr_str = f"{TAS_str}\n{CAS_str}\n{EAS_str}\n{atm_str}"
        else:
            repr_str = (f"{TAS_str}\n{CAS_str}\n{EAS_str}\n{M_str}\n{q_str}\n"
                        f"{atm_long_str}")
        return repr_str

    def __repr__(self):
        """Output string representation of class object. Default for __str__
        :returns: string output

        """
        return self.tostring(short_repr=True)

    @to_base_units_wrapper
    def M_inf_from_TAS(self, TAS):
        """Compute Mach number from true airspeed.

        :TAS: true airspeed
        :returns: Mach number

        """
        a_inf = self.atm.a
        mach = TAS/a_inf
        return mach

    @to_base_units_wrapper
    def TAS_from_M_inf(self, mach):
        """Compute true airspeed from Mach number.

        :mach: Mach number
        :returns: true airspeed

        """
        a_inf = self.atm.a
        TAS = mach*a_inf
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
        y = Phys.gamma_air
        a_0 = self._atm0.a
        p_0 = self._atm0.p
        # Account for compressibility with the isentropic flow equation:
        # CAS ~= IAS = a_0*sqrt( (2/(y-1))*((q_c/p_0 + 1)**((y-1)/y) - 1 ))
        # Solve for impact pressure q_c:
        q_c = p_0*(-1 + (1 + ((y-1)/2)*(CAS/a_0)**2)**(y/(y-1)))
        return q_c

    @to_base_units_wrapper
    def CAS_from_q_c(self, q_c):
        """Compute calibrated airspeed from impact pressure (accounting for
           compressibility).

        :CAS: calibrated airspeed
        :returns: impact pressure

        """
        y = Phys.gamma_air
        a_0 = self._atm0.a
        p_0 = self._atm0.p
        # Account for compressibility with the isentropic flow equation:
        CAS = a_0*sqrt((2/(y-1))*((q_c/p_0 + 1)**((y-1)/y) - 1))
        return CAS

    @to_base_units_wrapper
    def M_inf_from_q_c(self, q_c):
        """Compute Mach number from impact pressure.

        :q_c: impact pressure
        :returns: Mach number

        """
        y = Phys.gamma_air
        p_inf = self.atm.p
        mach = sqrt((2/(y-1))*((q_c/p_inf + 1)**((y-1)/y) - 1))

        return mach

    @to_base_units_wrapper
    def q_c_from_M_inf(self, mach):
        """Compute impact pressure from Mach number.

        :mach: Mach number
        :returns: impact pressure

        """
        y = Phys.gamma_air
        p_inf = self.atm.p
        q_c = p_inf*(-1 + (1 + ((y-1)/2)*mach**2)**(y/(y-1)))
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
