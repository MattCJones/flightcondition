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

import numpy as np

from flightcondition.atmosphere import Atmosphere
from flightcondition.constants import PhysicalConstants as Phys
from flightcondition.common import AliasAttributes, DimensionalData,\
    _len1array_to_scalar, _property_decorators
from flightcondition.isentropicflow import IsentropicFlow
from flightcondition.nondimensional import NonDimensional
from flightcondition.units import unit, dimless, check_area_dimensioned,\
    check_dimensioned, check_length_dimensioned, check_US_length_units,\
    to_base_units_wrapper


class Velocity(DimensionalData):
    """Class to hold airspeed data. """

    varnames = {
        'TAS': 'true_airspeed',
        'CAS': 'calibrated_airspeed',
        'EAS': 'equivalent_airspeed',
        'M': 'mach_number',
        'mu_M': 'mach_angle',
        'q_inf': 'dynamic_pressure',
        'q_c': 'impact_pressure',
        'p0': 'stagnation_pressure',
        'T0': 'stagnation_temperature',
        'Tr_lamr': 'recovery_temperature_laminar',
        'Tr_turb': 'recovery_temperature_turbulent',
        'Re_by_L': 'reynolds_per_length',
    }

    def __init__(self, byalt):
        """Initialize.

        Args:
            byalt (object): Altitude object
        """
        # Link to Atmosphere data
        self._byalt = byalt
        h0 = 0 * unit('kft')
        self._atm0 = Atmosphere(h0)

    def tostring(self, full_output=True, unit_system=None, max_var_chars=0,
                 pretty_print=True):
        """String representation of data structure.

        Args:
            full_output (bool): Set to True for full output
            unit_system (str): Set to 'US' for US units or 'SI' for SI
            max_var_chars (int): Maximum number of characters in unit string
            pretty_print (bool): Pretty print output

        Returns:
            str: String representation
        """
        # Set default unit system
        if unit_system is None:
            unit_system = self._byalt.unit_system

        if unit_system == 'US':
            EAS_units   = 'knots'
            TAS_units   = 'knots'
            CAS_units   = 'knots'
            q_inf_units     = 'lbf/ft^2'
            q_c_units     = 'lbf/ft^2'
            p0_units    = 'lbf/ft^2'
            T0_units    = 'degR'
            Tr_lamr_units = 'degR'
            Tr_turb_units = 'degR'
            Re_by_L_units = '1/in'
        else:  # default to SI units
            TAS_units   = 'm/s'
            CAS_units   = 'm/s'
            EAS_units   = 'm/s'
            q_inf_units     = 'Pa'
            q_c_units     = 'Pa'
            p0_units    = 'Pa'
            T0_units    = 'degK'
            Tr_lamr_units = 'degK'
            Tr_turb_units = 'degK'
            Re_by_L_units = '1/mm'

        # Insert longer variable name into output
        max_var_chars = max([
            max([len(v) for v in self.varnames.values()]),
            max_var_chars
        ])
        M_str = self._vartostr(
            var=self.M, var_str='M', to_units='',
            max_var_chars=max_var_chars, fmt_val="10.5g",
            pretty_print=pretty_print) + "\n"
        TAS_str = self._vartostr(
            var=self.TAS, var_str='TAS', to_units=TAS_units,
            max_var_chars=max_var_chars, fmt_val="10.5g",
            pretty_print=pretty_print) + "\n"
        CAS_str = self._vartostr(
            var=self.CAS, var_str='CAS', to_units=CAS_units,
            max_var_chars=max_var_chars, fmt_val="10.5g",
            pretty_print=pretty_print) + "\n"
        EAS_str = self._vartostr(
            var=self.EAS, var_str='EAS', to_units=EAS_units,
            max_var_chars=max_var_chars, fmt_val="10.5g",
            pretty_print=pretty_print) + "\n"
        q_inf_str = self._vartostr(
            var=self.q_inf, var_str='q_inf', to_units=q_inf_units,
            max_var_chars=max_var_chars, fmt_val="10.5g",
            pretty_print=pretty_print) + "\n"
        q_c_str = self._vartostr(
            var=self.q_c, var_str='q_c', to_units=q_c_units,
            max_var_chars=max_var_chars, fmt_val="10.5g",
            pretty_print=pretty_print) + "\n"
        p0_str = self._vartostr(
            var=self.p0, var_str='p0', to_units=p0_units,
            max_var_chars=max_var_chars, fmt_val="10.5g",
            pretty_print=pretty_print) + "\n"
        T0_str = self._vartostr(
            var=self.T0, var_str='T0', to_units=T0_units,
            max_var_chars=max_var_chars, fmt_val="10.5g",
            pretty_print=pretty_print) + "\n"
        Tr_lamr_str = self._vartostr(
            var=self.Tr_lamr, var_str='Tr_lamr', to_units=Tr_lamr_units,
            max_var_chars=max_var_chars, fmt_val="10.5g",
            pretty_print=pretty_print) + "\n"
        Tr_turb_str = self._vartostr(
            var=self.Tr_turb, var_str='Tr_turb', to_units=Tr_turb_units,
            max_var_chars=max_var_chars, fmt_val="10.5g",
            pretty_print=pretty_print) + "\n"
        Re_by_L_str = self._vartostr(
            var=self.Re_by_L, var_str='Re_by_L', to_units=Re_by_L_units,
            max_var_chars=max_var_chars, fmt_val="10.4e",
            pretty_print=pretty_print)

        if np.isnan(np.atleast_1d(self.mu_M)).all():
            mu_M_str = ""
        else:
            mu_M_str = self._vartostr(
                var=self.mu_M, var_str='mu_M', to_units='deg',
                max_var_chars=max_var_chars, fmt_val="10.5g",
                pretty_print=pretty_print) + "\n"

        if full_output:
            repr_str = (f"{TAS_str}{CAS_str}{EAS_str}{M_str}{mu_M_str}"
                        f"{q_inf_str}{q_c_str}{p0_str}{T0_str}{Tr_lamr_str}"
                        f"{Tr_turb_str}{Re_by_L_str}")
        else:
            repr_str = (f"{TAS_str}{CAS_str}{EAS_str}{M_str}{Re_by_L_str}")

        return repr_str

    def _check_and_size_input(self, inpvar, default_dim=None):
        """Check that input is correctly typed, then size input array. If
        scalar, leave as scalar.

        Args:
            inpvar (object): Input array (or scalar)
            default_dim (dimension): default dimension if input is unitless

        Returns:
            object: Sized array (or scalar)
        """
        if default_dim is not None:
            inpvar *= default_dim
        check_dimensioned(inpvar)

        # Force local unit registry
        inpvar = inpvar.magnitude * unit(str(inpvar.units))

        # Require altitude and speed arrays to be equal or singular speed array
        # TODO 2022-08-12: get this to work with length
        wrongsize_msg = ("Non-singular input arrays must be equal in size.")
        if np.shape(self._byalt.h):  # if altitude is an array
            if np.shape(inpvar):  # if inpvar is an array
                if not inpvar.size == self._byalt.h.size:
                    if not inpvar.size == 1 and not self._byalt.h.size == 1:
                        raise TypeError(wrongsize_msg)

            sizedarr = np.ones(np.shape(self._byalt.h))*inpvar
        else:
            sizedarr = inpvar

        tofloat = 1.0
        sizedarr = (sizedarr * tofloat)

        return sizedarr

    @staticmethod
    def _impact_pressure(M, p):
        """Equation for impact pressure

        Args:
            M (dimless): Mach number
            p (pressure): Static pressure

        Returns:
            pressure: Impact pressure
        """
        # p0 = p + q_c
        p0 = Velocity._stagnation_pressure(M=M, p=p)
        q_c = p0 - p
        return q_c

    @to_base_units_wrapper
    def _M_from_TAS(self, TAS):
        """Compute Mach number from true airspeed.

        Args:
            TAS (speed): True airspeed

        Returns:
            dimless: Mach number
        """
        a_inf = self._byalt.a
        M = NonDimensional.mach_number(U=TAS, a=a_inf)
        return M

    @to_base_units_wrapper
    def _TAS_from_M(self, M):
        """Compute true airspeed from Mach number.

        Args:
            M (dimless): Mach number

        Returns:
            speed: True airspeed
        """
        a_inf = self._byalt.a
        TAS = NonDimensional.mach_airspeed(M, a_inf)
        return TAS

    @to_base_units_wrapper
    def _M_from_EAS(self, EAS):
        """Computer Mach number from equivalent airspeed.

        Args:
            EAS (speed): Equivalent airspeed

        Returns:
            dimless: Mach number
        """
        a_inf_h0 = self._atm0.a
        p_inf_h0 = self._atm0.p
        p_inf = self._byalt.p
        delta = p_inf/p_inf_h0
        M = EAS/(a_inf_h0*np.sqrt(delta))
        return M

    @to_base_units_wrapper
    def _EAS_from_M(self, M):
        """Computer Mach number from equivalent airspeed.

        Args:
            M (dimless): Mach number

        Returns:
            speed: Equivalent airspeed
        """
        a_inf_h0 = self._atm0.a
        p_inf_h0 = self._atm0.p
        p_inf = self._byalt.p
        delta = p_inf/p_inf_h0
        EAS = M*(a_inf_h0*np.sqrt(delta))
        return EAS

    @to_base_units_wrapper
    def _q_c_from_CAS(self, CAS):
        """Compute impact pressure from calibrated airspeed (accounting for
           compressibility).

        Args:
            CAS (speed): Calibrated airspeed

        Returns:
            pressure: Impact pressure
        """
        a_h0 = self._atm0.a
        p_h0 = self._atm0.p
        # Account for compressibility with the isentropic flow equation
        M_ = NonDimensional.mach_number(U=CAS, a=a_h0)
        q_c = __class__._impact_pressure(M=M_, p=p_h0)
        return q_c

    @to_base_units_wrapper
    def _CAS_from_q_c(self, q_c):
        """Compute calibrated airspeed from impact pressure (accounting for
           compressibility).

        Args:
            CAS (speed): Calibrated airspeed

        Returns:
            pressure: Impact pressure

        """
        a_h0 = self._atm0.a
        p_h0 = self._atm0.p
        # Account for compressibility with the isentropic flow equation
        # M_should = __class__._isentropic_mach(p0=q_c, p=p_h0)  # DEPRECATED
        M = IsentropicFlow.M_from_p0_by_p((q_c+p_h0)/p_h0)
        CAS = NonDimensional.mach_airspeed(M, a_h0)  # subsonic
        return CAS

    @to_base_units_wrapper
    def _M_from_q_c(self, q_c):
        """Compute Mach number from impact pressure.

        Args:
            q_c (pressure): Impact pressure

        Returns:
            dimless: Mach number
        """
        p_inf = self._byalt.p
        # Isentropic flow equation
        # M_should = __class__._isentropic_mach(p0=q_c, p=p_inf)  # DEPRECATED
        M = IsentropicFlow.M_from_p0_by_p((q_c+p_inf)/p_inf)

        return M

    @to_base_units_wrapper
    def _q_c_from_M(self, M):
        """Compute impact pressure from Mach number.

        Args:
            M: Mach number

        Returns:
            impact pressure
        """
        p_inf = self._byalt.p
        # Solve for impact pressure from isentropic flow equation:
        q_c = __class__._impact_pressure(M=M, p=p_inf)
        return q_c

    @to_base_units_wrapper
    def _EAS_from_TAS(self, TAS, M):
        """Convert airspeed to true calibrated airspeed.

        Args:
            TAS (speed): True airspeed
            M (dimless): Mach number

        Returns:
            speed: Equivalent airspeed
        """
        a_inf_h0 = self._atm0.a
        p_inf_h0 = self._atm0.p
        p_inf = self._byalt.p
        delta = p_inf/p_inf_h0
        EAS = M*(a_inf_h0*np.sqrt(delta))
        return EAS

    @staticmethod
    @_len1array_to_scalar
    def _mach_angle(M):
        """Compute Mach angle

        Args:
            M (dimless): Mach number

        Returns:
            angle: mach angle
        """
        # Create new array and split computation to subsonic and supersonic
        M = np.atleast_1d(M)
        mu_M = M*0.0

        sub = M < 1  # subsonic filter
        if sub.any():
            mu_M[sub] = np.nan

        sup = M >= 1  # supersonic filter
        mu_M[sup] = np.arcsin(1.0/M[sup])

        mu_M.ito('deg')
        return mu_M

    @staticmethod
    @to_base_units_wrapper
    def _q_inf_from_TAS(TAS, rho):
        """Compute dynamic pressure from true airspeed.

        Args:
            TAS (speed): True airspeed
            rho (density): Ambient density

        Returns:
            pressure: Dynamic pressure
        """
        q_inf = 0.5*rho*TAS**2
        return q_inf

    @staticmethod
    def _stagnation_pressure(M, p, y=Phys.gamma_air):
        """Equation for stagnation pressure, i.e. bringing flow to rest
        isentropically.

        Assumes:
            Adiabatic
            Reversible

        Args:
            M (dimless): Mach number
            p (pressure): Static pressure
            y (dimless): ratio of specific heats

        Returns:
            pressure: Stagnation pressure
        """
        p0 = IsentropicFlow.p0_by_p(M, y)*p
        return p0

    @staticmethod
    def _stagnation_temperature(M, T, y=Phys.gamma_air):
        """Adiabatic flow equation for stagnation temperature, i.e. temperature
        brought to rest isentropically (note that the equation does not assume
        isentropic flow).

        Assumes:
            Adiabatic

        Args:
            M (pressure): Mach number
            T (temperature): Static (ambient) temperature
            y (dimless): ratio of specific heats

        Returns:
            temperature: Stagnation temperature
        """

        T0 = IsentropicFlow.T0_by_T(M, y)*T
        return T0

    @staticmethod
    def _recovery_temperature_laminar(M, T, y=Phys.gamma_air, Pr=Phys.Pr_air):
        """Adiabiatic wall temperature on infinite flat plate in laminar flow.

        Assumes:
            Adiabatic

        Args:
            M (pressure): Mach number
            T (temperature): Static (ambient) temperature
            y (dimless): ratio of specific heats
            Pr (dimless): Prandtl number

        Returns:
            temperature: Stagnation temperature
        """
        r = Pr**(1/2)
        T0 = IsentropicFlow.T0_by_T(M, y, r)*T
        return T0

    @staticmethod
    def _recovery_temperature_turbulent(M, T, y=Phys.gamma_air,
                                        Pr=Phys.Pr_air):
        """Adiabiatic wall temperature on infinite flat plate in turbulent
        flow.

        Assumes:
            Adiabatic

        Args:
            M (pressure): Mach number
            T (temperature): Static (ambient) temperature
            y (dimless): ratio of specific heats
            Pr (dimless): Prandtl number

        Returns:
            temperature: Stagnation temperature
        """
        r = Pr**(1/3)
        T0 = IsentropicFlow.T0_by_T(M, y, r)*T
        return T0

    @_property_decorators
    def M(self):
        """Get Mach number :math:`M` """
        return self._M

    @M.setter
    def M(self, M):
        """Mach number :math:`M` """
        self._M = self._check_and_size_input(M, default_dim=dimless)
        self._TAS = self._TAS_from_M(self.M)
        self._EAS = self._EAS_from_TAS(self.TAS, self.M)
        self._q_c = self._q_c_from_M(self.M)
        self._CAS = self._CAS_from_q_c(self.q_c)

    @_property_decorators
    def TAS(self):
        """Get true airspeed. """
        return self._TAS

    @TAS.setter
    def TAS(self, TAS):
        """Set true airspeed. """
        self._TAS = self._check_and_size_input(TAS)
        self._M = self._M_from_TAS(TAS)
        self._EAS = self._EAS_from_TAS(self.TAS, self.M)
        self._q_c = self._q_c_from_M(self.M)
        self._CAS = self._CAS_from_q_c(self.q_c)

    @_property_decorators
    def CAS(self):
        """Get calibrated airspeed. """
        return self._CAS

    @CAS.setter
    def CAS(self, CAS):
        """Calibrated airspeed. """
        self._CAS = self._check_and_size_input(CAS)
        self._q_c = self._q_c_from_CAS(self.CAS)
        self._M = self._M_from_q_c(self.q_c)
        self._TAS = self._TAS_from_M(self.M)
        self._EAS = self._EAS_from_TAS(self.TAS, self.M)

    @_property_decorators
    def EAS(self):
        """Get equivalent airspeed. """
        return self._EAS

    @EAS.setter
    def EAS(self, EAS):
        """Set equivalent airspeed. """
        self._EAS = self._check_and_size_input(EAS)
        self._M = self._M_from_EAS(self.EAS)
        self._TAS = self._TAS_from_M(self.M)
        self._q_c = self._q_c_from_M(self.M)
        self._CAS = self._CAS_from_q_c(self.q_c)

    @_property_decorators
    def mu_M(self):
        """Get Mach angle :math:`\\mu_M` """
        return __class__._mach_angle(self.M)

    @_property_decorators
    def q_c(self):
        """Get impact pressure :math:`q_c`"""
        return self._q_c

    @_property_decorators
    def q_inf(self):
        """Get dynamic pressure :math:`q_\\infty`"""
        q_inf = __class__._q_inf_from_TAS(TAS=self.TAS, rho=self._byalt.rho)
        return q_inf

    @_property_decorators
    def p0(self):
        """Get stagnation pressure :math:`p_0`"""
        p0 = self._q_inf_from_TAS(TAS=self.TAS, rho=self._byalt.rho)
        p0 = __class__._stagnation_pressure(M=self.M, p=self._byalt.p)
        return p0

    @_property_decorators
    def T0(self):
        """Get stagnation temperature :math:`T_0`"""
        T0 = __class__._stagnation_temperature(M=self.M, T=self._byalt.T)
        return T0

    @_property_decorators
    def Tr_lamr(self):
        """Get recovery temperature (laminar) :math:`Tr_{laminar}`"""
        Tr_lamr = __class__._recovery_temperature_laminar(M=self.M,
                                                          T=self._byalt.T)
        return Tr_lamr

    @_property_decorators
    def Tr_turb(self):
        """Get recovery temperature (turbulent) :math:`Tr_{turbulent}`"""
        Tr_turb = __class__._recovery_temperature_turbulent(M=self.M,
                                                            T=self._byalt.T)
        return Tr_turb

    @_property_decorators
    def Re_by_L(self):
        """Get Reynolds number per unit length :math:`Re_L`"""
        length_unit = 'in' if self._byalt.unit_system == 'US' else None
        Re_by_L = NonDimensional.reynolds_per_length(U=self.TAS,
                                                     nu=self._byalt.nu,
                                                     length_unit=length_unit)
        return Re_by_L


class Length(DimensionalData):
    """Class to hold length data. """

    varnames = {
        'L': 'length_scale',
        'Re': 'reynolds_number',
    }

    def __init__(self, byalt, byvel):
        """Initialize.

        Args:
            byalt (object): Altitude object
            byvel (object): Velocity object
        """
        # Link to Atmosphere and Velocity data
        self._byalt = byalt
        self._byvel = byvel

    def tostring(self, full_output=True, unit_system=None, max_var_chars=0,
                 pretty_print=True):
        """String representation of data structure.

        Args:
            full_output (bool): Set to True for full output
            unit_system (str): Set to 'US' for US units or 'SI' for SI
            max_var_chars (int): Maximum number of characters in unit string
            pretty_print (bool): Pretty print output

        Returns:
            str: String representation
        """
        # Set default unit system
        if unit_system is None:
            unit_system = self._byalt.unit_system

        if unit_system == 'US':
            L_units = 'ft'
        else:  # default to SI units
            L_units = 'm'

        # Insert longer variable name into output
        max_var_chars = max([
            max([len(v) for v in self.varnames.values()]),
            max_var_chars
        ])
        L_str = self._vartostr(
            var=self.L, var_str='L', to_units=L_units,
            max_var_chars=max_var_chars, fmt_val="10.5g",
            pretty_print=pretty_print)
        Re_str = self._vartostr(
            var=self.Re, var_str='Re', to_units='',
            max_var_chars=max_var_chars, fmt_val="10.4e",
            pretty_print=pretty_print)

        if full_output:
            repr_str = (f"{L_str}\n{Re_str}")
        else:
            repr_str = (f"{L_str}\n{Re_str}")

        return repr_str

    @_property_decorators
    def L(self):
        """Get length scale :math:`L`"""
        return self._L

    @L.setter
    def L(self, L):
        """Set length scale :math:`L`"""
        # Verify that length input is dimensional length quantity
        L = self._byvel._check_and_size_input(L)
        check_length_dimensioned(L)
        self._L = L

    @_property_decorators
    def Re(self):
        """Get Reynolds number :math:`Re`"""
        Re = NonDimensional.reynolds_number(
            U=self._byvel.TAS, L=self.L, nu=self._byalt.nu)
        return Re

    @_property_decorators
    def Kn(self):
        """Get Knudsen number :math:`K_n`"""
        Kn = NonDimensional.knudsen_number(self.L, self._byalt.MFP)
        return Kn


class FlightCondition(DimensionalData):
    """Easily convert between Mach number, true airspeed (TAS), calibrated
    airspeed (CAS), and equivalent airspeed (EAS) for given altitude(s).
    Additional flight condition data and atmospheric data is computed.

    All inputs must be dimensional unit quantities.

    Usage:
        from flightcondition import FlightCondition, unit, dimless

        # Compute flight condition at 3 km, Mach 0.5
        fc = FlightCondition(3*unit('km'), M=0.5)

        # Uncomment to print summary of flight condition quantities:
        #print(f"{fc}")

        # Uncomment to print abbreviated output in US units:
        #print(f"\n{fc.tostring(full_output=False, unit_system='US')}")

        # Access true, calibrated, equivalent airspeeds
        KTAS = fc.byvel.TAS.to('knots')
        KCAS = fc.byvel.CAS.to('knots')
        KEAS = fc.byvel.EAS.to('knots')
        print(f"Flying at {KTAS.magnitude:.4g} KTAS,"
            f" which is {KCAS.magnitude:.4g} KCAS,"
            f" or {KEAS.magnitude:.4g} KEAS")
        # >>> Flying at 319.4 KTAS, which is 277.7 KCAS, or 275.1 KEAS

        # Access atmospheric data (see Atmosphere class for more)
        atm = fc.byalt  # access Atmosphere object
        h, p, T, rho, nu, a = atm.h, atm.p, atm.T, atm.rho, atm.nu, atm.a
        print(f"The ambient temperature at {h.to('km'):.4g} is {T:.4g}")
        # >>> The ambient temperature at 3 km is 268.7 K

        # Compute again instead using true airspeed and altitude in km
        fc = FlightCondition(3.048*unit('km'), TAS=401.7*unit('mph'))
        #print(f"{fc}")  # uncomment to print output

        # Compute for a range of altitudes at 275.14 knots-equivalent
        # airspeed with a characteristic length scale of 10 meters
        fc = FlightCondition([0, 9.8425, 20]*unit('kft'),
                            EAS=275.14*unit('kt'),
                            L=10*unit('m'))

        # Compute additional derived quantities
        # Explore the class data structure for all options
        print(f"\nThe dynamic pressure in psi is "
            f"{fc.byvel.q_inf.to('psi'):.3g}")
        # >>> The dynamic pressure in psi is [1.78 1.78 1.78] psi
        print(f"The Reynolds number is {fc.bylen.Re:.3g}")
        # >>> The Reynolds number is [9.69e+07 8.82e+07 7.95e+07]

        # Alternatively access quantities by their full name
        print(fc.byvel.TAS == fc.byname.true_airspeed)
        # >>> [ True  True  True]
    """

    def __init__(
        self, h=None, M=None, TAS=None, CAS=None, EAS=None, L=None,
        unit_system="", **kwargs,
    ):
        """Constructor based on altitude and input airspeed in terms of Mach
        number, TAS, CAS, or EAS.  Input at least one format of airspeed at the
        desired altitude.  All inputs must be dimensional unit quantities.

        Args:
            h (length): Geometric altitude - aliases are 'alt', 'altitude'
            M (dimless): Mach number - aliases are 'mach', 'Mach', 'M_inf',
                'mach_number'
            TAS (speed): True airspeed - aliases are 'tas', 'true_airspeed',
                'U_inf', 'V_inf'
            CAS (speed): Calibrated airspeed - aliases are 'cas',
                'calibrated_airspeed'
            EAS (speed): Equivalent airspeed - aliases are 'eas',
                'equivalent_airspeed'
            L (length): Length scale - aliases are 'ell', 'bylen', 'length',
                'length_scale', 'l'
            unit_system (str): Set to 'US' for US units or 'SI' for SI
        """

        # Preprocess needed altitude-based quantities
        # Automatically process altitude through Atmosphere class
        self.byalt = Atmosphere(h=h, unit_system=unit_system, **kwargs)

        # Check for hidden aliases
        M_aliases = ['mach', 'Mach', 'M_inf', 'mach_number']
        if M is None:
            M = __class__._arg_from_alias(M_aliases, kwargs)
        TAS_aliases = ['tas', 'true_airspeed', 'U_inf', 'V_inf', 'VTAS',
                       'vtas']
        if TAS is None:
            TAS = __class__._arg_from_alias(TAS_aliases, kwargs)
        CAS_aliases = ['cas', 'calibrated_airspeed', 'VCAS', 'vcas']
        if CAS is None:
            CAS = __class__._arg_from_alias(CAS_aliases, kwargs)
        EAS_aliases = ['eas', 'equivalent_airspeed', 'VEAS', 'veas']
        if EAS is None:
            EAS = __class__._arg_from_alias(EAS_aliases, kwargs)
        L_aliases = ['ell', 'len', 'length', 'length_scale', 'l']
        if L is None:
            L = __class__._arg_from_alias(L_aliases, kwargs)

        # Check if KTAS, KCAS, or KEAS input and append knots unit if so
        KTAS_aliases = ['KTAS', 'ktas', 'knots_true_airspeed']
        if TAS is None:
            KTAS = __class__._arg_from_alias(KTAS_aliases, kwargs)
            TAS = None if KTAS is None else KTAS * unit('knots')
        KCAS_aliases = ['KCAS', 'kcas', 'knots_calibrated_airspeed']
        if CAS is None:
            KCAS = __class__._arg_from_alias(KCAS_aliases, kwargs)
            CAS = None if KCAS is None else KCAS * unit('knots')
        KEAS_aliases = ['KEAS', 'keas', 'knots_equivalent_airspeed']
        if EAS is None:
            KEAS = __class__._arg_from_alias(KEAS_aliases, kwargs)
            EAS = None if KEAS is None else KEAS * unit('knots')

        # Compute airspeed-based quantities
        self.byvel = Velocity(self.byalt)

        # Use Mach=0 if no airspeed is input
        if M is None and TAS is None and CAS is None and EAS is None:
            M = 0*dimless

        if M is not None:
            self.byvel.M = M
        elif TAS is not None:
            self.byvel.TAS = TAS
        elif CAS is not None:
            self.byvel.CAS = CAS
        elif EAS is not None:
            self.byvel.EAS = EAS
        else:
            raise TypeError("Input M, TAS, CAS, or EAS")

        # Check that computations are within valid Mach number limits
        M = np.atleast_1d(self.byvel.M)
        self._mach_min = 0 * dimless
        self._mach_max = 30 * dimless
        if (M < self._mach_min).any() or (self._mach_max < M).any():
            raise ValueError(
                f"Mach number is out of bounds "
                f"({self._mach_min:.5g} < M < {self._mach_max:.5g})"
            )

        # Compute length-scale-based quantities
        self.bylen = Length(self.byalt, self.byvel)

        # If length scale is not input, default to unity with dimentionals unit
        # based on US or SI determination
        if L is None:
            L_unit = unit('ft') if self.unit_system == 'US' else unit('m')
            self.bylen.L = 1.0 * L_unit
        else:  # length scale is input by user
            self.bylen.L = L

            # If altitude is 0, but length scale is input, determine if US
            # units based on the input length scale
            if self.byalt.h.magnitude.all() == 0:
                if check_US_length_units(L):
                    self.unit_system = 'US'

        # Initialize access by full quantity name through .byname.<name>
        self.byvel.byname = AliasAttributes(
            varsobj_arr=[self.byvel, ], varnames_dict_arr=[Velocity.varnames, ]
        )
        self.bylen.byname = AliasAttributes(
            varsobj_arr=[self.bylen, ], varnames_dict_arr=[Length.varnames, ]
        )

        # Set references at base level too
        varsobj_arr = [self.byalt, self.byvel, self.bylen, ]
        varnames_dict_arr = [Atmosphere.varnames,
                             Velocity.varnames,
                             Length.varnames, ]
        own_dict_arr = [
            {key: key for key in dict_.keys()} for dict_ in varnames_dict_arr
        ]

        self.byvar = AliasAttributes(
            varsobj_arr=varsobj_arr,
            varnames_dict_arr=own_dict_arr
        )

        self.varnames = {}
        self.varnames.update(Atmosphere.varnames)
        self.varnames.update(Velocity.varnames)
        self.varnames.update(Length.varnames)

        self.byname = AliasAttributes(
            varsobj_arr=varsobj_arr,
            varnames_dict_arr=varnames_dict_arr
        )

    def tostring(self, full_output=True, unit_system=None, pretty_print=True):
        """String representation of data structure.

        Args:
            full_output (bool): Set to True for full output
            unit_system (str): Set to 'US' for US units or 'SI' for SI
            pretty_print (bool): Pretty print output

        Returns:
            str: String representation
        """
        if unit_system is not None:
            self.unit_system = unit_system

        max_var_chars = max([
            max([len(v) for v in Atmosphere.varnames.values()]),
            max([len(v) for v in Velocity.varnames.values()]),
            max([len(v) for v in Length.varnames.values()]),
        ])
        alti_str = self.byalt.tostring(full_output, self.unit_system,
                                       pretty_print=pretty_print,
                                       max_var_chars=max_var_chars)
        spee_str = self.byvel.tostring(full_output, self.unit_system,
                                       max_var_chars=max_var_chars,
                                       pretty_print=pretty_print)
        leng_str = self.bylen.tostring(full_output, self.unit_system,
                                       max_var_chars=max_var_chars,
                                       pretty_print=pretty_print)

        unit_str = self.unit_system
        ext_str = "full_output=True" if full_output else "full_output=False"
        top_hdr = f"   Flight Condition (unit_system={unit_str}, {ext_str})"
        lin_str = "==========================================================="
        alt_hdr = "------------------  Altitude Quantities  ------------------"
        vel_hdr = "------------------  Velocity Quantities  ------------------"
        len_hdr = "------------------   Length Quantities   ------------------"

        repr_str = (f"{lin_str}\n{top_hdr}\n{lin_str}"
                    f"\n{alt_hdr}" f"\n{alti_str}"
                    f"\n{vel_hdr}" f"\n{spee_str}"
                    f"\n{len_hdr}" f"\n{leng_str}"
                    )

        return repr_str

    @property
    def unit_system(self):
        """Unit system to use: 'SI', 'US', etc.  Available unit systems given
        by dir(unit.sys).

        Returns:
            str: Unit system
        """
        return self.byalt.unit_system

    @unit_system.setter
    def unit_system(self, unit_system):
        """Unit system to use: 'SI', 'US', etc.  Available unit systems given
        by dir(unit.sys).

        Args:
            unit_system (str): Unit system
        """
        self.byalt.unit_system = unit_system

    def redim_force(self, S_ref):
        """Factor to redimensionalize force from force coefficient e.g.,

            .. math::

                L & = C_L q_\\infty S_ref
                  & = C_L \\text{redim_force}

        Args:
            S_ref (area): Reference area

        Returns:
            force: Redimensionalization factor
        """
        check_area_dimensioned(S_ref)
        redim_force = self.byvel.q_inf * S_ref

        # Set to force unit since to_base_units() gives mass*length/time^2
        if self.unit_system == 'US':
            redim_force.ito('lbf')
        elif self.unit_system == 'SI':
            redim_force.ito('N')

        return redim_force

    def redim_moment(self, S_ref, L):
        """Factor to redimensionalize moment from moment coefficient e.g.,

            .. math::

                M & = C_M q_\\infty S_ref L
                  & = C_M \\text{redim_moment}

        Args:
            S_ref (area): Reference area
            L (length): Reference length

        Returns:
            moment: Redimensionalization factor
        """
        check_area_dimensioned(S_ref)
        check_length_dimensioned(L)
        redim_moment = self.byvel.q_inf * S_ref * L

        # Set to force unit since to_base_units() gives mass*length/time^2
        if self.unit_system == 'US':
            redim_moment.ito('ft lbf')
        elif self.unit_system == 'SI':
            redim_moment.ito('m N')

        return redim_moment
