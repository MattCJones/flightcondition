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

import warnings

import numpy as np

from flightcondition.atmosphere import Atmosphere
from flightcondition.boundarylayer import BoundaryLayer
from flightcondition.constants import PhysicalConstants as Phys
from flightcondition.common import AliasAttributes, DimensionalData,\
    _len1array_to_scalar, _property_decorators
from flightcondition.isentropicflow import IsentropicFlow
from flightcondition.nondimensional import NonDimensional
from flightcondition.units import unit, dimless, check_area_dimensioned,\
    check_dimensioned, check_dimensionless, check_length_dimensioned,\
    check_US_length_units, to_base_units_wrapper


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

    def tostring(self, full_output=None, units=None, max_var_chars=0,
                 pretty_print=True):
        """String representation of data structure.

        Args:
            full_output (bool): Set to True for full output
            units (str): Set to 'US' for US units or 'SI' for SI
            max_var_chars (int): Maximum number of characters in unit string
            pretty_print (bool): Pretty print output

        Returns:
            str: String representation
        """
        # Set default unit system
        if units is None:
            units = self._byalt.units

        if units == 'US':
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

        # Assemble output string
        if full_output:
            repr_str = (f"{TAS_str}{CAS_str}{EAS_str}{M_str}{mu_M_str}"
                        f"{q_inf_str}{q_c_str}{p0_str}{T0_str}{Tr_lamr_str}"
                        f"{Tr_turb_str}{Re_by_L_str}")
        else:
            repr_str = (f"{TAS_str}{CAS_str}{EAS_str}{M_str}{Re_by_L_str}")

        return repr_str

    @staticmethod
    def _compare_input_size(input_a, input_b):
        """Enforce non-singular input arrays to be equal in size

        Args:
            input_a (object): Input A
            input_b (object): Input B

        Returns:
            object: sized_arr
        """
        wrongsize_msg = ("Non-singular input arrays must be equal in size.")

        if np.shape(input_b):  # if array
            if np.shape(input_a):  # if array
                if not input_a.size == 1 and not input_b.size == 1:
                    if not input_a.size == input_b.size:
                        raise TypeError(wrongsize_msg)
            sizedarr = np.ones(np.shape(input_b))*input_a
        else:
            sizedarr = input_a
        return sizedarr

    def _check_and_size_input(self, input_var, input_alt=None, input_vel=None):
        """Check that input is correctly typed, then size input array. If
        scalar, leave as scalar.

        Args:
            input_var (object): Input array (or scalar)
            input_alt (object): Input array of altitude level
            input_vel (object): Input array of velocity level (when assessing
                length scale as input_var)

        Returns:
            object: Sized array (or scalar)
        """
        check_dimensioned(input_var)

        # Force local unit registry
        input_var = input_var.magnitude * unit(str(input_var.units))

        # Require input, base, and secondary arrays to be singular or
        # non-singular but equal in size
        sizedarr = __class__._compare_input_size(input_var, input_alt)
        if input_vel is not None:
            sizedarr = __class__._compare_input_size(input_var, input_vel)

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
        y = Phys.gamma_air
        p0 = IsentropicFlow.p0_by_p(M, y)*p
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
        TAS = NonDimensional.mach_velocity(M, a_inf)
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
        CAS = NonDimensional.mach_velocity(M, a_h0)  # subsonic
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
    @unit.wraps(unit.rad, dimless.units)
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

    @_property_decorators
    def M(self):
        """Get Mach number :math:`M` """
        return self._M

    @M.setter
    def M(self, M):
        """Mach number :math:`M` """
        M *= dimless  # add dimless for raw float input
        self._M = self._check_and_size_input(M, input_alt=self._byalt.h)
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
        self._TAS = self._check_and_size_input(TAS, input_alt=self._byalt.h)
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
        self._CAS = self._check_and_size_input(CAS, input_alt=self._byalt.h)
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
        self._EAS = self._check_and_size_input(EAS, input_alt=self._byalt.h)
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
        M = self.M
        p = self._byalt.p
        y = Phys.gamma_air
        p0 = IsentropicFlow.p0_by_p(M, y)*p
        return p0

    @_property_decorators
    def T0(self):
        """Get stagnation temperature :math:`T_0`"""
        M = self.M
        T = self._byalt.T
        y = Phys.gamma_air
        T0 = IsentropicFlow.T0_by_T(M, y)*T
        return T0

    @_property_decorators
    def Tr_lamr(self):
        """Get recovery temperature (laminar) :math:`Tr_{laminar}`"""
        Tr_lamr = BoundaryLayer.recovery_temperature_laminar(
            M=self.M, T=self._byalt.T)
        return Tr_lamr

    @_property_decorators
    def Tr_turb(self):
        """Get recovery temperature (turbulent) :math:`Tr_{turbulent}`"""
        Tr_turb = BoundaryLayer.recovery_temperature_turbulent(
            M=self.M, T=self._byalt.T)
        return Tr_turb

    @_property_decorators
    def Re_by_L(self):
        """Get Reynolds number per unit length :math:`Re_L`"""
        length_unit = 'in' if self._byalt.units == 'US' else None
        Re_by_L = NonDimensional.reynolds_per_length(U=self.TAS,
                                                     nu=self._byalt.nu,
                                                     length_unit=length_unit)
        return Re_by_L


class Length(DimensionalData):
    """Class to hold length data. """

    varnames = {
        'L': 'length_scale',
        'Re': 'reynolds_number',
        # 'Kn': 'knudsen_number',
        'h_BL_lamr': 'boundary_thickness_laminar',
        'h_BL_turb': 'boundary_thickness_turbulent',
        'Cf_lamr': 'friction_coefficient_laminar',
        'Cf_turb': 'friction_coefficient_turbulent',
        'h_yplus1': 'wall_distance_yplus1',
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

    def tostring(self, full_output=True, units=None, max_var_chars=0,
                 pretty_print=True):
        """String representation of data structure.

        Args:
            full_output (bool): Set to True for full output
            units (str): Set to 'US' for US units or 'SI' for SI
            max_var_chars (int): Maximum number of characters in unit string
            pretty_print (bool): Pretty print output

        Returns:
            str: String representation
        """
        # Set default unit system
        if units is None:
            units = self._byalt.units

        if units == 'US':
            L_units = 'ft'
            h_BL_units = 'in'
        else:  # default to SI units
            L_units = 'm'
            h_BL_units = 'mm'

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
            max_var_chars=max_var_chars, fmt_val="10.5g",
            pretty_print=pretty_print)
        h_BL_lamr_str = self._vartostr(
            var=self.h_BL_lamr, var_str='h_BL_lamr', to_units=h_BL_units,
            max_var_chars=max_var_chars, fmt_val="10.5g",
            pretty_print=pretty_print)
        h_BL_turb_str = self._vartostr(
            var=self.h_BL_turb, var_str='h_BL_turb', to_units=h_BL_units,
            max_var_chars=max_var_chars, fmt_val="10.5g",
            pretty_print=pretty_print)
        Cf_lamr_str = self._vartostr(
            var=self.Cf_lamr, var_str='Cf_lamr', to_units=dimless,
            max_var_chars=max_var_chars, fmt_val="10.5g",
            pretty_print=pretty_print)
        Cf_turb_str = self._vartostr(
            var=self.Cf_turb, var_str='Cf_turb', to_units=dimless,
            max_var_chars=max_var_chars, fmt_val="10.5g",
            pretty_print=pretty_print)
        h_yplus1_str = self._vartostr(
            var=self.h_yplus1, var_str='h_yplus1', to_units=h_BL_units,
            max_var_chars=max_var_chars, fmt_val="10.5g",
            pretty_print=pretty_print)

        if full_output:
            repr_str = (f"{L_str}\n{Re_str}\n{h_BL_lamr_str}\n"
                        f"{h_BL_turb_str}\n{Cf_lamr_str}\n{Cf_turb_str}\n"
                        f"{h_yplus1_str}")
        else:
            repr_str = (f"{L_str}\n{Re_str}")

        return repr_str

    @to_base_units_wrapper
    def wall_distance_from_yplus(self, yplus):
        """Compute wall distance from yplus value.  Assume turbulent flow over
        a flat plate.

        :math:`y^+ = \\frac{y u_\\tau}{\\nu}`

        where

        :math:`u_\\tau = \\sqrt{\\frac{\\tau_w}{\\rho}}

        Args:
            yplus (dimless): yplus value

        Returns:
            length: Distance from wall
        """
        yplus *= dimless
        check_dimensionless(yplus)

        rho = self._byalt.rho
        nu = self._byalt.nu
        q_inf = self._byvel.q_inf
        with warnings.catch_warnings():  # catch warning for divide by 0
            warnings.simplefilter("ignore")
            Cf_turb = self.Cf_turb

            tau_wall = Cf_turb*q_inf
            u_tau = np.sqrt(tau_wall/rho)
            y = yplus*nu/u_tau

        return y

    @_property_decorators
    def L(self):
        """Get length scale :math:`L`"""
        return self._L

    @L.setter
    def L(self, L):
        """Set length scale :math:`L`"""
        # Verify that length input is dimensional length quantity
        L = self._byvel._check_and_size_input(L, input_alt=self._byalt.h,
                                              input_vel=self._byvel.M)
        check_length_dimensioned(L)
        self._L = L

    @_property_decorators
    def Re(self):
        """Get Reynolds number :math:`Re`"""
        Re = NonDimensional.reynolds_number(
            U=self._byvel.TAS, L=self.L, nu=self._byalt.nu)
        return Re

    @Re.setter
    def Re(self, Re):
        """Set Reynolds number :math:`Re`
        Set the true airspeed based on Reynolds number and length scale. """
        Re *= dimless  # add dimless for raw float input
        Re = self._byvel._check_and_size_input(
            Re, input_alt=self._byalt.h, input_vel=self._byvel.M)
        self._byvel.TAS = NonDimensional.reynolds_number_velocity(
            Re_L=Re, L=self.L, nu=self._byalt.nu)

    # @_property_decorators  # disable for now
    # def Kn(self):
    #     """Get Knudsen number :math:`K_n`"""
    #     Kn = NonDimensional.knudsen_number(self.L, self._byalt.MFP)
    #     return Kn

    @_property_decorators
    def h_BL_lamr(self):
        """Get laminar Boundary Layer Thickness :math:`\\h_BL_{lamr}` """
        x = self.L
        Re_x = NonDimensional.reynolds_number(
            U=self._byvel.TAS, L=x, nu=self._byalt.nu)
        h_BL_lamr = BoundaryLayer.flat_plate_boundary_layer_lamr(
            x=self.L, Re_x=Re_x)
        return h_BL_lamr

    @_property_decorators
    def h_BL_turb(self):
        """Get turbulent Boundary Layer Thickness :math:`\\h_BL_{turb}` """
        x = self.L
        Re_x = NonDimensional.reynolds_number(
            U=self._byvel.TAS, L=x, nu=self._byalt.nu)
        h_BL_turb = BoundaryLayer.flat_plate_boundary_layer_turb(
            x=self.L, Re_x=Re_x)
        return h_BL_turb

    @_property_decorators
    def Cf_lamr(self):
        """Get laminar skin friction coefficient :math:`\\Cf_{lamr}` """
        x = self.L
        Re_x = NonDimensional.reynolds_number(
            U=self._byvel.TAS, L=x, nu=self._byalt.nu)
        Cf_lamr = BoundaryLayer.flat_plate_skin_friction_coeff_lamr(Re_x=Re_x)
        return Cf_lamr

    @_property_decorators
    def Cf_turb(self):
        """Get turbulent skin friction coefficient :math:`\\Cf_{turb}` """
        Cf_turb = BoundaryLayer.flat_plate_skin_friction_coeff_turb(
            Re_x=self.Re, M=self._byvel.M)
        return Cf_turb

    @_property_decorators
    def h_yplus1(self):
        """Get height from flat plate wall in turbulent flow where yplus=1. """
        h_yplus1 = self.wall_distance_from_yplus(1)
        return h_yplus1


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
        #print(f"\n{fc.tostring(full_output=False, units='US')}")

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
        self, h=None, M=None, TAS=None, CAS=None, EAS=None, L=None, Re=None,
        units="", full_output=None, **kwargs,
    ):
        """Constructor based on altitude and input velocity in terms of Mach
        number, TAS, CAS, or EAS.  Input altitude, one format of velocity, and
        length scale.  Reynolds number can be input as an alternative to either
        velocity or length scale but not both.  All inputs must be dimensional
        unit quantities.

        Args:
            h (length): Geometric altitude - aliases are 'alt', 'altitude'
            M (dimless): Velocity as Mach number - aliases are 'mach', 'Mach',
                'M_inf', 'mach_number'
            TAS (speed): Velocity as true airspeed - aliases are 'tas',
                'true_airspeed', 'U_inf', 'V_inf'
            CAS (speed): Velocity as calibrated airspeed - aliases are 'cas',
                'calibrated_airspeed'
            EAS (speed): Velocity as equivalent airspeed - aliases are 'eas',
                'equivalent_airspeed'
            L (length): Length scale - aliases are 'ell', 'bylen', 'length',
                'length_scale', 'l'
            Re (dimless): Reynolds number - alternative to velocity or length
                scale but not both - aliases are 'Re_L', 'reynolds_number'
            units (str): Set to 'US' for US units or 'SI' for SI
            full_output (bool): Set to True for full output
        """

        # Preprocess needed altitude-based quantities
        # Automatically process altitude through Atmosphere class
        self.byalt = Atmosphere(h=h, units=units, full_output=full_output,
                                **kwargs)

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
        Re_aliases = ['Re_L', 'reynolds_number']
        if Re is None:
            Re = __class__._arg_from_alias(Re_aliases, kwargs)

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

        # Compute velocity-based quantities
        self.byvel = Velocity(self.byalt)

        if M is not None:
            self.byvel.M = M
        elif TAS is not None:
            self.byvel.TAS = TAS
        elif CAS is not None:
            self.byvel.CAS = CAS
        elif EAS is not None:
            self.byvel.EAS = EAS
        elif Re is not None:
            # Velocity is set based on Re and L - set dummy speed for now
            self.byvel.M = 0
        else:
            # Use Mach=0 if no velocity is input
            self.byvel.M = 0*dimless

        # Check that computations are within valid Mach number limits
        M_ = np.atleast_1d(self.byvel.M)
        self._mach_min = 0 * dimless
        self._mach_max = 30 * dimless
        if (M_ < self._mach_min).any() or (self._mach_max < M_).any():
            raise ValueError(
                f"Mach number is out of bounds "
                f"({self._mach_min:.5g} < M_ < {self._mach_max:.5g})"
            )

        # Compute length-scale-based quantities
        self.bylen = Length(self.byalt, self.byvel)

        # If length scale is not input, default to unity with dimentionals unit
        # based on US or SI determination
        if L is None:
            if Re is None:
                L_unit = unit('ft') if self.units == 'US' else unit('m')
                self.bylen.L = 1.0 * L_unit
            else:
                Re *= dimless
                check_dimensionless(Re)
                self.bylen.L = NonDimensional.reynolds_number_length(
                    Re, self.byvel.TAS, self.byalt.nu)
        else:  # length scale is input by user
            self.bylen.L = L

            # If altitude is 0, but length scale is input, determine if US
            # units based on the input length scale
            if units == "":
                if self.byalt.h.magnitude.all() == 0:
                    if check_US_length_units(L):
                        self.units = 'US'

        # Velocity is set based on Re and L
        if Re is not None:
            if (M is not None
                    or TAS is not None
                    or CAS is not None
                    or EAS is not None):
                msg = "Overriding velocity based on Reynolds number."
                warnings.warn(msg)
            self.bylen.Re = Re  # velocity is set based on Re and L

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

    def tostring(self, full_output=True, units=None, pretty_print=True):
        """String representation of data structure.

        Args:
            full_output (bool): Set to True for full output
            units (str): Set to 'US' for US units or 'SI' for SI
            pretty_print (bool): Pretty print output

        Returns:
            str: String representation
        """
        # Determine full output flag
        if full_output is None:
            if self.full_output is None:
                full_output = True
            else:
                full_output = self.full_output

        # Determine units
        if units is not None:
            self.units = units

        # Determine maximum characters to add spaces for and assemble string
        max_var_chars = max([
            max([len(v) for v in Atmosphere.varnames.values()]),
            max([len(v) for v in Velocity.varnames.values()]),
            max([len(v) for v in Length.varnames.values()]),
        ])
        alti_str = self.byalt.tostring(full_output, self.units,
                                       pretty_print=pretty_print,
                                       max_var_chars=max_var_chars)
        spee_str = self.byvel.tostring(full_output, self.units,
                                       max_var_chars=max_var_chars,
                                       pretty_print=pretty_print)
        leng_str = self.bylen.tostring(full_output, self.units,
                                       max_var_chars=max_var_chars,
                                       pretty_print=pretty_print)

        unit_str = self.units
        ext_str = "full_output=True" if full_output else "full_output=False"
        top_hdr = f"   Flight Condition (units={unit_str}, {ext_str})"
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
    def units(self):
        """Unit system to use: 'SI', 'US', etc.  Available unit systems given
        by dir(unit.sys).

        Returns:
            str: Unit system
        """
        return self.byalt.units

    @units.setter
    def units(self, units):
        """Unit system to use: 'SI', 'US', etc.  Available unit systems given
        by dir(unit.sys).

        Args:
            units (str): Unit system
        """
        self.byalt.units = units

    @property
    def full_output(self):
        """Enable or disable full output of data by default.

        Returns:
            bool: Full output flag
        """
        return self.byalt.full_output

    @full_output.setter
    def full_output(self, full_output):
        """Unit system to use: 'SI', 'US', etc.  Available unit systems given
        by dir(unit.sys).

        Args:
            full_output (bool): Full output flag
        """
        self.byalt.full_output = full_output

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
        if self.units == 'US':
            redim_force.ito('lbf')
        elif self.units == 'SI':
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
        if self.units == 'US':
            redim_moment.ito('ft lbf')
        elif self.units == 'SI':
            redim_moment.ito('m N')

        return redim_moment
