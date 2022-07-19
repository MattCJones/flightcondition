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

from flightcondition.atmosphere import Atmosphere, AccessByName,\
    DimensionalData
from flightcondition.constants import PhysicalConstants as Phys
from flightcondition.units import unit, dimless, check_area_dimensioned,\
    check_dimensioned, check_length_dimensioned, check_US_length_units,\
    to_base_units_wrapper


class Airspeed(DimensionalData):
    """Class to hold airspeed data. """

    M: unit.Quantity
    TAS: unit.Quantity
    CAS: unit.Quantity
    EAS: unit.Quantity
    q_inf: unit.Quantity
    q_c: unit.Quantity

    varnames = {
        'M': 'mach_number',
        'TAS': 'true_airspeed',
        'CAS': 'calibrated_airspeed',
        'EAS': 'equivalent_airspeed',
        'q_inf': 'dynamic_pressure',
        'q_c': 'impact_pressure',
        'Re_by_L': 'reynolds_per_length',
    }

    def tostring(self, full_output=True, US_units=None, max_var_chars=0,
                 pretty_print=True):
        """String representation of data structure.

        Args:
            full_output (bool): Set to True for full output
            US_units (bool): Set to True for US units and False for SI
            max_var_chars (int): Maximum number of characters in unit string
            pretty_print (bool): Pretty print output

        Returns:
            str: String representation
        """
        pp_ = '~P' if pretty_print else ''

        M_str = f"M       = {self.M:10.5g{pp_}}"

        if US_units:
            TAS_str = f"TAS     = {self.TAS.to('knots'):10.5g{pp_}}"
            CAS_str = f"CAS     = {self.CAS.to('knots'):10.5g{pp_}}"
            EAS_str = f"EAS     = {self.EAS.to('knots'):10.5g{pp_}}"
            q_str   = f"q_inf   = {self.q_inf.to('lbf/ft^2'):10.5g{pp_}}"
            Re_by_L_str = f"Re_by_L = {self.Re_by_L.to('1/in'):10.5g~P}"
        else:  # SI units
            TAS_str = f"TAS     = {self.TAS.to('m/s'):10.5g{pp_}}"
            CAS_str = f"CAS     = {self.CAS.to('m/s'):10.5g{pp_}}"
            EAS_str = f"EAS     = {self.EAS.to('m/s'):10.5g{pp_}}"
            q_str   = f"q_inf   = {self.q_inf.to('Pa'):10.5g{pp_}}"
            Re_by_L_str = f"Re_by_L = {self.Re_by_L.to('1/mm'):10.5g~P}"

        # Insert longer variable name into output
        max_var_chars = max([
            max([len(v) for v in self.varnames.values()]),
            max_var_chars
        ])
        M_str   = f"{self.varnames['M']:{max_var_chars}s} {M_str}"
        TAS_str = f"{self.varnames['TAS']:{max_var_chars}s} {TAS_str}"
        CAS_str = f"{self.varnames['CAS']:{max_var_chars}s} {CAS_str}"
        EAS_str = f"{self.varnames['EAS']:{max_var_chars}s} {EAS_str}"
        q_str   = f"{self.varnames['q_inf']:{max_var_chars}s} {q_str}"
        Re_by_L_str = (f"{self.varnames['Re_by_L']:{max_var_chars}s} "
                       f"{Re_by_L_str}")

        if full_output:
            repr_str = (f"{M_str}\n{TAS_str}\n{CAS_str}\n{EAS_str}\n{q_str}"
                        f"\n{Re_by_L_str}")
        else:
            repr_str = (f"{M_str}\n{TAS_str}\n{CAS_str}\n{EAS_str}"
                        f"\n{Re_by_L_str}")

        return repr_str


class Length(DimensionalData):
    """Class to hold length data. """

    L: unit.Quantity

    varnames = {
        'L': 'length_scale',
        'Re': 'reynolds_number',
    }

    def tostring(self, full_output=True, US_units=None, max_var_chars=0,
                 pretty_print=True):
        """String representation of data structure.

        Args:
            full_output (bool): Set to True for full output
            US_units (bool): Set to True for US units and False for SI
            max_var_chars (int): Maximum number of characters in unit string
            pretty_print (bool): Pretty print output

        Returns:
            str: String representation
        """

        Re_str = f"Re      = {self.Re:10.5g~P}"
        if US_units:
            L_str = f"L       = {self.L.to('ft'):10.5g~P}"
        else:  # SI units
            L_str = f"L       = {self.L.to('m'):10.5g~P}"

        # Insert longer variable name into output
        max_var_chars = max([
            max([len(v) for v in self.varnames.values()]),
            max_var_chars
        ])
        L_str = f"{self.varnames['L']:{max_var_chars}s} {L_str}"
        Re_str = f"{self.varnames['Re']:{max_var_chars}s} {Re_str}"

        if full_output:
            repr_str = (f"{L_str}\n{Re_str}")
        else:
            repr_str = (f"{L_str}\n{Re_str}")

        return repr_str


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
        #print(f"\n{fc.tostring(full_output=False, US_units=True)}")

        # Access true, calibrated, equivalent airspeeds
        KTAS = fc.vel.TAS.to('knots')
        KCAS = fc.vel.CAS.to('knots')
        KEAS = fc.vel.EAS.to('knots')
        print(f"Flying at {KTAS.magnitude:.4g} KTAS,"
            f" which is {KCAS.magnitude:.4g} KCAS,"
            f" or {KEAS.magnitude:.4g} KEAS")
        # >>> Flying at 319.4 KTAS, which is 277.7 KCAS, or 275.1 KEAS

        # Access atmospheric data (see Atmosphere class for more)
        atm = fc.atm  # access Atmosphere object
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
            f"{fc.vel.q_inf.to('psi'):.3g}")
        # >>> The dynamic pressure in psi is [1.78 1.78 1.78] psi
        print(f"The Reynolds number is {fc.len.Re:.3g}")
        # >>> The Reynolds number is [9.69e+07 8.82e+07 7.95e+07]

        # Alternatively access quantities by their full name
        print(fc.vel.TAS == fc.byname.true_airspeed)
        # >>> [ True  True  True]
    """

    def __init__(
        self, h=0*unit('ft'), M=None, TAS=None, CAS=None, EAS=None, L=None
    ):
        """Constructor based on altitude and input airspeed in terms of Mach
        number, TAS, CAS, or EAS.  Input at least one format of airspeed at the
        desired altitude.  All inputs must be dimensional unit quantities.

        Args:
            h (length): Geometric altitude
            M (dimless): Mach number
            TAS (speed): True airspeed
            CAS (speed): Calibrated airspeed
            EAS (speed): Equivalent airspeed
            L (length): Length scale
        """
        # Set up classes for airspeed quantities and length quantities
        self.vel = Airspeed()
        self.len = Length()

        # Preprocess needed altitude-based quantities
        # Automatically process altitude through Atmosphere class
        self.atm = Atmosphere(h)

        h0 = 0 * unit('kft')
        self._atm0 = Atmosphere(h0)
        p_inf = self.atm.p
        p_inf_h0 = self._atm0.p
        self._delta = p_inf/p_inf_h0

        # Determine if US units or not
        self.US_units = self.atm.US_units

        # Compute airspeed-based quantities
        # Use Mach=0 if no airspeed is input
        if M is None and TAS is None and CAS is None and EAS is None:
            M = 0*dimless

        if M is not None:
            self.vel.M = self._checkandsize(M, default_dim=dimless)
            self.vel.TAS = self._TAS_from_M(self.vel.M)
            self.vel.EAS = self._EAS_from_TAS(self.vel.TAS, self.vel.M)
            self.vel.q_c = self._q_c_from_M(self.vel.M)
            self.vel.CAS = self._CAS_from_q_c(self.vel.q_c)
        elif TAS is not None:
            self.vel.TAS = self._checkandsize(TAS)
            self.vel.M = self._M_from_TAS(TAS)
            self.vel.EAS = self._EAS_from_TAS(self.vel.TAS, self.vel.M)
            self.vel.q_c = self._q_c_from_M(self.vel.M)
            self.vel.CAS = self._CAS_from_q_c(self.vel.q_c)
        elif CAS is not None:
            self.vel.CAS = self._checkandsize(CAS)
            self.vel.q_c = self._q_c_from_CAS(self.vel.CAS)
            self.vel.M = self._M_from_q_c(self.vel.q_c)
            self.vel.TAS = self._TAS_from_M(self.vel.M)
            self.vel.EAS = self._EAS_from_TAS(self.vel.TAS, self.vel.M)
        elif EAS is not None:
            self.vel.EAS = self._checkandsize(EAS)
            self.vel.M = self._M_from_EAS(self.vel.EAS)
            self.vel.TAS = self._TAS_from_M(self.vel.M)
            self.vel.q_c = self._q_c_from_M(self.vel.M)
            self.vel.CAS = self._CAS_from_q_c(self.vel.q_c)
        else:
            raise TypeError("Input M, TAS, CAS, or EAS")
        self.vel.q_inf = self._q_inf_from_TAS(self.vel.TAS)

        # Check that computations are within valid Mach number range
        M = atleast_1d(self.vel.M)
        self._mach_min = 0 * dimless
        self._mach_max = 1 * dimless
        if (M < self._mach_min).any() or (self._mach_max < M).any():
            raise ValueError(
                f"Mach number is out of bounds "
                f"({self._mach_min:.5g} < M < {self._mach_max:.5g})"
            )

        # Compute Reynolds Number per length
        self.vel.Re_by_L = self._reynolds_per_length()

        # Compute length-scale-based quantities
        # If length scale is not input, default to unity with dimentionals unit
        # based on US or SI determination
        if L is None:
            L_unit = unit('ft') if self.US_units else unit('m')
            self.len.L = 1.0 * L_unit
        else:  # length scale is input by user
            self.len.L = L

            # Verify that length input is dimensional length quantity
            check_length_dimensioned(self.len.L)

            # If altitude is 0, but length scale is input, determine if US
            # units based on the input length scale
            if self.atm.h.magnitude.all() == 0:
                self.US_units = check_US_length_units(L)

        # Assign lengths scale and compute quantities
        self.len.Re = self._reynolds_number(self.len.L)

        # Allow access to variables by their full names in vel and len objects
        super().__init__()
        self.vel.byname._populate_data(varsobj=self.vel,
                                       varnames_dict=self.vel.varnames)
        self.len.byname._populate_data(varsobj=self.len,
                                       varnames_dict=self.len.varnames)

        # Access full names from base object
        self.byname._populate_data(varsobj=self.atm,
                                   varnames_dict=self.atm.varnames)
        self.byname._populate_data(varsobj=self.vel,
                                   varnames_dict=self.vel.varnames)
        self.byname._populate_data(varsobj=self.len,
                                   varnames_dict=self.len.varnames)

        # Access all variables from .byvar. at base level
        self.byvar = AccessByName()

        byvar_atm = {}
        for var, _ in self.atm.varnames.items():
            byvar_atm[var] = var
        self.byvar._populate_data(self.atm, byvar_atm)

        byvar_vel = {}
        for var, _ in self.vel.varnames.items():
            byvar_vel[var] = var
        self.byvar._populate_data(self.vel, byvar_vel)

        byvar_len = {}
        for var, _ in self.len.varnames.items():
            byvar_len[var] = var
        self.byvar._populate_data(self.len, byvar_len)

    @staticmethod
    def _isentropic_mach(p_0, p):
        """Isentropic flow equation for Mach number

        Args:
            p_0 (pressure): Stagnation pressure
            p (pressure): Static pressure

        Returns:
            dimless: Mach number
        """
        y = Phys.gamma_air
        M = sqrt((2/(y-1))*((p_0/p + 1)**((y-1)/y) - 1))
        return M

    @staticmethod
    def _isentropic_stagnation_pressure(M, p):
        """Isentropic flow equation for stagnation pressure

        Args:
            M (dimless): Mach number
            p (pressure): Static pressure

        Returns:
            pressure: Stagnation pressure
        """
        y = Phys.gamma_air
        p_0 = p*(-1 + (1 + ((y-1)/2)*M**2)**(y/(y-1)))
        return p_0

    @staticmethod
    def _mach_number(U, a):
        """Compute Mach number

        Args:
            U (speed): Airspeed
            a (speed): Sound speed

        Returns:
            dimless: Mach number
        """
        return U/a

    @staticmethod
    def _mach_airspeed(M, a):
        """Compute airspeed from Mach number

        Args:
            M (dimless): Mach number
            a (speed): Sound speed

        Returns:
            speed: Airspeed
        """
        return M*a

    def _checkandsize(self, inpvar, default_dim=None):
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
        # Require altitude and speed arrays to be equal or singular speed array
        wrongsize_msg = ("Input arrays for altitude and airspeed must either "
                         "equal in size or one must be singular.")
        if shape(self.atm.h):  # if altitude is an array
            if shape(inpvar):  # if inpvar is an array
                if not inpvar.size == self.atm.h.size:
                    if not inpvar.size == 1 and not self.atm.h.size == 1:
                        raise TypeError(wrongsize_msg)

            sizedarr = ones(shape(self.atm.h))*inpvar
        else:
            sizedarr = inpvar

        return sizedarr

    def print(self, *args, **kwargs):
        """Print tostring() function to stdout. """
        print(self.tostring(*args, **kwargs))

    def tostring(self, full_output=True, US_units=None, pretty_print=True):
        """String representation of data structure.

        Args:
            full_output (bool): Set to True for full output
            US_units (bool): Set to True for US units and False for SI
            pretty_print (bool): Pretty print output

        Returns:
            str: String representation
        """
        US_units = self.US_units if US_units is None else US_units
        alti_str = self.atm.tostring(full_output, US_units,
                                     pretty_print=pretty_print)
        max_var_chars = max([
            max([len(v) for v in self.vel.varnames.values()]),
            max([len(v) for v in self.atm.varnames.values()])
        ])
        spee_str = self.vel.tostring(full_output, US_units,
                                     max_var_chars=max_var_chars,
                                     pretty_print=pretty_print)
        leng_str = self.len.tostring(full_output, US_units,
                                     max_var_chars=max_var_chars,
                                     pretty_print=pretty_print)

        unit_str = "US" if US_units else "SI"
        ext_str = "extended" if full_output else "abbreviated"
        head_str = (f"    Flight Condition ({unit_str} units, {ext_str})")
        line_str = "========================================================="
        alti_hdr = "-----------------  Altitude Quantities  -----------------"
        spee_hdr = "-----------------  Airspeed Quantities  -----------------"
        leng_hdr = "-----------------   Length Quantities   -----------------"

        repr_str = (f"{line_str}\n{head_str}\n{line_str}"
                    f"\n{alti_hdr}" f"\n{alti_str}"
                    f"\n{spee_hdr}" f"\n{spee_str}"
                    f"\n{leng_hdr}" f"\n{leng_str}"
                    )

        return repr_str

    @to_base_units_wrapper
    def _M_from_TAS(self, TAS):
        """Compute Mach number from true airspeed.

        Args:
            TAS (speed): True airspeed

        Returns:
            dimless: Mach number
        """
        a_inf = self.atm.a
        M = __class__._mach_number(U=TAS, a=a_inf)
        return M

    @to_base_units_wrapper
    def _TAS_from_M(self, M):
        """Compute true airspeed from Mach number.

        Args:
            M (dimless): Mach number

        Returns:
            speed: True airspeed
        """
        a_inf = self.atm.a
        TAS = __class__._mach_airspeed(M, a_inf)
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
        M = EAS/(a_inf_h0*sqrt(self._delta))
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
        EAS = M*(a_inf_h0*sqrt(self._delta))
        return EAS

    @to_base_units_wrapper
    def _q_inf_from_TAS(self, TAS):
        """Compute dynamic pressure from true airspeed.

        Args:
            TAS (speed): True airspeed

        Returns:
            pressure: Dynamic pressure
        """
        rho_inf = self.atm.rho
        q_inf = 0.5*rho_inf*TAS**2
        return q_inf

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
        M_ = __class__._mach_number(U=CAS, a=a_h0)
        q_c = __class__._isentropic_stagnation_pressure(M=M_, p=p_h0)
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
        M_ = __class__._isentropic_mach(p_0=q_c, p=p_h0)
        CAS = __class__._mach_airspeed(M_, a_h0)
        return CAS

    @to_base_units_wrapper
    def _M_from_q_c(self, q_c):
        """Compute Mach number from impact pressure.

        Args:
            q_c (pressure): Impact pressure

        Returns:
            dimless: Mach number
        """
        p_inf = self.atm.p
        # Isentropic flow equation
        M = __class__._isentropic_mach(p_0=q_c, p=p_inf)

        return M

    @to_base_units_wrapper
    def _q_c_from_M(self, M):
        """Compute impact pressure from Mach number.

        Args:
            M: Mach number

        Returns:
            impact pressure
        """
        p_inf = self.atm.p
        # Solve for impact pressure from isentropic flow equation:
        q_c = __class__._isentropic_stagnation_pressure(M=M, p=p_inf)
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
        EAS = M*(a_inf_h0*sqrt(self._delta))
        return EAS

    @to_base_units_wrapper
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
        return self.q_inf * S_ref

    @to_base_units_wrapper
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
        return self.q_inf * S_ref * L

    @to_base_units_wrapper
    def _reynolds_number(self, L):
        """Compute Reynolds number for the flight condition.

        Args:
            L (length): Length scale

        Returns:
            dimless: Reynolds number
        """
        nu = self.atm.nu
        TAS = self.vel.TAS
        Re_L = TAS*L/nu
        return Re_L

    def _reynolds_per_length(self, length_unit='in'):
        """Compute Reynolds number divided by length unit.

        Reynolds number is,
            Re = U_inf * L / nu

        For some fluid dynamics solvers the user must input the Reynolds number
        in terms of Reynolds number per-unit-length,
            Re_by_L = Re_L / L_in_grid_units
                      = U_inf * (L/L_in_grid_units) / nu

        Assume the length scale is 5 feet.  When grid units are the same as
        standard length units, such as feet vs. feet, the term
        (L_in_grid_units/L) is unity, and does not change the Reynolds
        number magnitude.  For example,
            L/L_in_grid_units = (5 ft)/(5 ft) = 1 ft/ft

        When grid units differ from standard length units, such as inches vs.
        feet, the term (L_in_grid_units/L) becomes,
            L/L_in_grid_units = (5 ft)/(60 in) = 1/12 ft/in
        which is the conversion factor between feet and inches.

        So, the Re_by_L term becomes
            Re_by_L = Re_L * (1/12 1/in)

        In this case, it is helpful to compute Reynolds number per-unit-length
        in inches (length_unit='in').

        Args:
            length_unit (length): Desired length unit as string, ('in', 'mm',
                'cm')

        Returns:
            dimless: Reynolds number
        """
        nu_inf = self.atm.nu
        TAS = self.vel.TAS
        Re_by_length_unit = TAS/nu_inf
        Re_by_length_unit.ito(f"1/{length_unit}")
        return Re_by_length_unit
