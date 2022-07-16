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

from flightcondition.atmosphere import Atmosphere, DimensionalData
from flightcondition.constants import PhysicalConstants as Phys
from flightcondition.units import unit, dimless, check_area_dimensioned,\
    check_dimensioned, check_length_dimensioned, check_US_length_units,\
    to_base_units_wrapper


class Speed(DimensionalData):
    """Class to hold speed data. """

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

        M_str = f"M      = {self.M:10.5g{pp_}}"

        if US_units:
            TAS_str = f"TAS    = {self.TAS.to('knots'):10.5g{pp_}}"
            CAS_str = f"CAS    = {self.CAS.to('knots'):10.5g{pp_}}"
            EAS_str = f"EAS    = {self.EAS.to('knots'):10.5g{pp_}}"
            q_str   = f"q_inf  = {self.q_inf.to('lbf/ft^2'):10.5g{pp_}}"
        else:  # SI units
            TAS_str = f"TAS    = {self.TAS.to('m/s'):10.5g{pp_}}"
            CAS_str = f"CAS    = {self.CAS.to('m/s'):10.5g{pp_}}"
            EAS_str = f"EAS    = {self.EAS.to('m/s'):10.5g{pp_}}"
            q_str   = f"q_inf  = {self.q_inf.to('Pa'):10.5g{pp_}}"

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

        if full_output:
            repr_str = (f"{M_str}\n{TAS_str}\n{CAS_str}\n{EAS_str}\n{q_str}")
        else:
            repr_str = (f"{M_str}\n{TAS_str}\n{CAS_str}\n{EAS_str}")

        return repr_str


class Length(DimensionalData):
    """Class to hold length data. """

    ell: unit.Quantity

    varnames = {
        'ell': 'length_scale',
        'Re': 'reynolds_number',
        'Re_by_ell': 'reynolds_per_length',
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

        Re_str = f"Re     = {self.Re:10.5g~P}"
        if US_units:
            ell_str = f"ell    = {self.ell.to('ft'):10.5g~P}"
            Re_by_ell_str = f"Re/ell = {self.Re_by_ell.to('1/in'):10.5g~P}"
        else:  # SI units
            ell_str = f"ell    = {self.ell.to('m'):10.5g~P}"
            Re_by_ell_str = f"Re/ell = {self.Re_by_ell.to('1/mm'):10.5g~P}"

        # Insert longer variable name into output
        ell_str = f"{self.varnames['ell']:{max_var_chars}s} {ell_str}"
        Re_str = f"{self.varnames['Re']:{max_var_chars}s} {Re_str}"
        Re_by_ell_str = (f"{self.varnames['Re_by_ell']:{max_var_chars}s} "
                         f"{Re_by_ell_str}")

        if full_output:
            repr_str = (f"{ell_str}\n{Re_str}\n{Re_by_ell_str}")
        else:
            repr_str = (f"{ell_str}\n{Re_str}\n{Re_by_ell_str}")

        return repr_str


class FlightCondition:
    """Easily convert between Mach number, true airspeed (TAS), calibrated
    airspeed (CAS), and equivalent airspeed (EAS) for given altitude(s).
    Additional flight condition data and atmospheric data is computed.

    All inputs must be dimensional unit quantities.

    Usage:
        from flightcondition import FlightCondition, unit, dimless

        # Compute flight conditions for a scalar or array of altitudes
        altitudes = [0, 10, 33.5] * unit('kft')
        fc = FlightCondition(altitudes, EAS=300*unit('knots'))

        # Print all flight condition data:
        print(f"{fc}")

        # Print while specifying abbreviated output in SI units:
        print(f"\n{fc.tostring(full_output=False, US_units=False)}")

        # Access flight speed formats individually
        M = fc.speed.M
        U_inf, U_CAS, U_EAS = fc.speed.TAS, fc.speed.CAS, fc.speed.EAS

        # Access atmospheric data alone (see Atmosphere class for options)
        alt = fc.altitude  # access Atmosphere object 'altitude'
        p, T, rho, nu, a = alt.p, alt.T, alt.rho, alt.nu, alt.a

        # Input true/calibrated/equivalent airspeed or Mach number
        fc_TAS = FlightCondition(altitudes, TAS=300*unit('knots'))
        fc_CAS = FlightCondition(altitudes, CAS=300*unit('knots'))
        fc_EAS = FlightCondition(altitudes, EAS=300*unit('knots'))
        fc_M = FlightCondition(altitudes, M=0.4535*dimless)

        # Specify desired units on input and output
        altitudes_in_km = [0, 3.048, 10.2108] * unit('km')
        lengthscale = 60 * unit('in')  # arbitrary length scale of interest
        fc_other_units = FlightCondition(altitudes, EAS=154.33*unit('m/s'),
                                        ell=lengthscale)
        U_TAS = fc_other_units.speed.TAS
        print(f"\nThe true airspeed in m/s is {U_TAS.to('m/s'):.5g}")
        print(f"The true airspeed in km/s is {U_TAS.to('km/s'):.5g}")

        # Compute additional derived quantities (see class for all options)
        print(f"\nThe dynamic pressure in psi is "
            f"{fc.speed.q_inf.to('psi'):.5g}")
        print(f"The Reynolds number is {fc.length.Re:.5g}")
        print(f"The Reynolds number per-unit-length [1/in] is "
            f"{fc.length.Re_by_ell.to('1/in'):.5g}")
    """

    def __init__(
        self, h=0*unit('ft'), M=None, TAS=None, CAS=None, EAS=None, ell=None
    ):
        """Constructor based on altitude and input speed in terms of Mach
        number, TAS, CAS, or EAS.  Input at least one speed at the desired
        altitude.  All inputs must be dimensional unit quantities.

        Args:
            h (length): Geometric altitude
            M (dimless): Mach number
            TAS (speed): True airspeed
            CAS (speed): Calibrated airspeed
            EAS (speed): Equivalent airspeed
            ell (length): Length scale
        """
        # Set up classes for speed quantities and length quantities
        self.speed = Speed()
        self.length = Length()

        # Preprocess needed altitude-based quantities
        # Automatically process altitude through Atmosphere class
        self.altitude = Atmosphere(h)

        h0 = 0 * unit('kft')
        self._altitude0 = Atmosphere(h0)
        p_inf = self.altitude.p
        p_inf_h0 = self._altitude0.p
        self._delta = p_inf/p_inf_h0

        # Determine if US units or not
        self.US_units = self.altitude.US_units

        # Compute speed-based quantities
        # Use Mach=0 if no speed is input
        if M is None and TAS is None and CAS is None and EAS is None:
            M = 0*dimless

        if M is not None:
            self.speed.M = self._checkandsize(M, default_dim=dimless)
            self.speed.TAS = self._TAS_from_M(self.speed.M)
            self.speed.EAS = self._EAS_from_TAS(self.speed.TAS, self.speed.M)
            self.speed.q_c = self._q_c_from_M(self.speed.M)
            self.speed.CAS = self._CAS_from_q_c(self.speed.q_c)
        elif TAS is not None:
            self.speed.TAS = self._checkandsize(TAS)
            self.speed.M = self._M_from_TAS(TAS)
            self.speed.EAS = self._EAS_from_TAS(self.speed.TAS, self.speed.M)
            self.speed.q_c = self._q_c_from_M(self.speed.M)
            self.speed.CAS = self._CAS_from_q_c(self.speed.q_c)
        elif CAS is not None:
            self.speed.CAS = self._checkandsize(CAS)
            self.speed.q_c = self._q_c_from_CAS(self.speed.CAS)
            self.speed.M = self._M_from_q_c(self.speed.q_c)
            self.speed.TAS = self._TAS_from_M(self.speed.M)
            self.speed.EAS = self._EAS_from_TAS(self.speed.TAS, self.speed.M)
        elif EAS is not None:
            self.speed.EAS = self._checkandsize(EAS)
            self.speed.M = self._M_from_EAS(self.speed.EAS)
            self.speed.TAS = self._TAS_from_M(self.speed.M)
            self.speed.q_c = self._q_c_from_M(self.speed.M)
            self.speed.CAS = self._CAS_from_q_c(self.speed.q_c)
        else:
            raise TypeError("Input M, TAS, CAS, or EAS")
        self.speed.q_inf = self._q_inf_from_TAS(self.speed.TAS)

        # Check that computations are within valid Mach number range
        M = atleast_1d(self.speed.M)
        self._mach_min = 0 * dimless
        self._mach_max = 1 * dimless
        if (M < self._mach_min).any() or (self._mach_max < M).any():
            raise ValueError(
                f"Mach number is out of bounds "
                f"({self._mach_min:.5g} < M < {self._mach_max:.5g})"
            )

        # Compute length-scale-based quantities
        # If length scale is not input, default to unity with dimentionals unit
        # based on US or SI determination
        if ell is None:
            ell_unit = unit('ft') if self.US_units else unit('m')
            self.length.ell = 1.0 * ell_unit
        else:  # length scale is input by user
            self.length.ell = ell

            # Verify that length input is dimensional length quantity
            check_length_dimensioned(self.length.ell)

            # If altitude is 0, but length scale is input, determine if US
            # units based on the input length scale
            if self.altitude.h.magnitude.all() == 0:
                self.US_units = check_US_length_units(ell)

        # Assign lengths scale and compute quantities
        self.length.Re = self._reynolds_number(self.length.ell)
        self.length.Re_by_ell = self._reynolds_per_length()

        # Allow access to variables using full names
        self.speed._init_byvarname(self.speed.varnames)
        self.length._init_byvarname(self.length.varnames)

    def __str__(self):
        """Output string when object is printed.

        Returns:
            str: Full string output
        """
        return self.tostring(full_output=True)

    def __repr__(self):
        """Output string representation of class object.

        Returns:
            str: Full string output
        """
        return self.tostring(full_output=False)

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
        if shape(self.altitude.h):  # if altitude is an array
            if shape(inpvar):  # if inpvar is an array
                if not inpvar.size == self.altitude.h.size:
                    if not inpvar.size == 1 and not self.altitude.h.size == 1:
                        raise TypeError(wrongsize_msg)

            sizedarr = ones(shape(self.altitude.h))*inpvar
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
        alti_str = self.altitude.tostring(full_output, US_units,
                                          pretty_print=pretty_print)
        max_var_chars = max([
            max([len(v) for v in self.speed.varnames.values()]),
            max([len(v) for v in self.altitude.varnames.values()])
        ])
        spee_str = self.speed.tostring(full_output, US_units,
                                       max_var_chars=max_var_chars,
                                       pretty_print=pretty_print)
        leng_str = self.length.tostring(full_output, US_units,
                                        max_var_chars=max_var_chars,
                                        pretty_print=pretty_print)

        unit_str = "US" if US_units else "SI"
        ext_str = "extended" if full_output else "abbreviated"
        head_str = (f"    Flight Condition ({unit_str} units, {ext_str})")
        line_str = "========================================================="
        alti_hdr = "------------------ Altitude Quantities ------------------"
        spee_hdr = "------------------ Speed Quantities ---------------------"
        leng_hdr = "------------------ Length Quantities --------------------"

        repr_str = (f"{line_str}\n{head_str}\n{line_str}"
                    f"\n{alti_hdr}"
                    f"\n{alti_str}"
                    f"\n{spee_hdr}"
                    f"\n{spee_str}"
                    f"\n{leng_hdr}"
                    f"\n{leng_str}"
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
        a_inf = self.altitude.a
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
        a_inf = self.altitude.a
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
        a_inf_h0 = self._altitude0.a
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
        a_inf_h0 = self._altitude0.a
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
        rho_inf = self.altitude.rho
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
        a_h0 = self._altitude0.a
        p_h0 = self._altitude0.p
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
        a_h0 = self._altitude0.a
        p_h0 = self._altitude0.p
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
        p_inf = self.altitude.p
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
        p_inf = self.altitude.p
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
        a_inf_h0 = self._altitude0.a
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
    def redim_moment(self, S_ref, ell):
        """Factor to redimensionalize moment from moment coefficient e.g.,

            .. math::

                M & = C_M q_\\infty S_ref ell
                  & = C_M \\text{redim_moment}

        Args:
            S_ref (area): Reference area
            ell (length): Reference length

        Returns:
            moment: Redimensionalization factor
        """
        check_area_dimensioned(S_ref)
        check_length_dimensioned(ell)
        return self.q_inf * S_ref * ell

    @to_base_units_wrapper
    def _reynolds_number(self, ell):
        """Compute Reynolds number for the flight condition.

        Args:
            ell (length): Length scale

        Returns:
            dimless: Reynolds number
        """
        nu = self.altitude.nu
        TAS = self.speed.TAS
        Re_ell = TAS*ell/nu
        return Re_ell

    def _reynolds_per_length(self, length_unit='in'):
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

        Args:
            length_unit (length): Desired length unit as string, ('in', 'mm',
                'cm')

        Returns:
            dimless: Reynolds number
        """
        nu_inf = self.altitude.nu
        TAS = self.speed.TAS
        Re_by_length_unit = TAS/nu_inf
        Re_by_length_unit.ito(f"1/{length_unit}")
        return Re_by_length_unit
