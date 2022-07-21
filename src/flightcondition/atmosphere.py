#!/usr/bin/env python
"""Compute atmospheric properties from ICAO 1993 standard atmosphere model.

Dependencies: numpy, pint

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

from functools import wraps
from numpy import atleast_1d, exp, ndarray, pi, shape, size, str_, sqrt,\
    zeros_like

from flightcondition.constants import PhysicalConstants as Phys
from flightcondition.constants import AtmosphereConstants as Atmo
from flightcondition.units import unit, check_dimensioned,\
    check_length_dimensioned, check_US_length_units, to_base_units_wrapper


def _len1array_to_scalar(func):
    """Decorator to output scalar if array is length 1.

    Args:
        func (object): Function to wrap

    Returns:
        object: Scalar or array
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        output = func(*args, **kwargs)
        if isinstance(output.magnitude, ndarray):
            if size(output) == 1:
                return output[0]
        return output
    return wrapper


def _formatarr(unitstr):
    """Decorator to format arrays from base atmospheric model.

    Args:
        unitstr (str): Unit string

    Returns:
        callable: Function after wrapped operations complete
    """
    def inner_decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            arr = func(*args, **kwargs)
            arr *= unit(unitstr)
            if size(arr) == 1:
                arr = arr[0]
            return arr
        return wrapper
    return inner_decorator


def _property_decorators(func):
    """ Combine multiple decorators for properties of the class.

    Equivalent to:

        @property
        @to_base_units_wrapper
        @_len1array_to_scalar
    """
    return property(to_base_units_wrapper(_len1array_to_scalar(func)))
    # return _len1array_to_scalar(to_base_units_wrapper(property(func)))


class AccessByName():
    """Nested class to reference quantities by their full name. """

    def _populate_data(self, varsobj, varnames_dict):
        """Populate full names and link to variable
        Args:
            varnames_dict: Dictionary that maps variables to their longer
                names
            varsobj: Object that holds all of the variables
        """
        for var, varname in varnames_dict.items():
            setattr(self, varname, getattr(varsobj, var))


class DimensionalData:
    """Parent class to hold dimensional data"""

    def __init__(self):
        """Initialize object. """
        self.byname = AccessByName()

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


class Atmosphere(DimensionalData):
    """Compute quantities from International Civil Aviation Organization (ICAO)
    1993, which extends the US 1976 Standard Atmospheric Model to 80 km.

    Usage:
        from flightcondition import Atmosphere, unit

        # Compute atmospheric data for a scalar or array of altitudes
        h = [0.0, 44.2, 81.0] * unit('km')
        atm = Atmosphere(h)

        # Uncomment to print all atmospheric quantities:
        #print(f"\n{atm}")

        # Uncomment to print while specifying abbreviated output in US units:
        #print(f"\n{atm.tostring(full_output=False, US_units=True)}")

        # See also the linspace() function from numpy, e.g.
        # h = linspace(0, 81.0, 82) * unit('km')

        # Access individual properties and convert to desired units: "
        p, T, rho, nu, a, k = atm.p, atm.T, atm.rho, atm.nu, atm.a, atm.k
        print(f"\nThe pressure in psi is {p.to('psi'):.3g}")
        # >>> The pressure in psi is [14.7 0.024 0.000129] psi

        # Compute additional properties such as mean free path
        # Explore the class data structure for all options
        print( f"\nThe mean free path = {atm.MFP:.3g}")
        # >>> The mean free path = [7.25e-08 4.04e-05 0.00564] yd
    """

    varnames = {
        'h': 'geometric_altitude',
        'H': 'geopotential_altitude',
        'p': 'pressure',
        'T': 'temperature',
        'rho': 'density',
        'a': 'sound_speed',
        'mu': 'dynamic_viscosity',
        'nu': 'kinematic_viscosity',
        'k': 'thermal_conductivity',
        'g': 'gravity',
        'MFP': 'mean_free_path',
    }

    class Layer():
        """Nested class to compute and store layer data. """

        def __init__(self, H_arr):
            """Initialize Layer nested class.

            Args:
                H_arr (length): Geopotential altitude

            """
            H_arr = atleast_1d(H_arr)

            self._layer_name = [""]*size(H_arr)
            self._H_base = zeros_like(H_arr) * unit('m')
            self._T_base = zeros_like(H_arr) * unit('K')
            self._T_grad = zeros_like(H_arr) * unit('K/m')
            self._p_base = zeros_like(H_arr) * unit('Pa')

            for idx, H in enumerate(H_arr):
                jdx = __class__._layer_idx(H)
                self._layer_name[idx] = Atmo.layer_names[jdx]
                self._H_base[idx] = Atmo.H_base[jdx]
                self._T_base[idx] = Atmo.T_base[jdx]
                self._T_grad[idx] = Atmo.T_grad[jdx]
                self._p_base[idx] = Atmo.p_base[jdx]

        @staticmethod
        def _layer_idx(H):
            """Find index for layer data.

            Args:
                H_arr (length): Geopotential altitude

            Returns:
                int: Index of layer
            """
            idx = None
            for idx, H_base in enumerate(Atmo.H_base[:-1]):
                H_base_np1 = Atmo.H_base[idx+1]
                if H_base <= H < H_base_np1:
                    return idx
            else:  # no break
                raise ValueError("H out of bounds.")

        @property
        def name(self):
            """Layer name """
            if size(self._layer_name) == 1:
                layername = self._layer_name[0]
            else:
                layername = self._layer_name
            return layername

        @_property_decorators
        def H_base(self):
            """Layer base geopotential altitude :math:`H_{base}` """
            return self._H_base

        @_property_decorators
        def T_base(self):
            """Layer base temperature :math:`T_{base}` """
            return self._T_base

        @_property_decorators
        def T_grad(self):
            """Layer base temperature gradient :math:`T_{grad}` """
            return self._T_grad

        @_property_decorators
        def p_base(self):
            """Layer base pressure :math:`p_{base}` """
            return self._p_base

    def __init__(self, h=0*unit('kft')):
        """Input geometric altitude - object contains the corresponding
        atmospheric quantities.

        Args:
            h (length): Geometric altitude
        """
        # Compute altitude bounds
        self._H_min = Atmo.H_base[0]
        self._H_max = Atmo.H_base[-1]
        self._h_min = __class__.h_from_H(self._H_min)
        self._h_max = __class__.h_from_H(self._H_max)

        # Process altitude input
        self._h = self._process_input_altitude(h)
        self._H = self.H_from_h(self.h)
        self.layer = __class__.Layer(self.H)

        # Allow access to variables using full names
        super().__init__()
        self.byname._populate_data(varsobj=self, varnames_dict=self.varnames)

    def _process_input_altitude(self, alt):
        """Check that input is of type Quantity from pint package. Check that
        input is length dimension.  Check bounds.  Format as array even if
        scalar input.

        Args:
            alt (length): Input scalar or array of altitudes

        Returns:
            length: Geometric altitude
        """
        tofloat = 1.0
        h = atleast_1d(alt) * tofloat

        check_dimensioned(h)
        check_length_dimensioned(h)

        if len(shape(h)) > 1:
            raise TypeError("Input must be scalar or 1-D array.")

        if (h < self._h_min).any() or (self._h_max < h).any():
            raise ValueError(
                f"Input altitude is out of bounds "
                f"({self._h_min:.5g} < h < {self._h_max:.5g})"
            )

        self.US_units = check_US_length_units(h)
        if self.US_units:
            unit.default_system = 'US'  # 'US', 'imperial'

        return h

    @to_base_units_wrapper
    def knudsen_number(self, ell):
        """Compute the Knudsen number :math:`K_n`

        Args:
            ell (length): Length scale

        Returns:
            dimless: Knudsen number
        """
        Kn = self.MFP/ell
        return Kn

    def tostring(self, full_output=True, US_units=None, max_var_chars=0,
                 pretty_print=True):
        """String representation of data structure.

        Args:
            full_output (bool): Set to True for full output
            US_units (bool): Set to True for US units and False for SI
            max_var_chars (int): Maximum number of characters in unit string
            pretty_print (bool): Pretty print output

        Returns:
            str: String representation of class object
        """
        US_units = self.US_units if US_units is None else US_units
        pp_ = '~P' if pretty_print else ''
        if US_units:
            h_str   = f"h       = {self.h.to('kft'):10.5g{pp_}}"
            H_str   = f"H       = {self.H.to('kft'):10.5g{pp_}}"
            p_str   = f"p       = {self.p.to('lbf/ft^2'):10.5g{pp_}}"
            T_str   = f"T       = {self.T.to('degR'):10.5g{pp_}}"
            rho_str = f"rho     = {self.rho.to('slug/ft^3'):10.4e{pp_}}"
            a_str   = f"a       = {self.a.to('ft/s'):10.5g{pp_}}"
            mu_str  = f"mu      = {self.mu.to('lbf/ft^2 s'):10.4e{pp_}}"
            nu_str  = f"nu      = {self.nu.to('ft^2/s'):10.4e{pp_}}"
            k_str   = f"k       = {self.k.to('slug ft/s^3/degR'):10.4e{pp_}}"
            g_str   = f"g       = {self.g.to('ft/s^2'):10.5g{pp_}}"
            MFP_str = f"MFP     = {self.MFP.to('ft'):10.4e{pp_}}"
        else:  # SI units
            h_str   = f"h       = {self.h.to('km'):10.5g{pp_}}"
            H_str   = f"H       = {self.H.to('km'):10.5g{pp_}}"
            p_str   = f"p       = {self.p.to('Pa'):10.5g{pp_}}"
            T_str   = f"T       = {self.T.to('degK'):10.5g{pp_}}"
            rho_str = f"rho     = {self.rho.to('kg/m^3'):10.4e{pp_}}"
            a_str   = f"a       = {self.a.to('m/s'):10.5g{pp_}}"
            mu_str  = f"mu      = {self.mu.to('Pa s'):10.4e{pp_}}"
            nu_str  = f"nu      = {self.nu.to('m^2/s'):10.4e{pp_}}"
            k_str   = f"k       = {self.k.to('W/m/K'):10.4e{pp_}}"
            g_str   = f"g       = {self.g.to('m/s^2'):10.5g{pp_}}"
            MFP_str = f"MFP     = {self.MFP.to('m'):10.4e{pp_}}"

        # Insert longer variable name into output
        max_var_chars = max([
            max([len(v) for v in self.varnames.values()]),
            max_var_chars
        ])
        h_str   = f"{self.varnames['h']:{max_var_chars}s} {h_str}"
        H_str   = f"{self.varnames['H']:{max_var_chars}s} {H_str}"
        p_str   = f"{self.varnames['p']:{max_var_chars}s} {p_str}"
        T_str   = f"{self.varnames['T']:{max_var_chars}s} {T_str}"
        rho_str = f"{self.varnames['rho']:{max_var_chars}s} {rho_str}"
        a_str   = f"{self.varnames['a']:{max_var_chars}s} {a_str}"
        mu_str  = f"{self.varnames['mu']:{max_var_chars}s} {mu_str}"
        nu_str  = f"{self.varnames['nu']:{max_var_chars}s} {nu_str}"
        k_str   = f"{self.varnames['k']:{max_var_chars}s} {k_str}"
        g_str   = f"{self.varnames['g']:{max_var_chars}s} {g_str}"
        MFP_str = f"{self.varnames['MFP']:{max_var_chars}s} {MFP_str}"

        if full_output:
            if type(self.layer.name) is str_:  # singular string
                trunc_layer_name = self.layer.name
            else:
                trunc_layer_name = "[" + " ".join([
                    f"{s[:10]}" for s in self.layer.name
                ]) + "]"
            layer_str = (f"{'atmospheric_layer ':{max_var_chars}} name    = "
                         f"{trunc_layer_name}")

            repr_str = (f"{h_str}\n{H_str}\n{p_str}\n{T_str}\n{rho_str}\n"
                        f"{a_str}\n{mu_str}\n{nu_str}\n{k_str}\n{g_str}\n"
                        f"{MFP_str}\n{layer_str}")
        else:
            repr_str = (f"{h_str}\n{p_str}\n{T_str}\n{rho_str}\n{a_str}\n"
                        f"{nu_str}")
        return repr_str

    @staticmethod
    def H_from_h(h):
        """Convert geometric to geopotential altitude.

        :math:`H = \\frac{R_{earth} h}{R_{earth} + h}`

        Args:
            h (length): Geometric altitude

        Returns:
            length: Geopotential altitude
        """
        R_earth = Phys.R_earth
        H = R_earth*h/(R_earth + h)
        return H

    @staticmethod
    def h_from_H(H):
        """Convert geopotential to geometric altitude.

        :math:`h = \\frac{R_{earth} H}{R_{earth} - H}`

        Args:
            H (length): Geopotential altitude

        Returns:
            length: Geometric altitude
        """
        R_earth = Phys.R_earth
        h = R_earth*H/(R_earth - H)
        return h

    @_property_decorators
    def h(self):
        """Geometric altitude :math:`h` """
        return self._h

    @_property_decorators
    def H(self):
        """Geopotential altitude :math:`H` """
        return self._H

    @_property_decorators
    def p(self):
        """Air pressure :math:`p` """
        H_base = atleast_1d(self.layer.H_base)
        T_base = atleast_1d(self.layer.T_base)
        T_grad = atleast_1d(self.layer.T_grad)
        p_base = atleast_1d(self.layer.p_base)

        H = atleast_1d(self.H)
        T = atleast_1d(self.T)
        g_0 = Phys.g
        R_air = Phys.R_air

        p = zeros_like(H) * unit('Pa')

        # Pressure equation changes between T_grad == 0 and T_grad != 0
        s = T_grad == 0
        p[s] = p_base[s]*exp((-g_0/(R_air*T[s]))*(H[s] - H_base[s]))

        s = T_grad != 0
        p[s] = p_base[s]*(
            1 + (T_grad[s]/T_base[s])*(H[s] - H_base[s])
        )**((1/T_grad[s])*(-g_0/R_air))

        return p

    @_property_decorators
    def T(self):
        """Ambient air temperature :math:`T` """
        T_grad = atleast_1d(self.layer.T_grad)
        H_base = atleast_1d(self.layer.H_base)
        T_base = atleast_1d(self.layer.T_base)
        T = T_base + T_grad*(self.H - H_base)
        return T

    @_property_decorators
    def rho(self):
        """Ambient air density :math:`\\rho` """
        p = self.p
        T = self.T
        R_air = Phys.R_air
        rho = p/(R_air*T)
        return rho

    @_property_decorators
    def a(self):
        """Ambient speed of sound :math:`a` """
        T = self.T
        gamma_air = Phys.gamma_air
        R_air = Phys.R_air
        a = sqrt(gamma_air*R_air*T)
        return a

    @_property_decorators
    def mu(self):
        """Ambient dynamic viscosity :math:`\\mu` """
        T = self.T
        beta_s = Atmo.beta_s
        S = Atmo.S
        mu = beta_s*T**1.5/(T + S)
        return mu

    @_property_decorators
    def nu(self):
        """Ambient kinematic viscosity :math:`\\nu` """
        mu = self.mu
        rho = self.rho
        nu = mu/rho
        return nu

    @_property_decorators
    def k(self):
        """Ambient thermal conductivity :math:`k` """
        # Empirical formula requires T in Kelvin
        T_K = self.T.to('K').magnitude
        c1 = 0.002648151
        k_ = c1*T_K**1.5 / (T_K + (245.4*10**(-12/T_K)))

        # Re-assign to dimensioned variable
        k = k_ * unit('W/m/K')
        return k

    @_property_decorators
    def g(self):
        """Gravitational acceleration at altitude :math:`g` """
        h = self.h
        g_0 = Phys.g
        R_earth = Phys.R_earth
        g = g_0*(R_earth/(R_earth + h))**2
        return g

    @_property_decorators
    def MFP(self):
        """Mean free path """
        p = self.p
        T = self.T
        N_A = Phys.N_A
        R = Phys.R
        sigma = Phys.collision_diam_air

        n = N_A*p/(R*T)  # number density
        MFP = 1/(sqrt(2)*pi*sigma**2*n)
        return MFP
