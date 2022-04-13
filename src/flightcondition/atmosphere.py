#!/usr/bin/env python
"""Compute atmospheric properties from ICAO 1993 standard atmosphere model.

Dependencies: numpy, pint

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

from functools import wraps
from numpy import atleast_1d, array, exp, ndarray, pi, shape, size, sqrt,\
        zeros_like

from flightcondition.constants import PhysicalConstants as Phys
from flightcondition.constants import AtmosphereConstants as Atmo
from flightcondition.units import unit, check_dimensioned,\
    check_US_length_units, check_length_dimensioned,\
    to_base_units_wrapper


def _atleast_1d(arr):
    """DEPRECATED: My version of numpy.atleast_1d that supports Python 3.7.
    This function will be deleted in future versions.

    :arr: array
    :returns: TODO

    """
    if len(shape(arr)) == 0:  # scalar, non-array
        return array([arr.magnitude]) * arr.units
    else:  # already an array
        return arr


def _len1array_to_scalar(func):
    """Decorator to output scalar if array is length 1."""

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

    :arr: input array from base atmospheric model
    :unitstr: unit string
    :returns: dimensionalized, formatted array

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


class Atmosphere():

    """Compute quantities from International Civil Aviation Organization (ICAO)
    1993, which extends the US 1976 Standard Atmospheric Model to 80 km.

    Usage:

        from flightcondition import Atmosphere, unit

        # Compute atmospheric data for a scalar or array of altitudes
        h = [0.0, 12.7, 44.2, 81.0] * unit('km')
        atm = Atmosphere(h)

        # Print abbreviated output:
        print(f"\n{atm}")

        # Print extended output in US units:
        print(f"\n{atm.tostring(full_output=True, US_units=True)}")

        # See also the linspace() function from numpy, e.g.
        # h = linspace(0, 81.0, 82) * unit('km')

        # Access individual properties and convert to desired units: "
        p, T, rho, nu, a = atm.p, atm.T, atm.rho, atm.nu, atm.a
        print(f"\nThe pressure in psi is {p.to('psi'):.5g}")

        # Compute additional properties such as thermal conductivity,
        # mean free path, and more (see class for all options)
        print(f"\nThe thermal conductivity is {atm.k:.5g}"
            f"\nThe mean free path = {atm.mean_free_path:.5g}")
    """

    class Layer():
        """Nested class to compute and store layer data."""

        def __init__(self, H_arr):
            """Initialize Layer nested class.

            :H_arr: geopotential altitude
            """
            H_arr = _atleast_1d(H_arr)

            self._name = [""]*size(H_arr)
            self._H_base = zeros_like(H_arr) * unit('m')
            self._T_base = zeros_like(H_arr) * unit('K')
            self._T_grad = zeros_like(H_arr) * unit('K/m')
            self._p_base = zeros_like(H_arr) * unit('Pa')

            for idx, H in enumerate(H_arr):
                jdx = __class__._layer_idx(H)
                self._name[idx] = Atmo.layer_names[jdx]
                self._H_base[idx] = Atmo.H_base[jdx]
                self._T_base[idx] = Atmo.T_base[jdx]
                self._T_grad[idx] = Atmo.T_grad[jdx]
                self._p_base[idx] = Atmo.p_base[jdx]

        @staticmethod
        def _layer_idx(H):
            """Find index for layer data.

            :H_arr: geopotential altitude
            :returns: index

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
            """Layer name"""
            if size(self._name) == 1:
                layer_name = self._name[0]
            else:
                layer_name = self._name
            return layer_name

        @property
        @to_base_units_wrapper
        @_len1array_to_scalar
        def H_base(self):
            """Layer base geopotential height :math:`H_{base}`"""
            return self._H_base

        @property
        @to_base_units_wrapper
        @_len1array_to_scalar
        def T_base(self):
            """Layer base temperature :math:`T_{base}`"""
            return self._T_base

        @property
        @to_base_units_wrapper
        @_len1array_to_scalar
        def T_grad(self):
            """Layer base temperature gradient :math:`T_{grad}`"""
            return self._T_grad

        @property
        @to_base_units_wrapper
        @_len1array_to_scalar
        def p_base(self):
            """Layer base pressure :math:`p_{base}`"""
            return self._p_base

    def __init__(self, altitude):
        """Input altitude - object contains the corresponding atmospheric
        quantities.

        :altitude: geometric altitude input in pint dimensional units

        """

        # Compute altitude bounds
        self._H_min = Atmo.H_base[0]
        self._H_max = Atmo.H_base[-1]
        self._h_min = __class__.h_from_H(self._H_min)
        self._h_max = __class__.h_from_H(self._H_max)

        # Process altitude input
        self._h = self._process_input_altitude(altitude)
        self._H = self.H_from_h(self.h)
        self.layer = __class__.Layer(self.H)

    def __str__(self):
        """Output string representation of class object. Default for __str__
        :returns: string output

        """
        return self.tostring(full_output=False)

    def _process_input_altitude(self, alt):
        """Check that input is of type Quantity from pint package. Check that
        input is length dimension.  Check bounds.  Format as array even if
        scalar input.

        :alt: input scalar or array of altitudes
        :returns: geometric altitude

        """
        tofloat = 1.0
        h = _atleast_1d(alt) * tofloat

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
        :returns: Knudsen number
        """
        Kn = self.mean_free_path/ell
        return Kn

    def tostring(self, full_output=True, US_units=None):
        """String representation of data structure.

        :full_output: set to True for full output
        :US_units: set to True for US units and False for SI
        :returns: string representation

        """
        US_units = self.US_units if US_units is None else US_units
        layer_str = f"layer = {self.layer.name}"
        if US_units:
            h_str = f"h = {self.h.to('kft'):8.5g~P}"
            H_str = f"H = {self.H.to('kft'):8.5g~P}"
            p_str = f"p = {self.p.to('lbf/ft^2'):8.5g~P}"
            T_str = f"T = {self.T.to('degR'):8.5g~P}"
            rho_str = f"rho = {self.rho.to('slug/ft^3'):8.5g~P}"
            a_str = f"a = {self.a.to('ft/s'):8.5g~P}"
            mu_str = f"mu = {self.mu.to('lbf/ft^2 s'):8.5g~P}"
            nu_str = f"nu = {self.nu.to('ft^2/s'):8.5g~P}"
            k_str = f"k = {self.k.to('slug ft/s^3/degR'):8.5g~P}"
            g_str = f"g = {self.g.to('ft/s^2'):8.5g~P}"
            MFP_str = f"mean_free_path = {self.mean_free_path.to('ft'):8.5g~P}"
        else:  # SI units
            h_str = f"h = {self.h.to('km'):8.5g~P}"
            H_str = f"H = {self.H.to('km'):8.5g~P}"
            p_str = f"p = {self.p.to('Pa'):8.5g~P}"
            T_str = f"T = {self.T.to('degK'):8.5g~P}"
            rho_str = f"rho = {self.rho.to('kg/m^3'):8.5g~P}"
            a_str = f"a = {self.a.to('m/s'):8.5g~P}"
            mu_str = f"mu = {self.mu.to('Pa s'):8.5g~P}"
            nu_str = f"nu = {self.nu.to('m^2/s'):8.5g~P}"
            k_str = f"k = {self.k.to('W/m/K'):8.5g~P}"
            g_str = f"g = {self.g.to('m/s^2'):8.5g~P}"
            MFP_str = f"mean_free_path = {self.mean_free_path.to('m'):8.5g~P}"

        if full_output:
            repr_str = (f"{h_str}\n{H_str}\n{layer_str}\n{p_str}\n{T_str}\n"
                        f"{rho_str}\n{a_str}\n{mu_str}\n{nu_str}\n{k_str}\n"
                        f"{g_str}\n{MFP_str}")
        else:
            repr_str = (f"{h_str}\n{p_str}\n{T_str}\n{rho_str}\n{a_str}\n"
                        f"{nu_str}")
        return repr_str

    @staticmethod
    def H_from_h(h):
        """Convert geometric to geopotential altitude.

        :math:`H = \\frac{R_{earth} h}{R_{earth} + h}`
        """
        R_earth = Phys.R_earth
        H = R_earth*h/(R_earth + h)
        return H

    @staticmethod
    def h_from_H(H):
        """Convert geopotential to geometric altitude.

        :math:`h = \\frac{R_{earth} H}{R_{earth} - H}`
        """
        R_earth = Phys.R_earth
        h = R_earth*H/(R_earth - H)
        return h

    @property
    @to_base_units_wrapper
    @_len1array_to_scalar
    def h(self):
        """Geometric height :math:`h`"""
        return self._h

    @property
    @to_base_units_wrapper
    @_len1array_to_scalar
    def H(self):
        """Geopotential height :math:`H`"""
        return self._H

    @property
    @to_base_units_wrapper
    @_len1array_to_scalar
    def p(self):
        """Air pressure :math:`p`"""
        H_base = _atleast_1d(self.layer.H_base)
        T_base = _atleast_1d(self.layer.T_base)
        T_grad = _atleast_1d(self.layer.T_grad)
        p_base = _atleast_1d(self.layer.p_base)

        H = _atleast_1d(self.H)
        T = _atleast_1d(self.T)
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

    @property
    @to_base_units_wrapper
    @_len1array_to_scalar
    def T(self):
        """Ambient air temperature :math:`T`"""
        T_grad = _atleast_1d(self.layer.T_grad)
        H_base = _atleast_1d(self.layer.H_base)
        T_base = _atleast_1d(self.layer.T_base)
        T = T_base + T_grad*(self.H - H_base)
        return T

    @property
    @to_base_units_wrapper
    @_len1array_to_scalar
    def rho(self):
        """Ambient air density :math:`\\rho`"""
        p = self.p
        T = self.T
        R_air = Phys.R_air
        rho = p/(R_air*T)
        return rho

    @property
    @to_base_units_wrapper
    @_len1array_to_scalar
    def a(self):
        """Ambient speed of sound :math:`a`"""
        T = self.T
        gamma_air = Phys.gamma_air
        R_air = Phys.R_air
        a = sqrt(gamma_air*R_air*T)
        return a

    @property
    @to_base_units_wrapper
    @_len1array_to_scalar
    def mu(self):
        """Ambient dynamic viscosity :math:`\\mu`"""
        T = self.T
        beta_s = Atmo.beta_s
        S = Atmo.S
        mu = beta_s*T**1.5/(T + S)
        return mu

    @property
    @to_base_units_wrapper
    @_len1array_to_scalar
    def nu(self):
        """Ambient kinematic viscosity :math:`\\nu`"""
        mu = self.mu
        rho = self.rho
        nu = mu/rho
        return nu

    @property
    @to_base_units_wrapper
    @_len1array_to_scalar
    def k(self):
        """Ambient thermal conductivity :math:`k`"""
        # Empirical formula requires T in Kelvin
        T_K = self.T.to('K').magnitude
        c1 = 0.002648151
        k_ = c1*T_K**1.5 / (T_K + (245.4*10**(-12/T_K)))

        # Re-assign to dimensioned variable
        k = k_ * unit('W/m/K')
        return k

    @property
    @to_base_units_wrapper
    @_len1array_to_scalar
    def g(self):
        """Gravitational acceleration at altitude :math:`g`"""
        h = self.h
        g_0 = Phys.g
        R_earth = Phys.R_earth
        g = g_0*(R_earth/(R_earth + h))**2
        return g

    @property
    @to_base_units_wrapper
    @_len1array_to_scalar
    def mean_free_path(self):
        """Mean free path"""
        p = self.p
        T = self.T
        N_A = Phys.N_A
        R = Phys.R
        sigma = Phys.collision_diam_air

        n = N_A*p/(R*T)  # number density
        mean_free_path = 1/(sqrt(2)*pi*sigma**2*n)
        return mean_free_path
