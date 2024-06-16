#!/usr/bin/env python
"""Compute atmospheric properties from ICAO 1993 standard atmosphere model.

Dependencies: numpy, pint

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

import re
import warnings

import numpy as np

from pymsis import msis

from flightcondition.constants import PhysicalConstants as Phys
from flightcondition.constants import AtmosphereConstants as Atmo
from flightcondition.common import AliasAttributes, DimensionalData,\
    _property_decorators
from flightcondition.units import unit, check_dimensioned,\
    check_length_dimensioned, check_US_length_units,\
    check_pressure_dimensioned, check_temperature_dimensioned


class Layer(DimensionalData):
    """Class to compute and store layer data. """

    names_dict = {
        'name': 'layer_name',
        'H_base': 'base_geopotential_height',
        'T_base': 'base_geopotential_temperature',
        'T_grad': 'temperature_gradient',
        'p_base': 'base_static_pressure',
    }

    def __init__(self, H_arr, full_output=True, units=None):
        """Initialize Layer nested class.

        Args:
            H_arr (length): Geopotential altitude
            full_output (bool): Set to True for full output
            units (str): Set to 'US' for US units or 'SI' for SI
        """
        H_arr = np.atleast_1d(H_arr)
        self.full_output = full_output
        self.units = units

        self._layer_name = [""]*np.size(H_arr)
        self._H_base = np.zeros_like(H_arr) * unit('m')
        self._T_base = np.zeros_like(H_arr) * unit('K')
        self._T_grad = np.zeros_like(H_arr) * unit('K/m')
        self._p_base = np.zeros_like(H_arr) * unit('Pa')

        for idx, H in enumerate(H_arr):
            jdx = __class__._layer_idx(H)
            self._layer_name[idx] = Atmo.layer_names[jdx]
            self._H_base[idx] = Atmo.H_base[jdx]
            self._T_base[idx] = Atmo.T_base[jdx]
            self._T_grad[idx] = Atmo.T_grad[jdx]
            self._p_base[idx] = Atmo.p_base[jdx]

        # Initialize access by full quantity name through .byname.<name>
        self.byname = AliasAttributes(
            varsobj_arr=[self, ], names_dict_arr=[__class__.names_dict, ])

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

    def tostring(self, full_output=None, units=None, max_sym_chars=None,
                 max_name_chars=None, pretty_print=True):
        """Output string representation of class object.

        Args:
            full_output (bool): Set to True for full output
            units (str): Set to 'US' for US units or 'SI' for SI
            max_sym_chars (int): Maximum characters in symbol name
            max_name_chars (int): Maximum characters iin full name
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

        if units is not None:
            self.units = units

        if self.units == 'US':
            H_base_units   = 'kft'
            T_base_units   = 'degR'
            T_grad_units   = 'degR/kft'
            p_base_units   = 'lbf/ft^2'
        else:  # SI units
            H_base_units   = 'km'
            T_base_units   = 'degK'
            T_grad_units   = 'degK/km'
            p_base_units   = 'Pa'

        # Insert longer variable name into output
        if max_sym_chars is None:
            max_sym_chars = max([len(v) for v in self.names_dict.keys()])
        if max_name_chars is None:
            max_name_chars = max([len(v) for v in self.names_dict.values()])

        H_base_str = self._vartostr(var=self.H_base, var_str='H_base',
                                    to_units=H_base_units,
                                    max_sym_chars=max_sym_chars,
                                    max_name_chars=max_name_chars,
                                    fmt_val="10.5g", pretty_print=pretty_print)
        T_base_str = self._vartostr(var=self.T_base, var_str='T_base',
                                    to_units=T_base_units,
                                    max_sym_chars=max_sym_chars,
                                    max_name_chars=max_name_chars,
                                    fmt_val="10.5g", pretty_print=pretty_print)
        T_grad_str = self._vartostr(var=self.T_grad, var_str='T_grad',
                                    to_units=T_grad_units,
                                    max_sym_chars=max_sym_chars,
                                    max_name_chars=max_name_chars,
                                    fmt_val="10.5g", pretty_print=pretty_print)
        p_base_str = self._vartostr(var=self.p_base, var_str='p_base',
                                    to_units=p_base_units,
                                    max_sym_chars=max_sym_chars,
                                    max_name_chars=max_name_chars,
                                    fmt_val="10.5g", pretty_print=pretty_print)

        # Create layer string
        if type(self.name) is np.str_:  # singular array
            trunc_layer_name = self.name
        else:  # truncate layer name to no more than 10 character
            trunc_layer_name = [f"{s[:10]}" for s in self.name]
        layer_str = self._vartostr(var=trunc_layer_name, var_str='name',
                                   to_units=None, max_sym_chars=max_sym_chars,
                                   max_name_chars=max_name_chars,
                                   fmt_val="", pretty_print=False)

        # Assemble output string
        if full_output:
            repr_str = (f"{layer_str}\n{H_base_str}\n{T_base_str}\n"
                        f"{T_grad_str}\n{p_base_str}")
        else:
            repr_str = (f"{layer_str}")
        return repr_str

    @property
    def name(self):
        """Layer name """
        if np.size(self._layer_name) == 1:
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


class Species(DimensionalData):
    """Class to store species data. """
    names_dict = {
        'N2': 'diatomic_nitrogen',
        'O2': 'diatomic_oxygen',
        'O': 'atomic_oxygen',
        'He': 'helium',
        'H': 'hydrogen',
        'Ar': 'argon',
        'N': 'nitrogen',
        'aO': 'anomalous_oxygen',
        'NO': 'nitric_oxide',
        'nrho': 'number_density',
        'R_gas': 'gas_constant',
        'y': 'specific_heat_ratio',
    }

    def __init__(self, species_arr, rho_mix, full_output=True,
                 units=None):
        """Initialize Species class.

        Args:
            species_arr (array): Array of species concentrations
            rho_mix (density): Density of the mixture
        """
        nN2 = species_arr[0]
        nO2 = species_arr[1]
        nO = species_arr[2]
        nHe = species_arr[3]
        nH = species_arr[4]
        nAr = species_arr[5]
        nN = species_arr[6]
        naO = species_arr[7]
        nNO = species_arr[8]

        self.full_output = full_output
        self.units = units

        self._nrho = nN2 + nO2 + nO + nHe + nH + nAr + nN\
                   + naO + nNO

        # Compute number concentrations of each species
        self._N2 = nN2/self._nrho
        self._O2 = nO2/self._nrho
        self._O = nO/self._nrho
        self._He = nHe/self._nrho
        self._H = nH/self._nrho
        self._Ar = nAr/self._nrho
        self._N = nN/self._nrho
        self._aO = naO/self._nrho
        self._NO = nNO/self._nrho

        # Compute individual gas constant for the mixture
        N_A = Phys.N_A
        _N = 14.0067 * unit("g/mol")
        _O = 15.9994 * unit("g/mol")
        _He = 4.002602 * unit("g/mol")
        _H = 1.00794 * unit("g/mol")
        _Ar = 39.948 * unit("g/mol")
        MW_N2 = _N*2
        MW_O2 = _O*2
        MW_O = _O
        MW_He = _He
        MW_H = _H
        MW_Ar = _Ar
        MW_N = _N
        MW_aO = _O
        MW_NO = _N + _O

        # Compute species densities
        rho_N2 = nN2*MW_N2/N_A
        rho_O2 = nO2*MW_O2/N_A
        rho_O = nO*MW_O/N_A
        rho_He = nHe*MW_He/N_A
        rho_H = nH*MW_H/N_A
        rho_Ar = nAr*MW_Ar/N_A
        rho_N = nN*MW_N/N_A
        rho_aO = naO*MW_aO/N_A
        rho_NO = nNO*MW_NO/N_A

        # Compute gas constants
        R = Phys.R
        R_N2 = R/MW_N2
        R_O2 = R/MW_O2
        R_O = R/MW_O
        R_He = R/MW_He
        R_H = R/MW_H
        R_Ar = R/MW_Ar
        R_N = R/MW_N
        R_aO = R/MW_aO
        R_NO = R/MW_NO

        self._R_gas = (
                R_N2*rho_N2 + R_O2*rho_O2 + R_O*rho_O + R_He*rho_He
                + R_H*rho_H + R_Ar*rho_Ar + R_N*rho_N + R_aO*rho_aO
                + R_NO*rho_NO
                ) / rho_mix

        # Compute ratio of specific heats
        # TODO 2024-06-15: probably more complexity, e.g., temperature
        y_N2 = 1.400
        y_O2 = 1.395
        y_O = 1.395  # TODO 2024-06-15: determine this
        y_He = 1.667
        y_H = 1.405
        y_Ar = 1.667
        y_N = 1.400  # TODO 2024-06-15: determine this
        y_aO = 1.395  # TODO 2024-06-15: determine this
        y_NO = 1.386

        self._y = (
                y_N2*rho_N2 + y_O2*rho_O2 + y_O*rho_O + y_He*rho_He
                + y_H*rho_H + y_Ar*rho_Ar + y_N*rho_N + y_aO*rho_aO
                + y_NO*rho_NO
                ) / rho_mix


    def tostring(self, full_output=None, units=None, max_sym_chars=None,
                 max_name_chars=None, pretty_print=True):
        """Output string representation of class object.

        Args:
            full_output (bool): Set to True for full output
            units (str): Set to 'US' for US units or 'SI' for SI
            max_sym_chars (int): Maximum characters in symbol name
            max_name_chars (int): Maximum characters iin full name
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

        if units is not None:
            self.units = units

        if self.units == 'US':
            nrho_units   = "1/ft^3"
            R_gas_units   = "ft lbf/(slug degR)"
        else:  # SI units
            nrho_units   = "1/m^3"
            R_gas_units   = "J/(kg degK)"
        species_units   = ""

        # Insert longer variable name into output
        if max_sym_chars is None:
            max_sym_chars = max([len(v) for v in self.names_dict.keys()])
        if max_name_chars is None:
            max_name_chars = max([len(v) for v in self.names_dict.values()])

        N2_str = self._vartostr(var=self.N2, var_str='N2',
                                to_units=species_units,
                                max_sym_chars=max_sym_chars,
                                max_name_chars=max_name_chars,
                                fmt_val="10.5g", pretty_print=pretty_print)
        O2_str = self._vartostr(var=self.O2, var_str='O2',
                                to_units=species_units,
                                max_sym_chars=max_sym_chars,
                                max_name_chars=max_name_chars,
                                fmt_val="10.5g", pretty_print=pretty_print)
        O_str = self._vartostr(var=self.O, var_str='O',
                               to_units=species_units,
                               max_sym_chars=max_sym_chars,
                               max_name_chars=max_name_chars,
                               fmt_val="10.5g", pretty_print=pretty_print)
        He_str = self._vartostr(var=self.He, var_str='He',
                                to_units=species_units,
                                max_sym_chars=max_sym_chars,
                                max_name_chars=max_name_chars,
                                fmt_val="10.5g", pretty_print=pretty_print)
        H_str = self._vartostr(var=self.H, var_str='H',
                               to_units=species_units,
                               max_sym_chars=max_sym_chars,
                               max_name_chars=max_name_chars,
                               fmt_val="10.5g", pretty_print=pretty_print)
        Ar_str = self._vartostr(var=self.Ar, var_str='Ar',
                                to_units=species_units,
                                max_sym_chars=max_sym_chars,
                                max_name_chars=max_name_chars,
                                fmt_val="10.5g", pretty_print=pretty_print)
        N_str = self._vartostr(var=self.N, var_str='N',
                               to_units=species_units,
                               max_sym_chars=max_sym_chars,
                               max_name_chars=max_name_chars,
                               fmt_val="10.5g", pretty_print=pretty_print)
        aO_str = self._vartostr(var=self.aO, var_str='aO',
                                to_units=species_units,
                                max_sym_chars=max_sym_chars,
                                max_name_chars=max_name_chars,
                                fmt_val="10.5g", pretty_print=pretty_print)
        NO_str = self._vartostr(var=self.NO, var_str='NO',
                                to_units=species_units,
                                max_sym_chars=max_sym_chars,
                                max_name_chars=max_name_chars,
                                fmt_val="10.5g", pretty_print=pretty_print)

        nrho_str = self._vartostr(var=self.nrho, var_str='nrho',
                                to_units=nrho_units,
                                max_sym_chars=max_sym_chars,
                                max_name_chars=max_name_chars,
                                fmt_val="10.5g", pretty_print=pretty_print)

        R_gas_str = self._vartostr(var=self.R_gas, var_str='R_gas',
                                   to_units=R_gas_units,
                                   max_sym_chars=max_sym_chars,
                                   max_name_chars=max_name_chars,
                                   fmt_val="10.5g", pretty_print=pretty_print)

        y_str = self._vartostr(var=self.y, var_str='y', to_units="",
                               max_sym_chars=max_sym_chars,
                               max_name_chars=max_name_chars,
                               fmt_val="10.5g", pretty_print=pretty_print)

        # Assemble output string
        if full_output:
            repr_str = (f"{N2_str}\n{O2_str}\n{O_str}\n{He_str}\n{H_str}"
                        f"\n{Ar_str}\n{N_str}\n{aO_str}\n{NO_str}\n{nrho_str}"
                        #f"\n{R_gas_str}\n{y_str}"
                       )
        else:
            repr_str = (f"{N2_str}\n{O2_str}\n{O_str}\n{He_str}\n{H_str}"
                        f"\n{Ar_str}\n{N_str}\n{aO_str}\n{NO_str}\n{nrho_str}")
        return repr_str


    @_property_decorators
    def N2(self):
        """Diatomic nitrogen. """
        return self._N2

    @_property_decorators
    def O2(self):
        """Diatomic oxygen. """
        return self._O2

    @_property_decorators
    def O(self):
        """Atomic oxygen. """
        return self._O

    @_property_decorators
    def He(self):
        """Helium. """
        return self._He

    @_property_decorators
    def H(self):
        """Atomic hydrogen. """
        return self._H

    @_property_decorators
    def Ar(self):
        """Argon. """
        return self._Ar

    @_property_decorators
    def N(self):
        """Atomic nitrogen. """
        return self._N

    @_property_decorators
    def aO(self):
        """Anomalous oxygen. """
        return self._aO

    @_property_decorators
    def NO(self):
        """Nitric oxide. """
        return self._NO

    @_property_decorators
    def nrho(self):
        """Number density. """
        return self._nrho

    @_property_decorators
    def R_gas(self):
        """Individual gas constant of gas mixture. """
        return self._R_gas

    @_property_decorators
    def y(self):
        """Individual gas constant of gas mixture. """
        return self._y


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
        #print(f"\n{atm.tostring(full_output=False, units='US')}")

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

    names_dict = {
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

    def __init__(self, h=None, units=None, full_output=None, model=None,
                 datetime=None, lon=None, lat=None, f107=None,
                 f107a=None, **kwargs):
        """Input geometric altitude - object contains the corresponding
        atmospheric quantities.

        See class definition for special hidden arguments such as `h_kft` to
        input a scalar that automatically converts to the kilo-foot dimension,
        or `h_km` for kilometers.

        Args:
            h (length): Geometric altitude - aliases are 'alt', 'altitude'
            units (str): Set to 'US' for US units or 'SI' for SI
            full_output (bool): Set to True for full output
            model (str): "standard" for Standard Atmosphere, "msis0.0" for NRL
                MSIS 0.0, "msis2.0" for MSIS 2.0, or "msis2.1" for MSIS 2.0
            datetime (str): UTC datetime in format "YYYY-MM-DDTHH:SS", e.g.,
                "2024-01-01T00:00" (MSIS models only)
            lon (float): Longitude in degrees (MSIS only)
            lat (float): Latitude in degrees (MSIS only)
            f107 (float): Daily F10.7 solar flux on the previous day (MSIS
                only)
            f107a (float): F10.7 solar flux 81-day average centered on the
                input datetime (MSIS only)
            ap (float): Daily magnetic index
        """
        self.full_output = full_output
        self._datetime = "2003-01-01T00:00" if datetime is None else datetime
        self._lon = 0.0 if lon is None else lon
        self._lat = 45.0 if lat is None else lat
        self._f107 = 150.0 if f107 is None else f107
        self._f107a = 150.0 if f107a is None else f107a
        self._species = None
        # Check units and set initially
        if units in dir(unit.sys):
            self.units = units
        else:
            self.units = 'SI'

        # Process altitude input
        # Check for hidden aliases
        h_aliases = ['alt', 'altitude', 'z']
        if h is None:
            h = __class__._arg_from_alias(h_aliases, kwargs)

        # Check if special h_kft syntactic sugar is used
        h_kft_aliases = ['h_kft', 'z_kft', 'kft']
        if h is None:
            h_kft = __class__._arg_from_alias(h_kft_aliases, kwargs)
            if h_kft is not None:
                h = h_kft * unit('kft')

        # Check if special h_km syntactic sugar is used
        h_km_aliases = ['h_km', 'z_km', 'km']
        if h is None:
            h_km = __class__._arg_from_alias(h_km_aliases, kwargs)
            if h_km is not None:
                h = h_km * unit('km')

        # Compute altitude bounds
        H_max_standard = Atmo.H_base[8]
        H_max_msis = Atmo.H_base[-1]
        h_max_msis = __class__._h_from_H(H_max_msis)
        h_max_standard = __class__._h_from_H(H_max_standard)

        # Default to 0 kft
        if h is None:
            h = 0 * unit('ft')

        # Set model type
        if model is None:
            self.model = "msis" if (h > h_max_standard).any() else "standard"
        else:
            self.model = model

        #  Set altitude bounds and set altitude
        self._H_min = Atmo.H_base[0]
        self._h_min = __class__._h_from_H(self._H_min)
        if "msis" in self.model:
            self._H_max = H_max_msis
            self._h_max = h_max_msis
        else:
            self._H_max = H_max_standard
            self._h_max = h_max_standard
        self.h = h

        # Further process unit system
        if units not in dir(unit.sys):  # check if usable system
            if check_US_length_units(h):
                self.units = 'US'
            else:
                self.units = 'SI'

        # Initialize access by full quantity name through .byname.<name>
        self.byname = AliasAttributes(
            varsobj_arr=[self, ], names_dict_arr=[__class__.names_dict, ])

    def tostring(self, full_output=None, units=None, max_sym_chars=None,
                 max_name_chars=None, pretty_print=True):
        """String representation of data structure.

        Args:
            full_output (bool): Set to True for full output
            units (str): Set to 'US' for US units or 'SI' for SI
            max_sym_chars (int): Maximum characters in symbol name
            max_name_chars (int): Maximum characters iin full name
            pretty_print (bool): Pretty print output

        Returns:
            str: String representation of class object
        """
        if units is not None:
            self.units = units

        if self.units == 'US':
            h_units   = 'kft'
            H_units   = 'kft'
            p_units   = 'lbf/ft^2'
            T_units   = 'degR'
            rho_units = 'slug/ft^3'
            a_units   = 'ft/s'
            mu_units  = 'lbf/ft^2 s'
            nu_units  = 'ft^2/s'
            k_units   = 'slug ft/s^3/degR'
            g_units   = 'ft/s^2'
            MFP_units = 'ft'
        else:  # SI units
            h_units   = 'km'
            H_units   = 'km'
            p_units   = 'Pa'
            T_units   = 'degK'
            rho_units = 'kg/m^3'
            a_units   = 'm/s'
            mu_units  = 'Pa s'
            nu_units  = 'm^2/s'
            k_units   = 'W/m/K'
            g_units   = 'm/s^2'
            MFP_units = 'm'

        # Insert longer variable name into output
        if max_sym_chars is None:
            max_sym_chars = max([len(v) for v in self.names_dict.keys()])
        if max_name_chars is None:
            max_name_chars = max([len(v) for v in self.names_dict.values()])

        h_str = self._vartostr(var=self.h, var_str='h', to_units=h_units,
                               max_sym_chars=max_sym_chars,
                               max_name_chars=max_name_chars,
                               fmt_val="10.5g", pretty_print=pretty_print)
        H_str = self._vartostr(var=self.H, var_str='H', to_units=H_units,
                               max_sym_chars=max_sym_chars,
                               max_name_chars=max_name_chars,
                               fmt_val="10.5g", pretty_print=pretty_print)
        p_str = self._vartostr(var=self.p, var_str='p', to_units=p_units,
                               max_sym_chars=max_sym_chars,
                               max_name_chars=max_name_chars,
                               fmt_val="10.5g", pretty_print=pretty_print)
        T_str = self._vartostr(var=self.T, var_str='T', to_units=T_units,
                               max_sym_chars=max_sym_chars,
                               max_name_chars=max_name_chars,
                               fmt_val="10.5g", pretty_print=pretty_print)
        fmt_val_rho = "10.4e" if units == 'US' else "10.5g"
        rho_str = self._vartostr(var=self.rho, var_str='rho',
                                 to_units=rho_units,
                                 max_sym_chars=max_sym_chars,
                                 max_name_chars=max_name_chars,
                                 fmt_val=fmt_val_rho,
                                 pretty_print=pretty_print)
        a_str = self._vartostr(var=self.a, var_str='a', to_units=a_units,
                               max_sym_chars=max_sym_chars,
                               max_name_chars=max_name_chars,
                               fmt_val="10.5g", pretty_print=pretty_print)
        mu_str = self._vartostr(var=self.mu, var_str='mu', to_units=mu_units,
                                max_sym_chars=max_sym_chars,
                                max_name_chars=max_name_chars,
                                fmt_val="10.4e", pretty_print=pretty_print)
        nu_str = self._vartostr(var=self.nu, var_str='nu', to_units=nu_units,
                                max_sym_chars=max_sym_chars,
                                max_name_chars=max_name_chars,
                                fmt_val="10.4e", pretty_print=pretty_print)
        k_str = self._vartostr(var=self.k, var_str='k', to_units=k_units,
                               max_sym_chars=max_sym_chars,
                               max_name_chars=max_name_chars,
                               fmt_val="10.4e", pretty_print=pretty_print)
        g_str = self._vartostr(var=self.g, var_str='g', to_units=g_units,
                               max_sym_chars=max_sym_chars,
                               max_name_chars=max_name_chars,
                               fmt_val="10.5g", pretty_print=pretty_print)
        MFP_str = self._vartostr(var=self.MFP, var_str='MFP',
                                 to_units=MFP_units,
                                 max_sym_chars=max_sym_chars,
                                 max_name_chars=max_name_chars,
                                 fmt_val="10.4e", pretty_print=pretty_print)

        # Create layer string
        layer_str = self.layer.tostring(full_output=False,
                                        max_sym_chars=max_sym_chars,
                                        max_name_chars=max_name_chars)

        # Determine full output flag
        if full_output is None:
            if self.full_output is None:
                full_output = True
            else:
                full_output = self.full_output

        # Assemble output string
        if full_output:
            repr_str = (f"{h_str}\n{H_str}\n{p_str}\n{T_str}\n{rho_str}\n"
                        f"{a_str}\n{mu_str}\n{nu_str}\n{k_str}\n{g_str}\n"
                        f"{MFP_str}\n{layer_str}")
        else:
            repr_str = (f"{h_str}\n{p_str}\n{T_str}\n{rho_str}\n{a_str}\n"
                        f"{nu_str}")
        return repr_str

    def _run_msis(self):
        """Run NRL MSIS model.

        Returns:
            pressure: Pressure
            temperature: Temperature
        """
        h_km = self._h.to('km').magnitude

        # Determine version
        pattern = r"msis\s*(\d+(\.\d+)?)"
        match = re.search(pattern, str(self.model).lower())
        _default_version = 2.1
        _version = float(match.group(1)) if match else _default_version

        # Format version to int if 0 or 2 to accommodate pymsis
        _version = 0 if _version == 0.0 else _version
        _version = 2 if _version == 2.0 else _version

        self.model = f"msis{_version:1.1f}"

        _out = msis.run(dates=self._datetime, lons=self._lon, lats=self._lat,
                        alts=h_km, f107s=self._f107, f107as=self._f107a,
                        version=_version)
        _out[np.isnan(_out)] = 0.0
        N_species = 9  # N2 O2 O He H Ar N aO NO
        N_fc = len(self._h)
        out = _out.reshape((N_fc, N_species+2))

        # Unpack number densities
        rho = out[:, 0] * unit('kg/m^3')
        nN2 = out[:, 1] * unit('1/m^3')
        nO2 = out[:, 2] * unit('1/m^3')
        nO = out[:, 3] * unit('1/m^3')
        nHe = out[:, 4] * unit('1/m^3')
        nH = out[:, 5] * unit('1/m^3')
        nAr = out[:, 6] * unit('1/m^3')
        nN = out[:, 7] * unit('1/m^3')
        naO = out[:, 8] * unit('1/m^3')
        nNO = out[:, 9] * unit('1/m^3')
        T = out[:, 10] * unit('degK')

        self._species = Species([nN2, nO2, nO, nHe, nH, nAr, nN, naO, nNO],
                                 rho)

        # Find pressure from temperature and density
        R_gas = self._species.R_gas
        p =  rho*(R_gas*T)

        return p, T

    @staticmethod
    @unit.wraps(unit.m, (unit.m, unit.m))
    def _H_from_h(h, R_earth=Phys.R_earth):
        """Convert geometric to geopotential altitude.

        :math:`H = \\frac{R_{earth} h}{R_{earth} + h}`

        Args:
            h (length): Geometric altitude
            R_earth (length): Radius of the Earth

        Returns:
            length: Geopotential altitude
        """
        H = R_earth*h/(R_earth + h)
        return H

    @staticmethod
    @unit.wraps(unit.m, (unit.m, unit.m))
    def _h_from_H(H, R_earth=Phys.R_earth):
        """Convert geopotential to geometric altitude.

        :math:`h = \\frac{R_{earth} H}{R_{earth} - H}`

        Args:
            H (length): Geopotential altitude
            R_earth (length): Radius of the Earth

        Returns:
            length: Geometric altitude
        """
        h = R_earth*H/(R_earth - H)

        return h

    @staticmethod
    @unit.wraps(unit('W/m/K').units, unit.K)
    def _thermal_conductivity(T_K):
        """Compute thermal conductivity.

        Args:
            T_K (temperature): Dimensional Temperature at altitude, which is
                automatically converted to Kelvin

        Returns:
            power/length/temperature: Geometric altitude in W/m/K
        """
        # Empirical formula requires T in Kelvin
        c1 = 0.002648151
        k = c1*T_K**1.5 / (T_K + (245.4*10**(-12/T_K)))
        return k

    @property
    def units(self):
        """Get unit system to use: 'SI', 'US', etc.  Available unit systems
        given by dir(unit.sys).

        Returns:
            str: Unit system
        """
        return self._units

    @units.setter
    def units(self, units):
        """Set unit system to use: 'SI', 'US', etc.  Available unit systems
        given by dir(unit.sys).

        Args:
            units (str): Unit system
        """
        if units not in dir(unit.sys):
            warnings.warn(f"'{units} is not available. Try one of the "
                          f"following: {dir(unit.sys)}")
            return
        else:
            self._units = units
            unit.default_system = units
            if hasattr(self, 'layer'):
                self.layer.units = units

    @property
    def full_output(self):
        """Enable or disable full output of data by default.

        Returns:
            bool: Full output flag
        """
        return self._full_output

    @full_output.setter
    def full_output(self, full_output):
        """Unit system to use: 'SI', 'US', etc.  Available unit systems given
        by dir(unit.sys).

        Args:
            full_output (bool): Full output flag
        """
        self._full_output = full_output

    @_property_decorators
    def h(self):
        """Get geometric altitude :math:`h`

        Returns:
            length: Geometric altitude
        """
        return self._h

    @h.setter
    def h(self, h):
        """Set geometric altitude :math:`h`

        Check that input is of type Quantity from pint package. Check that
        input is length dimension.  Check bounds.  Format as array even if
        scalar input.

        Args:
            h (length): Input scalar or array of altitudes
        """
        tofloat = 1.0
        h = np.atleast_1d(h) * tofloat

        check_dimensioned(h)
        h = h.magnitude * unit(str(h.units))  # force local unit registry
        check_length_dimensioned(h)

        if len(np.shape(h)) > 1:
            raise TypeError("Input must be scalar or 1-D array.")

        if (h < self._h_min).any() or (self._h_max < h).any():
            raise ValueError(
                f"Input altitude is out of bounds "
                f"({self._h_min:.5g} < h < {self._h_max:.5g})"
            )

        # Update quantities
        self._h = h
        self._H = self._H_from_h(self.h)
        self.layer = Layer(self.H, full_output=self.full_output,
                           units=self.units)
        self._T = None
        self._p = None
        if not self.model == "standard":
            p, T = self._run_msis()
            self._p = p
            self._T = T

    @_property_decorators
    def H(self):
        """Geopotential altitude :math:`H` """
        return self._H

    @_property_decorators
    def p(self):
        """Air pressure :math:`p` """
        # Only compute p if not set by user
        if self._p is None:
            if self.model == "standard":
                H_base = np.atleast_1d(self.layer.H_base)
                T_base = np.atleast_1d(self.layer.T_base)
                T_grad = np.atleast_1d(self.layer.T_grad)
                p_base = np.atleast_1d(self.layer.p_base)

                H = np.atleast_1d(self.H)
                T = np.atleast_1d(self.T)
                g_0 = Phys.g
                R_gas = Phys.R_air

                p = np.zeros_like(H) * unit('Pa')

                # Pressure equation changes between T_grad == 0 and T_grad != 0
                s = T_grad == 0
                p[s] = p_base[s]*np.exp((-g_0/(R_gas*T[s]))*(H[s] - H_base[s]))

                s = T_grad != 0
                p[s] = p_base[s]*(
                    1 + (T_grad[s]/T_base[s])*(H[s] - H_base[s])
                )**((1/T_grad[s])*(-g_0/R_gas))
            else:
                p, _, _ = self._run_msis()
        else:  # return user set value
            p = self._p

        return p

    @p.setter
    def p(self, p):
        """Override ambient air pressure """
        check_dimensioned(p)
        check_pressure_dimensioned(p)
        # Check that p is same size as h
        if np.size(p) != np.size(self._h):
            raise AttributeError("Input array must be same size as altitude")
        self._p = p

    @_property_decorators
    def T(self):
        """Ambient air temperature :math:`T` """
        # Only compute T if not user set
        if self._T is None:
            if self.model == "standard":
                T_grad = np.atleast_1d(self.layer.T_grad)
                H_base = np.atleast_1d(self.layer.H_base)
                T_base = np.atleast_1d(self.layer.T_base)
                T = T_base + T_grad*(self.H - H_base)
            else:
                _, T = self._run_msis()
        else:  # return user set value
            T = self._T

        return T

    @T.setter
    def T(self, T):
        """Override ambient air temperature """
        check_dimensioned(T)
        check_temperature_dimensioned(T)
        # Check that T is same size as h
        if np.size(T) != np.size(self._h):
            raise AttributeError("Input array must be same size as altitude")
        self._T = T

    @_property_decorators
    def rho(self):
        """Ambient air density :math:`\\rho` """
        p = self.p
        T = self.T
        if self.model == "standard":
            R_gas = Phys.R_air
        else:
            R_gas = self._species.R_gas
        rho = p/(R_gas*T)  # TODO 2022-08-07: separate out eqns
        return rho

    @_property_decorators
    def a(self):
        """Ambient speed of sound :math:`a` """
        T = self.T
        if self.model == "standard":
            y_air = Phys.y_air
            R_gas = Phys.R_air
        else:
            R_gas = self._species.R_gas
            y_air = self._species.y
        a = np.sqrt(y_air*R_gas*T)  # TODO 2022-08-07: separate out eqns
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
        k = __class__._thermal_conductivity(self.T)
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
        MFP = 1/(np.sqrt(2)*np.pi*sigma**2*n)
        return MFP

    @property
    def species(self):
        """Species concentrations. """
        return self._species
