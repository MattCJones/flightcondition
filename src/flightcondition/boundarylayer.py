#!/usr/bin/env python
"""Compute properties related to boundary layer flow.

Dependencies: numpy, pint

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

import numpy as np

from flightcondition.constants import PhysicalConstants as Phys
from flightcondition.isentropicflow import IsentropicFlow


class BoundaryLayer:
    """Compute properties related to boundary layer flow. """

    @staticmethod
    def recovery_temperature_laminar(M, T, y=Phys.gamma_air, Pr=Phys.Pr_air):
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
    def recovery_temperature_turbulent(M, T, y=Phys.gamma_air, Pr=Phys.Pr_air):
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

    @staticmethod
    def flat_plate_boundary_layer_lamr(x, Re_x):
        """Compute laminar boundary layer for laminar flow over a flat plate.

        Args:
            x (length): Distance along flat plate
            R_x (dimless): Reynolds number at distance along flat plate

        Returns:
            length: boundary layer thickness
        """
        delta_lamr = 5.0*x/Re_x**0.5
        return delta_lamr

    @staticmethod
    def flat_plate_boundary_layer_turb(x, Re_x):
        """Compute turbulent boundary layer for laminar flow over a flat plate.

        Args:
            x (length): Distance along flat plate
            R_x (dimless): Reynolds number at distance along flat plate

        Returns:
            length: boundary layer thickness
        """
        delta_turb = 0.37*x/Re_x**0.2
        return delta_turb

    @staticmethod
    def flat_plate_skin_friction_coeff_lamr(Re_x):
        """Compute integrated skin friction coefficient for laminar flow over a
        flat plate.  Derived from Blasius solution.

        :math:`C_f(x) = \\tau_w / 0.5 \\rho U^2`

        :math:`\\tau_w` is the wall shear stress, :math:`\\rho` is density,
        and :math:`U_\\infty^2` is freestream velocity.

        Args:
            R_x (dimless): Reynolds number at distance along flat plate

        Returns:
            length: skin friction coefficient
        """
        # Total, integrated, Cf is 2 times local Cf
        Cf_lamr = 2*0.664*Re_x**(-0.5)
        return Cf_lamr

    @staticmethod
    def flat_plate_skin_friction_coeff_turb(Re_x, source='granville1977'):
        """Compute integrated skin friction coefficient for turbulent flow over
        a flat plate.

        :math:`C_f(x) = \\tau_w / 0.5 \\rho U_\\infty^2`

        :math:`\\tau_w` is the wall shear stress, :math:`\\rho` is density,
        and :math:`U_\\infty^2` is freestream velocity.

        Args:
            R_x (dimless): Reynolds number at distance along flat plate

        Returns:
            length: skin friction coefficient
        """
        # See https://www.cfd-online.com/Wiki/Skin_friction_coefficient
        if source == 'granville1977':
            Cf_turb = 0.0776*(np.log10(Re_x) - 1.88)**(-2) + 60*Re_x**(-1)
        elif source == 'ittc1957':
            Cf_turb = 0.075*(np.log10(Re_x) - 2)**(-2)
        elif source == 'schultz_grunov1940':
            Cf_turb = 0.427*(np.log10(Re_x) - 0.407)**(-2.64)
        elif source == 'prandtl_schlichting1932':
            Cf_turb = 0.455*(np.log10(Re_x))**(-2.58)
        elif source == 'prandtl1927':
            Cf_turb = 0.074*Re_x**(-0.2)
        elif source == 'schlichting':  # only valid for Re_x < 1e9
            Cf_turb = (2*np.log10(Re_x) - 0.65)**(-2.3)
        elif source == 'powerlaw':  # only valid for 5e5 < Re_x < 1e7
            Cf_turb = 0.0592*Re_x**(-0.2)
        else:
            raise ValueError(
                f"Input reference is invalid: {source}\n\tSee source code."
            )
        return Cf_turb
