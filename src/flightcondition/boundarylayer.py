#!/usr/bin/env python
"""Compute properties related to boundary layer flow.

Dependencies: numpy, pint

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

import warnings

import numpy as np

from flightcondition.constants import PhysicalConstants as Phys
from flightcondition.isentropicflow import IsentropicFlow
from flightcondition.units import dimless


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
        with warnings.catch_warnings():  # catch warning for divide by 0
            warnings.simplefilter("ignore")
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
        # Derivation of Prandtl (1927),
        # Eq. 7.43 in Boundary Layer Analysis 2nd. Ed. by Schetz and Bowersox
        with warnings.catch_warnings():  # catch warning for divide by 0
            warnings.simplefilter("ignore")
            delta_turb = 0.375*x/Re_x**0.2
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
        with warnings.catch_warnings():  # catch warning for divide by 0
            warnings.simplefilter("ignore")
            Cf_lamr = 2*0.664*Re_x**(-0.5)
        return Cf_lamr

    @staticmethod
    def _roskam1987_skin_friction_coeff_turb(Re_x, M=None):
        """Compute integrated skin friction coefficient for turbulent flow over
        a flat plate.

        Args:
            Re_x (dimless): Reynolds number at distance along flat plate
            M (dimless): Mach number at distance along flat plate

        Returns:
            length: skin friction coefficient
        """
        # See references:
        # * Hoak, D.E, et al, USAF Stability and Control Datcom, Flight Control
        #   Division, Air Force Flight Dynamics Laboratory, WPAFB, Ohio,
        #   45433-0000, 1978
        # * Roskam, J, Airplane Design Part VI: Preliminary Calculation of
        #   Aerodynamic, Thrust and Power Characteristics, Roskam Aviation and
        #   Engineering Corporation, 1987
        #
        # Cf_turb = A*Re_x**B
        # where A and B are obtained from experiment for a given Mach number
        # Cf_turb, Re_x, and M must have the same size

        # Process data and format as array
        Re_x = np.atleast_1d(Re_x)
        Cf_turb = np.zeros_like(Re_x) * dimless
        if M is None:
            M = np.zeros_like(Re_x) * dimless
        else:
            M = np.atleast_1d(M) * dimless

        # If M is size 1, temporarily set it to the size of Re_x
        if np.size(M) == 1 and np.size(Re_x) > 1:
            M = np.ones_like(Re_x)*M

        if not (np.size(M) == np.size(Re_x)):
            raise AttributeError("Mach number and Reynolds number arrays must "
                                 "be the same size")

        # Empirical data curve coefficients
        M_data_arr = np.array([0, 0.3, 0.7, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0])
        A_data_arr = np.array([
            0.0391, 0.0399, 0.0392, 0.0376, 0.0381, 0.0371, 0.0329, 0.0286,
            0.0261,
        ])
        B_data_arr = np.array([
            -0.157, -0.159, -0.160, -0.159, -0.161, -0.164, -0.162, -0.161,
            -0.161,
        ])

        # Interpolate to find Cf_turb for the given Mach number
        for i, _ in enumerate(M_data_arr):
            if i == len(M_data_arr)-1:
                break
            # Select lower and upper values for index
            M_i = M_data_arr[i]
            M_ip1 = M_data_arr[i+1]
            A_i = A_data_arr[i]
            A_ip1 = A_data_arr[i+1]
            B_i = B_data_arr[i]
            B_ip1 = B_data_arr[i+1]

            # Compute filter for Cf_turb array
            sa = M_i <= M
            sb = M < M_ip1
            s = sa*sb

            Cf_i = A_i*Re_x[s]**B_i
            Cf_ip1 = A_ip1*Re_x[s]**B_ip1
            Cf_turb[s] = Cf_i + ((Cf_ip1 - Cf_i)/(M_ip1 - M_i))*M[s]

        return Cf_turb

    @staticmethod
    def flat_plate_skin_friction_coeff_turb(Re_x, M=None, source='roskam1987'):
        """Compute integrated skin friction coefficient for turbulent flow over
        a flat plate.

        :math:`C_f(x) = \\tau_w / 0.5 \\rho U_\\infty^2`

        :math:`\\tau_w` is the wall shear stress, :math:`\\rho` is density,
        and :math:`U_\\infty^2` is freestream velocity.

        Args:
            Re_x (dimless): Reynolds number at distance along flat plate

        Returns:
            length: skin friction coefficient
        """
        # See references:
        # * https://www.cfd-online.com/Wiki/Skin_friction_coefficient
        # * https://www.sciencedirect.com/topics/engineering/smooth-flat-plate
        with warnings.catch_warnings():  # catch warning for divide by 0
            warnings.simplefilter("ignore")
            if source == 'calibrated_powerlaw':
                # Only valid for 5e5 < Re_x < 1e7
                Cf_turb = 0.0592*Re_x**(-0.2)
            elif source == 'schlichting':
                # Only valid for Re_x < 1e9
                Cf_turb = (2*np.log10(Re_x) - 0.65)**(-2.3)
            elif source == 'prandtl1927':
                Cf_turb = 0.074*Re_x**(-0.2)
            elif source == 'prandtl_schlichting1932':
                Cf_turb = 0.455*(np.log10(Re_x))**(-2.58)
            elif source == 'schultz_grunov1940':
                Cf_turb = 0.427*(np.log10(Re_x) - 0.407)**(-2.64)
            elif source == 'ittc1957':
                Cf_turb = 0.075*(np.log10(Re_x) - 2)**(-2)
            elif source == 'granville1977':
                Cf_turb = 0.0776*(np.log10(Re_x) - 1.88)**(-2) + 60*Re_x**(-1)
            elif source == 'roskam1987':
                # 'roskam1987' is fitted to Re and Mach number
                Cf_turb = __class__._roskam1987_skin_friction_coeff_turb(
                    Re_x=Re_x, M=M)
            elif source == 'white2006':
                Cf_turb = 0.523/(np.log(0.06*Re_x))**2
            else:
                raise ValueError(
                    f"Input reference is invalid: {source}\n\tSee source code."
                )
        return Cf_turb
