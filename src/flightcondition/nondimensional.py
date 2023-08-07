#!/usr/bin/env python
"""Compute non-dimensional flow quantities.

Dependencies: numpy, pint

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

# flake8: noqa E201

import numpy as np

from flightcondition.constants import PhysicalConstants as Phys
from flightcondition.units import dimless, to_base_units_wrapper

class NonDimensional:
    """Compute non-dimensional flow quantities. """

    @staticmethod
    def knudsen_number(ell, MFP):
        """Compute the Knudsen number :math:`K_n`

        Args:
            ell (length): Length scale
            MFP (length): Mean free path

        Returns:
            dimless: Knudsen number
        """
        Kn = MFP/ell
        return Kn

    @staticmethod
    def mach_number(U, a):
        """Compute Mach number

        Args:
            U (speed): Velocity
            a (speed): Sound speed

        Returns:
            dimless: Mach number
        """
        return U/a

    @staticmethod
    def mach_velocity(M, a):
        """Compute velocity from Mach number

        Args:
            M (dimless): Mach number
            a (speed): Sound speed

        Returns:
            speed: Velocity
        """
        U = M*a if M is not None else None
        return U

    @staticmethod
    def reynolds_number(U, L, nu):
        """Compute Reynolds number :math:`Re`, the ratio of inertial to viscous
        forces.

        Args:
            U (speed): Velocity
            L (length): Length scale
            nu (length^2/time): Kinematic viscosity

        Returns:
            dimless: Reynolds number
        """
        Re_L = U*L/nu
        return Re_L

    @staticmethod
    def reynolds_number_velocity(Re_L, L, nu):
        """Compute velocity from Reynolds number, length, and kinematic
        viscosity.

        Args:
            Re_L (dimless): Reynolds number
            L (length): Length scale
            nu (length^2/time): Kinematic viscosity

        Returns:
            speed: Velocity
        """
        U = Re_L*nu/L
        return U

    @staticmethod
    def reynolds_number_length(Re_L, U, nu):
        """Compute length from Reynolds number, velocity, and kinematic
        viscosity.

        Args:
            Re_L (dimless): Reynolds number
            U (speed): Velocity
            nu (length^2/time): Kinematic viscosity

        Returns:
            speed: Velocity
        """
        L = Re_L*nu/U
        return L

    @staticmethod
    def reynolds_per_length(U, nu, length_unit='in'):
        """Compute Reynolds number divided by length unit.

        Reynolds number is,
            Re = U * L / nu

        For some fluid dynamics solvers the user must input the Reynolds number
        in terms of Reynolds number per-unit-length,
            Re_by_L = Re_L / L_in_grid_units
                      = U * (L/L_in_grid_units) / nu

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
            U (speed): Velocity
            nu (speed^2/length): Kinematic viscosity
            length_unit (length): Desired length unit as string, ('in', 'mm',
                'cm')

        Returns:
            dimless: Reynolds number
        """
        Re_by_length_unit = U/nu
        if length_unit is not None:
            Re_by_length_unit.ito(f"1/{length_unit}")
        return Re_by_length_unit
