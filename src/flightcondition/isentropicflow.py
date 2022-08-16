#!/usr/bin/env python
"""Compute isentropic flow relations.

Dependencies: numpy, pint

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

# flake8: noqa E201

import numpy as np

from flightcondition.constants import PhysicalConstants as Phys


class IsentropicFlow:
    """Compute isentropic flow relations. """

    # Critical values
    a_star_by_a0 = 0.913
    p_star_by_p0 = 0.528
    rho_star_by_rho0 = 0.634

    @staticmethod
    def p0_by_p(M, y=Phys.gamma_air):
        """Compute ratio of stagnation pressure (brought to rest
        isentropically) to static pressure.

        Args:
            M (float): Mach number
            y (float): ratio of specific heats

        Returns:
            float: ratio of stagnation pressure to static pressure
        """
        T0_by_T = __class__.T0_by_T(M, y)
        p0_by_p = T0_by_T**(y/(y-1))
        return p0_by_p

    @staticmethod
    def T0_by_T(M, y=Phys.gamma_air, r=1):
        """Compute ratio of stagnation temperature (brought to rest
        isentropically) to static temperature.

        Args:
            M (float): Mach number
            y (float): ratio of specific heats
            r (float): recovery factor (if computing recovery temperature)

        Returns:
            float: ratio of stagnation temperature to static temperature
        """
        T0_by_T = 1 + r*((y-1)/2)*M**2
        return T0_by_T

    @staticmethod
    def a0_by_a(M, y=Phys.gamma_air):
        """Compute ratio of stagnation sound speed (brought to rest
        isentropically) to ambient sound speed.

        Args:
            M (float): Mach number

        Returns:
            float: ratio of stagnation sound speed to static sound speed
        """
        a0_by_a = np.sqrt(__class__.T0_by_T(M, y))
        return a0_by_a

    @staticmethod
    def rho0_by_rho(M, y=Phys.gamma_air):
        """Compute ratio of stagnation density (brought to rest
        isentropically) to ambient density.

        Args:
            M (float): Mach number
            y (float): ratio of specific heats

        Returns:
            float: ratio of stagnation density to static density
        """
        T0_by_T = __class__.T0_by_T(M, y)
        rho0_by_rho = T0_by_T**(1/(y-1))
        return rho0_by_rho

    @staticmethod
    def A_star_by_A(M, y=Phys.gamma_air):
        """Compute ratio of area to area at which Mach equals 1.

        Args:
            M (float): Mach number

        Returns:
            float: ratio of area to area at which Mach equals 1
        """
        T0_by_T = __class__.T0_by_T(M, y)
        A_star_by_A = (
            (1/M**2) * ( (2/(y+1))*T0_by_T )**((y+1)/(y-1))
        )**(-1/2)
        return A_star_by_A

    @staticmethod
    def M_from_p0_by_p(p0_by_p, y=Phys.gamma_air):
        """Compute Mach number from ratio of stagnation pressure (brought to
        rest isentropically) to static pressure.

        Args:
            p0_by_p (float): ratio of stagnation pressure to pressure
            y (float): ratio of specific heats

        Returns:
            float: Mach number
        """
        T0_by_T = p0_by_p**((y-1)/y)
        M = __class__.M_from_T0_by_T(T0_by_T)
        return M

    @staticmethod
    def M_from_T0_by_T(T0_by_T, y=Phys.gamma_air):
        """Compute Mach number from ratio of stagnation temperature (brought to
        rest isentropically) to static temperature.

        Args:
            T0_by_T (float): ratio of stagnation temperature to temperature
            y (float): ratio of specific heats

        Returns:
            float: Mach number
        """
        # Derive from: T0_by_T = 1 + ((y-1)/2)*M**2
        M = np.sqrt( (T0_by_T - 1)/((y-1)/2) )
        return M
