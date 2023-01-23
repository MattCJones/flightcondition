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
            float: pressure ratio
        """
        T0_by_T = __class__.T0_by_T(M, y)
        p0_by_p = T0_by_T**(y/(y-1))
        return p0_by_p
    @staticmethod
    def p_by_p0(M, y=Phys.gamma_air):
        """Compute ratio of static pressure to stagnation pressure (brought to
        rest isentropically).

        Args:
            M (float): Mach number
            y (float): ratio of specific heats

        Returns:
            float: pressure ratio
        """
        return 1/__class__.p0_by_p(M, y=y)

    @staticmethod
    def T0_by_T(M, y=Phys.gamma_air, r=1):
        """Compute ratio of stagnation temperature (brought to rest
        isentropically) to static temperature.

        Args:
            M (float): Mach number
            y (float): ratio of specific heats
            r (float): recovery factor (if computing recovery temperature)

        Returns:
            float: temperature ratio
        """
        T0_by_T = 1 + r*((y-1)/2)*M**2
        return T0_by_T

    @staticmethod
    def T_by_T0(M, y=Phys.gamma_air, r=1):
        """Compute ratio of static temperature to stagnation temperature
        (brought to rest isentropically).

        Args:
            M (float): Mach number
            y (float): ratio of specific heats
            r (float): recovery factor (if computing recovery temperature)

        Returns:
            float: temperature ratio
        """
        return 1/__class__.T0_by_T(M, y=y, r=r)

    @staticmethod
    def a0_by_a(M, y=Phys.gamma_air):
        """Compute ratio of stagnation sound speed (brought to rest
        isentropically) to ambient sound speed.

        Args:
            M (float): Mach number

        Returns:
            float: sound speed ratio
        """
        a0_by_a = np.sqrt(__class__.T0_by_T(M, y))
        return a0_by_a
    @staticmethod
    def a_by_a0(M, y=Phys.gamma_air):
        return 1/__class__.a0_by_a(M, y=y)

    @staticmethod
    def rho0_by_rho(M, y=Phys.gamma_air):
        """Compute ratio of stagnation density (brought to rest
        isentropically) to ambient density.

        Args:
            M (float): Mach number
            y (float): ratio of specific heats

        Returns:
            float: density ratio
        """
        T0_by_T = __class__.T0_by_T(M, y)
        rho0_by_rho = T0_by_T**(1/(y-1))
        return rho0_by_rho

    @staticmethod
    def rho_by_rho0(M, y=Phys.gamma_air):
        """Compute ratio of ambient density to stagnation density (brought to
        rest isentropically).

        Args:
            M (float): Mach number
            y (float): ratio of specific heats

        Returns:
            float: density ratio
        """
        return 1/__class__.rho0_by_rho(M, y=y)

    @staticmethod
    def A_star_by_A(M, y=Phys.gamma_air):
        """Compute ratio of local area to area at which Mach equals 1.

        Args:
            M (float): Mach number

        Returns:
            float: area ratio
        """
        T0_by_T = __class__.T0_by_T(M, y)
        A_star_by_A = (
            (1/M**2) * ( (2/(y+1))*T0_by_T )**((y+1)/(y-1))
        )**(-1/2)
        return A_star_by_A

    @staticmethod
    def A_by_A_star(M, y=Phys.gamma_air):
        """Compute ratio of area at which Mach equals 1 to local area.

        Args:
            M (float): Mach number

        Returns:
            float: area ratio
        """
        return 1/__class__.A_star_by_A(M, y=y)

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
        M = __class__.M_from_T0_by_T(T0_by_T, y=y)
        return M

    @staticmethod
    def M_from_p_by_p0(p_by_p0, y=Phys.gamma_air):
        """Compute Mach number from ratio of static pressure to stagnation
        pressure (brought to rest isentropically).

        Args:
            p_by_p0 (float): ratio of pressure to stagnation pressure
            y (float): ratio of specific heats

        Returns:
            float: Mach number
        """
        return __class__.M_from_p0_by_p(1/p_by_p0, y=y)

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

    @staticmethod
    def M_from_T_by_T0(T_by_T0, y=Phys.gamma_air):
        """Compute Mach number from ratio of static temperature to stagnation
        temperature (brought to rest isentropically).

        Args:
            T_by_T0 (float): ratio of temperature to stagnation temperature
            y (float): ratio of specific heats

        Returns:
            float: Mach number
        """
        return __class__.M_from_T0_by_T(1/T_by_T0, y=y)

    @staticmethod
    def M_from_a0_by_a(a0_by_a, y=Phys.gamma_air):
        """Compute Mach number from ratio of stagnation sound speed (brought to
        rest isentropically) to ambient sound speed.

        Args:
            a0_by_a (float): ratio of stagnation sound speed to sound speed
            y (float): ratio of specific heats

        Returns:
            float: Mach number
        """
        T0_by_T = a0_by_a**2
        M = __class__.M_from_T0_by_T(T0_by_T, y=y)
        return M

    @staticmethod
    def M_from_a_by_a0(a_by_a0, y=Phys.gamma_air):
        """Compute Mach number from ratio of ambient sound speed to stagnation
        sound speed (brought to rest isentropically).

        Args:
            a_by_a0 (float): ratio of sound speed to stagnation sound speed
            y (float): ratio of specific heats

        Returns:
            float: Mach number
        """
        return __class__.M_from_a0_by_a(1/a_by_a0, y=y)

    @staticmethod
    def M_from_rho0_by_rho(rho0_by_rho, y=Phys.gamma_air):
        """Compute Mach number from ratio of stagnation density (brought to
        rest isentropically) to ambient density.

        Args:
            rho0_by_rho (float): ratio of stagnation density to density
            y (float): ratio of specific heats

        Returns:
            float: Mach number
        """
        T0_by_T = rho0_by_rho**(y-1)
        M = __class__.M_from_T0_by_T(T0_by_T, y=y)
        return M

    @staticmethod
    def M_from_rho_by_rho0(rho_by_rho0, y=Phys.gamma_air):
        """Compute Mach number from ratio of ambient density to stagnation
        density (brought to rest isentropically).

        Args:
            rho_by_rho0 (float): ratio of density to stagnation density
            y (float): ratio of specific heats

        Returns:
            float: Mach number
        """
        return __class__.M_from_rho0_by_rho(1/rho_by_rho0, y=y)
