#!/usr/bin/env python
"""Compute isentropic flow relations.

Dependencies: numpy, pint

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

# flake8: noqa E201

import warnings

import numpy as np

from flightcondition.constants import PhysicalConstants as Phys


class IsentropicFlow:
    """Compute isentropic flow relations. """

    # Critical values
    a_star_by_a0 = 0.913
    p_star_by_p0 = 0.528
    rho_star_by_rho0 = 0.634

    @staticmethod
    def p0_by_p(M, y=Phys.gamma_air.magnitude):
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
    def p_by_p0(M, y=Phys.gamma_air.magnitude):
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
    def T0_by_T(M, y=Phys.gamma_air.magnitude, r=1):
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
    def T_by_T0(M, y=Phys.gamma_air.magnitude, r=1):
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
    def a0_by_a(M, y=Phys.gamma_air.magnitude):
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
    def a_by_a0(M, y=Phys.gamma_air.magnitude):
        return 1/__class__.a0_by_a(M, y=y)

    @staticmethod
    def rho0_by_rho(M, y=Phys.gamma_air.magnitude):
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
    def rho_by_rho0(M, y=Phys.gamma_air.magnitude):
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
    def A_star_by_A(M, y=Phys.gamma_air.magnitude):
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
    def A_by_A_star(M, y=Phys.gamma_air.magnitude):
        """Compute ratio of area at which Mach equals 1 to local area.

        Args:
            M (float): Mach number

        Returns:
            float: area ratio
        """
        return 1/__class__.A_star_by_A(M, y=y)

    @staticmethod
    def M_from_p0_by_p(p0_by_p, y=Phys.gamma_air.magnitude):
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
    def M_from_p_by_p0(p_by_p0, y=Phys.gamma_air.magnitude):
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
    def M_from_T0_by_T(T0_by_T, y=Phys.gamma_air.magnitude):
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
    def M_from_T_by_T0(T_by_T0, y=Phys.gamma_air.magnitude):
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
    def M_from_a0_by_a(a0_by_a, y=Phys.gamma_air.magnitude):
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
    def M_from_a_by_a0(a_by_a0, y=Phys.gamma_air.magnitude):
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
    def M_from_rho0_by_rho(rho0_by_rho, y=Phys.gamma_air.magnitude):
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
    def M_from_rho_by_rho0(rho_by_rho0, y=Phys.gamma_air.magnitude):
        """Compute Mach number from ratio of ambient density to stagnation
        density (brought to rest isentropically).

        Args:
            rho_by_rho0 (float): ratio of density to stagnation density
            y (float): ratio of specific heats

        Returns:
            float: Mach number
        """
        return __class__.M_from_rho0_by_rho(1/rho_by_rho0, y=y)

    @staticmethod
    def _M_from_A_by_A_star(A_by_A_star, y=Phys.gamma_air.magnitude,
                            M_guess=5.0, tol=1e-10, max_iter=1e4):
        """Solve for subsonic or supersonic Mach number from ratio of local
        area to choked flow area.

        Args:
            A_by_A_star (float): area ratio
            y (float): ratio of specific heats
            M_guess (float): first guess when solving for Mach number
            tol (float): converge error below this tolerance and exit
            max_iter (int): maximum number of iterations

        Returns:
            float: Mach number
        """
        def zero_func(M):
            A_by_A_star_check = __class__.A_by_A_star(M=M, y=y)
            return abs(A_by_A_star - A_by_A_star_check)

        if A_by_A_star < 1:
            M = np.nan
        else:
            M = newton_fsolve(zero_func, x0=M_guess, tol=tol,
                              max_iter=max_iter)
        return M

    @staticmethod
    def subsonic_M_from_A_by_A_star(A_by_A_star, y=Phys.gamma_air.magnitude,
                                    M_guess=None, tol=1e-10, max_iter=1e4):
        """Solve for subsonic Mach number from ratio of local area to choked
        flow area.

        Args:
            A_by_A_star (float): area ratio
            y (float): ratio of specific heats
            M_guess (float): first guess when solving for Mach number
            tol (float): converge error below this tolerance and exit
            max_iter (int): maximum number of iterations

        Returns:
            float: Mach number
        """
        M_guess = 0.5/A_by_A_star if M_guess is None else M_guess
        M = __class__._M_from_A_by_A_star(A_by_A_star, y=y, M_guess=M_guess,
                                          tol=tol, max_iter=max_iter)
        return M

    @staticmethod
    def supersonic_M_from_A_by_A_star(A_by_A_star, y=Phys.gamma_air.magnitude,
                                      M_guess=None, tol=1e-10, max_iter=1e4):
        """Solve for supersonic Mach number from ratio of local area to choked
        flow area.

        Args:
            A_by_A_star (float): area ratio
            y (float): ratio of specific heats
            M_guess (float): first guess when solving for Mach number
            tol (float): converge error below this tolerance and exit
            max_iter (int): maximum number of iterations

        Returns:
            float: Mach number
        """
        M_guess = 10
        M = __class__._M_from_A_by_A_star(A_by_A_star, y=y, M_guess=M_guess,
                                          tol=tol, max_iter=max_iter)
        return M


def newton_fsolve(func, x0, x1=None, tol=1e-10, max_iter=1e4):
    """Custom implementation of Newton's Method to solve for the roots of a
    function, i.e. func(x)=0.

    This custom implementation reduces the need for external dependencies.

    Args:
        func (function): solving function f(x)=0 for x
        x0 (float): initial guess for x
        tol (float): converge error below this tolerance and exit
        max_iter (int): maximum number of iterations

    Returns: solution variable x for f(x)=0
    """
    if x1 is None:
        x1 = 1.2*x0
    xnm1 = x0
    x = x1
    f = func
    for i in range(int(max_iter)):
        fxnm1 = f(xnm1)
        fx = f(x)
        fpx = (fx-fxnm1)/(x-xnm1)
        xn = x - fx/fpx
        err = abs(xn - x)
        if err < tol:
            return xn
        xnm1 = x
        x = xn
    else:  # no break
        warnings.warn(f"Exceeded max_iter={max_iter}. Breaking.")
        return xn
