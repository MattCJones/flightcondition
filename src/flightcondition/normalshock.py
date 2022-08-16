#!/usr/bin/env python
"""Compute normal shock quantities.

Dependencies: numpy, pint

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

# flake8: noqa E201

import numpy as np

from flightcondition.constants import PhysicalConstants as Phys


class NormalShock:
    """Compute normal shock quantities. """

    @staticmethod
    def M2(M1, y=Phys.gamma_air):
        """Compute Mach number downstream of the shock.

        Args:
            M (float): Mach number upstream of shock

        Returns:
            float: Mach number downstream of shock
        """
        M2 = np.sqrt( (M1**2*(y-1) + 2) / (2*y*M1**2 - (y-1)) )
        return M2

    @staticmethod
    def p2_by_p1(M1, y=Phys.gamma_air):
        """Compute ratio of static pressure over shock, downstream/upstream.

        Args:
            M (float): Mach number upstream of shock

        Returns:
            float: ratio of static pressure over shock, downstream/upstream
        """
        p2_by_p1 = (2*y*M1**2)/(y+1) - (y-1)/(y+1)
        return p2_by_p1

    @staticmethod
    def rho2_by_rho1(M1, y=Phys.gamma_air):
        """Compute ratio of density over shock, downstream/upstream.

        Args:
            M (float): Mach number upstream of shock

        Returns:
            float: ratio of density over shock, downstream/upstream
        """
        rho2_by_rho1 = ((y+1)*M1**2) / ((y-1)*M1**2 + 2)
        return rho2_by_rho1

    @staticmethod
    def T2_by_T1(M1, y=Phys.gamma_air):
        """Compute ratio of static pressure over shock, downstream/upstream.

        Args:
            M (float): Mach number upstream of shock

        Returns:
            float: ratio of temperature over shock, downstream/upstream
        """
        T2_by_T1 = (
            (1 + ((y-1)/2)*M1**2)
            *
            (((2*y)/(y-1))*M1**2 - 1)
            /
            (M1**2*((2*y)/(y-1) + (y-1)/2))
        )
        return T2_by_T1

    @staticmethod
    def p02_by_p01(M1, y=Phys.gamma_air):
        """Compute ratio of stagnation pressure over shock, downstream/upstream.

        Args:
            M (float): Mach number upstream of shock

        Returns:
            float: ratio of temperature over shock, downstream/upstream
        """
        p02_by_p01 = (
            (
                (((y+1)/2)*M1**2)
                /
                (1 + ((y-1)/2)*M1**2)
            )**(y/(y-1))
            *
            (
                1
                /
                (((2*y)/(y+1))*M1**2 - (y-1)/(y+1))
            )**(1/(y-1))
        )
        return p02_by_p01
