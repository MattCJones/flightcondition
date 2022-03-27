#!/usr/bin/env python
"""
Dimensioned constants.

Dependencies: numpy, pint

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

from numpy import array

from aeroutils.units import unit, dimless


class PhysicalConstants():
    """Dimensionalized constants for general physics.

    :g: acceleration due to gravity
    :R: universal gas constant
    :R_air: gas constant for air
    :gamma_air: ratio of specific heats for air
    :collision_diam_air: effective collision diameter of an air molecule
    :N_A: Avogadro's constant
    :R_earth: radius of Earth
    """
    g = 9.80665 * unit('m/s^2')
    R = 8.31432 * unit('N m / (mol K)')
    R_air = 287.05287 * unit('J/(K*kg)')
    gamma_air = 1.4 * dimless
    collision_diam_air = 0.365e-9 * unit('m')
    N_A = 6.02257e23 * unit('1/mol')
    R_earth = 6.356766e3 * unit('km')


class AtmosphereConstants():
    """Dimensionalized constants for atmospheric modeling.

    :p_0: sea level ambient pressure
    :T_ice: sea level ice point temperature
    :T_0: sea level ambient temperature
    :rho_0: sea level ambient density
    :S: Sutherland's empirical constant in equation for dynamic viscosity
    :beta_s: Sutherland's empirical constant in equation for dynamic viscosity
    """

    # Sea level
    p_0 = 101.325e3 * unit('Pa')
    T_ice = 273.15 * unit('K')
    T_0 = 288.15 * unit('K')
    rho_0 = 1.225 * unit('kg/m^3')
    S = 110.4 * unit('K')
    beta_s = 1.458e-6 * unit('kg/(m*s*K^(1/2))')

    layer_names = array([
        "Troposphere",
        "Tropopause",
        "Stratosphere",
        "Stratosphere",
        "Stratopause",
        "Mesosphere",
        "Mesosphere",
        "Mesopause",
        ])

    # Layer base geopotential heights
    H_base = array([
        -5.0,
        0.0,
        11.0,
        20.0,
        32.0,
        47.0,
        51.0,
        71.0,
        80.0
        ]) * unit('km')

    # Layer temperature gradients
    T_grad = array([
        -6.5,
        -6.5,
        0.0,
        1.0,
        2.8,
        0.0,
        -2.8,
        -2.0,
        -2.0
        ]) * 1e-3 * unit('K/m')

    # Layer base temperatures
    T_base = array([
        320.65,
        288.15,
        216.65,
        216.65,
        228.65,
        270.65,
        270.65,
        214.65,
        196.65,
        ]) * unit('K')

    # Layer base pressures
    p_base = array([
        1.77687e+5,
        1.01325e+5,
        2.26320e+4,
        5.47487e+3,
        8.68014e+2,
        1.10906e+2,
        6.69384e+1,
        3.95639e+0,
        8.86272e-1,
        ]) * unit('Pa')
