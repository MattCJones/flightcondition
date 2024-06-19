#!/usr/bin/env python
"""Dimensioned constants.

Dependencies: numpy, pint

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

from numpy import array

from flightcondition.units import unit, dimless


class PhysicalConstants():
    """Dimensionalized constants for general physics.
    Air properties assume a calorically perfect gas in non-extreme conditions
    (approximately: T < 400 K and p < 20 psi).

    Attributes:
        g (acceleration): Acceleration due to gravity
        R (energy/temperature/mol): Universal gas constant
        R_air (energy/temperature/mass): Gas constant for air (non-extreme
            conditions)
        Cp_air (energy/temperature/mass): Specific heat at constant pressure
            for air (non-extreme conditions)
        Cv_air (energy/temperature/mass): Specific heat at constant volume
            for air (non-extreme conditions)
        y_air (dimless): Ratio of specific heats for air (low temperature
            and pressure conditions)
        collision_diam_air (length): Effective collision diameter of an air
            molecule
        N_A (1/mol): Avogadro's constant
        R_earth (length): Radius of Earth
    """
    g = 9.80665 * unit('m/s^2')
    R = 8.31432 * unit('N m / (mol K)')
    R_air = 287.05287 * unit('J/(K kg)')
    Cp_air = 1.005e3 * unit('J/(K kg)')
    Cv_air = 0.718e3 * unit('J/(K kg)')
    y_air = 1.4 * dimless
    Pr_air = 0.7 * dimless
    collision_diam_air = 0.365e-9 * unit('m')
    N_A = 6.02257e23 * unit('1/mol')
    k_B = R/N_A
    R_earth = 6.356766e3 * unit('km')  # 6378.1370 km eq, 6356.7523 km polar
    G = 6.67408e-11 * unit("m^3/kg/s^2")
    M_earth = 5.9722e24 * unit("kg")


class AtmosphereConstants():
    """Dimensionalized constants for atmospheric modeling.

    Attributes:
        p_0 (pressure): Sea level ambient pressure
        T_ice (temperature): Sea level ice point temperature
        T_0 (temperature): Sea level ambient temperature
        rho_0 (density): Sea level ambient density
        S (temperature): Sutherland's empirical constant in equation for
            dynamic viscosity
        beta_s (mass/length/time/temperature^0.5): Sutherland's empirical
            constant in equation for dynamic viscosity
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
        "Thermosphere",
        "Exosphere",
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
        84.8525,
        630.5631,  # h = 700 km
        3886.3221,  # h = 10000 km
    ]) * unit('km')

    # Layer temperature gradients
    T_grad = array([
        -6.5,   # H = -5 km
        -6.5,   # H =  0 km
        0.0,    # H = 11 km
        1.0,    # H = 20 km
        2.8,    # H = 32 km
        0.0,    # H = 47 km
        -2.8,   # H = 51 km
        -2.0,   # H = 71 km
        -2.0,   # h = 86 km
        -2.0,  # unused
        -2.0,  # unused
    ]) * 1e-3 * unit('K/m')

    # Layer base temperatures
    T_base = array([
        320.65,  # H = -5 km
        288.15,  # H =  0 km
        216.65,  # H = 11 km
        216.65,  # H = 20 km
        228.65,  # H = 32 km
        270.65,  # H = 47 km
        270.65,  # H = 51 km
        214.65,  # H = 71 km
        196.65,  # h = 86 km
        196.65,  # unused
        196.65,  # unused
    ]) * unit('K')

    # Layer base pressures
    p_base = array([
        1.77687e+5,  # H = -5 km
        1.01325e+5,  # H =  0 km
        2.26320e+4,  # H = 11 km
        5.47487e+3,  # H = 20 km
        8.68014e+2,  # H = 32 km
        1.10906e+2,  # H = 47 km
        6.69384e+1,  # H = 51 km
        3.95639e+0,  # H = 71 km
        8.86272e-1,  # h = 86 km
        8.86272e-1,  # unused
        8.86272e-1,  # unused
    ]) * unit('Pa')
