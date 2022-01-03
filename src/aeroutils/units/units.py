#!/usr/bin/env python
"""Provide dimensional units for mathematical operations.
Using the pint package.

Dependencies: pint

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

from functools import wraps

from pint import UnitRegistry

unit = UnitRegistry(system='mks')
unit.default_format = '~P'
dimless = unit('dimensionless')


def check_dimensioned(inp):
    """Check that input is dimensional (type 'Quantity' from pint package) and
    *not* from a different unit registry.  """

    if isinstance(inp, unit.Quantity):
        dummy = 1 * unit('m')
        # Registry equality check
        if inp._REGISTRY is not dummy._REGISTRY:
            msg = ("Cannot operate with {} and {} of different registries."
                   "For units, try:\n\tfrom aeroutils import unit ")
            raise ValueError(
                msg.format(inp.__class__.__name__, dummy.__class__.__name__)
            )
    else:
        raise TypeError("Input value is not correctly typed! Use"
                        " dimensional type 'Quantity' from pint package.")


def check_length_dimensioned(inp):
    """Check that input is dimensional type Quantity from pint package."""
    length_dimensionality = (1*unit('ft')).dimensionality
    if not (inp.dimensionality == length_dimensionality):
        raise TypeError("Input value is not correctly typed! Use length"
                        " dimensional unit.")


def to_base_units_wrapper(func):
    """Function decorator to convert output variable units to base units. """
    @wraps(func)
    def wrapper(*args, **kwargs):
        output = func(*args, **kwargs)
        output.ito_base_units()
        return output
    return wrapper
