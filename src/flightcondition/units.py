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
from inspect import currentframe

from pint import UnitRegistry

__all__ = ['unit', 'dimless', 'printv']

unit = UnitRegistry(system='SI')
unit.default_format = '~P'
dimless = unit('dimensionless')

_US_length_units = ('ft', 'feet', 'foot', 'kft', 'kilofoot', 'kilofeet'
                          'fts', 'feets', 'foots', 'kfts', 'kilofoots',
                          'kilofeets', 'mi', 'mile', 'miles')

def check_dimensioned(inp):
    """Check that input is dimensional (type 'Quantity' from pint package) and
    *not* from a different unit registry.  """

    if isinstance(inp, unit.Quantity):
        dummy = 1 * unit('m')
        # Registry equality check
        if inp._REGISTRY is not dummy._REGISTRY:
            msg = ("Cannot operate with {} and {} of different registries."
                   "For units, try:\n\tfrom flightcondition import unit ")
            raise ValueError(
                msg.format(inp.__class__.__name__, dummy.__class__.__name__)
            )
    else:
        raise TypeError("Input value is not correctly typed! Use"
                        " dimensional type 'Quantity' from pint package.")


def check_length_dimensioned(inp):
    """Check that input is length dimension type Quantity from pint package."""
    length_dimensionality = (1*unit('ft')).dimensionality
    if not (inp.dimensionality == length_dimensionality):
        raise TypeError("Input value is not correctly typed! Use length"
                        " dimensional unit.")


def check_US_length_units(ell):
    """Check if length unit type is an US unit

    :ell: length unit of type Quantity from pint package.
    :returns: True if US unit else False

    """
    return ell.units in _US_length_units


def check_area_dimensioned(inp):
    """Check that input is area type Quantity from pint package."""
    area_dimensionality = (1*unit('ft^2')).dimensionality
    if not (inp.dimensionality == area_dimensionality):
        raise TypeError("Input value is not correctly typed! Use area"
                        " dimensional unit.")


def name_of_var(var):
    """Find name of variables using local items.

    :var: dimensioned variable
    :returns: user-coded name of variable

    """
    try:
        local_vars = currentframe().f_back.f_back.f_locals.items()
        match = [name for name, val in local_vars if val is var]
    except AttributeError:
        local_vars = currentframe().f_back.f_locals.items()
        match = [name for name, val in local_vars if val is var]
    name = match[0] if match else "unknown"
    return name


def printv(var, to=None, var_name="", *args, **kwargs):
    """Print name and value of a Pint unit-specified variable.
    For example,

        distance = 99.9 * unit('m')
        printv(distance)
        # prints "distance = 99.9 m"

    *Note*: as of Python 3.8, simply use the f-string syntax, e.g.
        x=7
        print(f"{x=}")

    :var: variable to be printed
    :to: (str), convert to another unit
    :var_name: overwrite variable name
    :*args: additional arguments
    :**kwargs: additional keyword arguments
    :returns: None

    """
    formatted_output = (var.to(to) if to is not None else var.to_base_units())
    var_name = var_name if var_name else name_of_var(var)
    print(f"{var_name} = {formatted_output:.5g~P}", *args, **kwargs)


def to_base_units_wrapper(func):
    """Function decorator to convert output variable units to base units. """
    @wraps(func)
    def wrapper(*args, **kwargs):
        output = func(*args, **kwargs)
        output.ito_base_units()
        return output
    return wrapper
