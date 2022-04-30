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
                    'fts', 'feets', 'foots', 'kfts', 'kilofoots', 'kilofeets',
                    'mi', 'mile', 'miles')


def check_dimensioned(inp):
    """Check that input is dimensional (type 'Quantity' from pint package) and
    *not* from a different unit registry.

    Args:
        inp (object): Object to assert as dimensional type

    """
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
    """Check that input is length dimension type Quantity from pint package.

    Args:
        inp (object): Object to assert as length dimensional type

    """
    length_dimensionality = (1*unit('ft')).dimensionality
    if not (inp.dimensionality == length_dimensionality):
        raise TypeError("Input value is not correctly typed! Use length"
                        " dimensional unit.")


def check_US_length_units(ell):
    """Check if length unit type is an US unit

    Args:
        ell (length): Length unit of type Quantity from pint package.

    Returns:
        bool: True if US unit else False
    """
    return ell.units in _US_length_units


def check_area_dimensioned(inp):
    """Check that input is area type Quantity from pint package.

    Args:
        inp (object): Object to assert as area dimensional type

    """
    area_dimensionality = (1*unit('ft^2')).dimensionality
    if not (inp.dimensionality == area_dimensionality):
        raise TypeError("Input value is not correctly typed! Use area"
                        " dimensional unit.")


def name_of_var(var):
    """Find name of variables using local items.

    Args:
        var (unit): Dimensioned variable

    Returns:
        str: User-coded name of variable
    """
    try:
        local_vars = currentframe().f_back.f_back.f_locals.items()
        match = [name for name, val in local_vars if val is var]
    except AttributeError:
        local_vars = currentframe().f_back.f_locals.items()
        match = [name for name, val in local_vars if val is var]
    name = match[0] if match else "unknown"
    return name


def printv(var, to=None, name="", prec=".5g", fmt="~P", *args, **kwargs):
    """Print name and value of a Pint dimensional variable.
    For example,

        distance = 99.9 * unit('m')
        printv(distance)
        # prints "distance = 99.9 m"

    *Note*: As of Python 3.8, simply use the f-string syntax, e.g.
        x=7
        print(f"{x=}")

    Args:
        var (unit): Variable to be printed
        to (str): Convert to another unit
        name (str): Overwrite variable name
        prec (str): Precision formatter
        fmt (str): Additional formatter for pretty print or lack thereof
        *args: Additional arguments
        **kwargs: Additional keyword arguments

    """
    if isinstance(var, unit.Quantity):
        output = (var.to(to) if to is not None else var.to_base_units())
    else:
        output = var
    name = name if name else name_of_var(var)
    print(f"{name} = {output:{prec}{fmt}}", *args, **kwargs)


def to_base_units_wrapper(func):
    """Function decorator to convert output variable units to base units.

    Args:
        func (callable): Function to wrap

    Returns:
        callable: Called, wrapped function
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        output = func(*args, **kwargs)
        output.ito_base_units()
        return output
    return wrapper
