#!/usr/bin/env python
"""Provide dimensional units for mathematical operations using the pint
package.

Dependencies: pint

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

from functools import wraps
# from importlib.resources import files  # Python 3.9+
from pathlib import Path
from pkg_resources import resource_filename  # Python 3.8+
from inspect import currentframe
from warnings import warn

from pint import UnitRegistry, DefinitionSyntaxError

__all__ = ['unit', 'dimless', 'printv']

# For Python 3.9+:
# fc_units_file = files("flightcondition").joinpath("data/fc_units_en.txt")
# For Python 3.8+:
fc_units_file = Path(
    resource_filename("flightcondition", "data/fc_units_en.txt"))  # 3.8+
try:
    unit = UnitRegistry(str(fc_units_file))
except DefinitionSyntaxError:
    warn(f"Failed to load custom unit system:\n{fc_units_file}\n"
         "Loading default 'pint' unit system instead")
    unit = UnitRegistry(system='mks')

unit.default_format = '~P'
dimless = unit('dimensionless')

_US_length_units = ('ft', 'feet', 'foot', 'kft', 'kilofoot', 'kilofeet'
                    'fts', 'feets', 'foots', 'kfts', 'kilofoots', 'kilofeets',
                    'mi', 'mile', 'miles')


def check_same_registry(inp):
    """Check if input units are using this package's unit registry.

    Args:
        inp (object): Object to assert as dimensional type

    """
    if isinstance(inp, unit.Quantity):  # method to check if same registry
        dummy = 1 * unit('m')
        # Registry equality check
        return (inp._REGISTRY is dummy._REGISTRY)
    else:
        return False


def check_dimensioned(inp):
    """Check that input is dimensional (type 'Quantity' from pint package).

    Args:
        inp (object): Object to assert as dimensional type

    """
    try:
        inp.units.compatible_units()  # check that uses pint methods
    except AttributeError:
        raise TypeError("Input value is not correctly dimensioned! Use"
                        " dimensional type 'Quantity' from pint package.")


def check_dimensionless(inp):
    """Check that input is dimensionless type Quantity from pint package.

    Args:
        inp (object): Object to assert as length dimensional type

    """
    if not inp.check('[]'):
        raise TypeError("Input value is not correctly typed! Use dimensionless"
                        " dimensional unit.")


def check_length_dimensioned(inp):
    """Check that input is length dimension type Quantity from pint package.

    Args:
        inp (object): Object to assert as length dimensional type

    """
    if not inp.check('[length]'):
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
    if not inp.check('[length]^2'):
        raise TypeError("Input value is not correctly typed! Use area"
                        " dimensional unit.")


def check_pressure_dimensioned(inp):
    """Check that input is pressure type Quantity from pint package.

    Args:
        inp (object): Object to assert as pressure dimensional type

    """
    # p = F/A = ma/A = M L /T2 /L2 = M /L /T2
    if not inp.check('[mass] [length]^-1 [time]^-2'):
        raise TypeError("Input value is not correctly typed! Use pressure"
                        " dimensional unit.")


def check_temperature_dimensioned(inp):
    """Check that input is temperature type Quantity from pint package.

    Args:
        inp (object): Object to assert as temperature dimensional type

    """
    if not inp.check('[temperature]'):
        raise TypeError("Input value is not correctly typed! Use temperature"
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
        fmt = ""
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
        if output is not None:
            output.ito_base_units()
        return output
    return wrapper
