#!/usr/bin/env python
"""Common classes and functions.

Dependencies: numpy, pint

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

from functools import wraps

import numpy as np

from flightcondition.units import to_base_units_wrapper


class AliasAttributes():
    """Nested alias class to reference quantities by prescribed aliases. """

    def __init__(self, varsobj_arr, varnames_dict_arr):
        """Populate full names and link to variable
        Args:
            varsobj (list): List of objects that holds all of the variables
            varnames_dict (list): List of dictionaries that maps variables to
                their alias names
        """
        # Avoid infinite loop when setting properties
        super().__setattr__("varsobj_arr", varsobj_arr)
        super().__setattr__("varnames_dict_arr", varnames_dict_arr)

    def __dir__(self):
        """Add tab completion for alias names. """
        # Avoid infinite loop when loading varsobj and varnames_dict
        varnames_dict_arr = super().__getattribute__("varnames_dict_arr")
        varnames_arr = []
        for varnames_dict in varnames_dict_arr:
            for varname in varnames_dict.values():
                varnames_arr.append(varname)
        return varnames_arr

    def __getattribute__(self, attr):
        """Get referenced quantity

        Args:
            attr (str): attribute name
        """
        # Avoid infinite loop when loading varsobj and varnames_dict
        varsobj_arr = super().__getattribute__("varsobj_arr")
        varnames_dict_arr = super().__getattribute__("varnames_dict_arr")
        for varsobj, varnames_dict in zip(varsobj_arr, varnames_dict_arr):
            for var, varname in varnames_dict.items():
                if attr == varname:
                    return getattr(varsobj, var)

    def __setattr__(self, attr, attrval):
        """Set referenced quantity

        Args:
            attr (str): attribute name
            attrval (quantity): attribute value
        """
        # Avoid infinite loop when loading varsobj and varnames_dict
        varsobj_arr = super().__getattribute__("varsobj_arr")
        varnames_dict_arr = super().__getattribute__("varnames_dict_arr")
        for varsobj, varnames_dict in zip(varsobj_arr, varnames_dict_arr):
            for var, varname in varnames_dict.items():
                if attr == varname:
                    return setattr(varsobj, var, attrval)


class DimensionalData:
    """Parent class to hold dimensional data"""

    varnames = {}

    def __eq__(self, otherobj):
        """Check equality.
        Returns:
            dict: Dictionary representation of object
        """
        if otherobj.__class__ is not self.__class__:
            return NotImplemented

        return self.asdict() == otherobj.asdict()

    def __str__(self):
        """Output string when object is printed.

        Returns:
            str: Full string output
        """
        return self.tostring(full_output=True)

    def __repr__(self):
        """Output string representation of class object.

        Returns:
            str: Full string output
        """
        # Determine full output flag
        if self.full_output is None:
            full_output = False
        else:
            full_output = self.full_output

        return self.tostring(full_output=full_output)

    def asdict(self):
        """Return class data as dictionary.
        Returns:
            dict: Class data
        """
        obj_dict = {}
        for var, varname in self.varnames.items():
            obj_dict[var] = getattr(self, var)
        return obj_dict

    def print(self, *args, **kwargs):
        """Print tostring() function to stdout. """
        print(self.tostring(*args, **kwargs))

    def tostring(self, full_output=True):
        """Override this function to output string representation of class
        object

        Args:
            full_output (bool): Set to True for full output

        Returns:
            str: String representation
        """
        return ""

    def _vartostr(self, var, var_str, to_units, max_var_chars=0,
                  fmt_val="10.5g", pretty_print=False):
        """Formatted variable string with variable, full name, and value.
        Shortens array output to fit better on screen.

        Args:
            var (Quantity): Variable
            var_str (str): Variable string
            to_units (str): Units to convert to
            fmt_val (str): String for formatting value
            pretty_print (bool): Pretty print format

        Returns:
            str: formatted string
        """
        pp_ = "~P" if pretty_print else ""
        var_name = self.varnames[var_str]
        var.ito(to_units)
        if np.size(var) > 6:
            mag = var.magnitude
            var_val_str = (f"[{mag[0]:{fmt_val}} {mag[1]:{fmt_val}} "
                           f"{mag[2]:{fmt_val}} ... {mag[-3]:{fmt_val}} "
                           f"{mag[-2]:{fmt_val}} {mag[-1]:{fmt_val}}] "
                           f"{var.units:{pp_}}")
        else:
            var_val_str = f"{var:{fmt_val}{pp_}}"
        var_str = f"{var_name:{max_var_chars}s} {var_str:10s} = {var_val_str}"
        return var_str

    @staticmethod
    def _arg_from_alias(alias_list, kwargs_dict):
        """Determine argument from hidden alias list

        Args:
            alias_list (dict): List of possible argument aliases

        Returns:
            unit.Quantity: argument
        """
        for alias in alias_list:
            if alias in kwargs_dict.keys():
                return kwargs_dict[alias]
        else:  # no break
            return None


def _len1array_to_scalar(func):
    """Decorator to output scalar if array is length 1.

    Args:
        func (object): Function to wrap

    Returns:
        object: Scalar or array
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        output = func(*args, **kwargs)
        if isinstance(output.magnitude, np.ndarray):
            if np.size(output) == 1:
                return output[0]
        return output
    return wrapper


def _property_decorators(func):
    """ Combine multiple decorators for properties of the class.

    Equivalent to:

        @property
        @to_base_units_wrapper
        @_len1array_to_scalar
    """
    return property(to_base_units_wrapper(_len1array_to_scalar(func)))
