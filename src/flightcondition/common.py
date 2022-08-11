#!/usr/bin/env python
"""Common classes and functions.

Dependencies:

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
        return self.tostring(full_output=False)

    @staticmethod
    def _arg_from_alias(alias_list, kwargs_dict):
        """Determine argument from hidden alias list

        Args:
            alias_list (dict): list of possible argument aliases

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
