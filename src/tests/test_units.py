#!/usr/bin/env python
"""
Test unit capabilities.

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

# flake8: noqa E203

import difflib
import pytest

# from importlib_resources import files  # for older versions of Python
from importlib.resources import files
from pathlib import Path

import pint

from flightcondition import Atmosphere, FlightCondition, unit

def diff(file1, file2):
    """Find differences in two files. """
    with open(file1, 'r') as file:
        lines1 = file.readlines()
    with open(file2, 'r') as file:
        lines2 = file.readlines()

    diff_lines = []
    for line in difflib.unified_diff(lines1, lines2, lineterm=''):
        diff_lines.append(line)

    return diff_lines

def test_custom_units():
    """Test that custom units database is only different from default pint
    database by one line. """
    pint_dir = Path(pint.__file__).parent
    pint_constants_en = pint_dir / "constants_en.txt"
    pint_default_en = pint_dir / "default_en.txt"
    fc_constants_en = files("flightcondition").joinpath("data/constants_en.txt")
    fc_units_en = files("flightcondition").joinpath("data/fc_units_en.txt")

    constants_diff_lines = diff(pint_constants_en, fc_constants_en)
    assert constants_diff_lines == []

    truth_diff_lines = [
        '--- ', '+++ ', '@@ -867,6 +867,6 @@', ' @end\n', ' \n',
        (' @system US using USCSLiquidVolume, USCSDryVolume, USCSVolumeOther, '
        'USCSLengthInternational, USCSLengthSurvey, AvoirdupoisUS\n'),
        '-    yard\n', '+    foot\n', '     pound\n', ' @end\n']
    units_diff_lines = diff(pint_default_en, fc_units_en)
    assert units_diff_lines == truth_diff_lines

def test_feet_units():
    """Test that feet are properly set for custom 'US' units. """

    fc = FlightCondition(h=10*unit('kft'), M=0.55, L=2*unit('m'))
    assert str(fc.h.units) == "ft"
    assert str(fc.TAS.units) == "ft/s"

    fc = FlightCondition(h=10*unit('km'), M=0.55, L=2*unit('ft'))
    assert str(fc.h.units) == "m"

    fc.units = 'US'
    assert str(fc.h.units) == "ft"
    assert str(fc.TAS.units) == "ft/s"
