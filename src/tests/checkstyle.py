#!/usr/bin/env python3
"""
Check that PEP8 format is followed

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

import os

from pathlib import Path
from shlex import split
from shutil import which
from subprocess import run


def check_format(py_file_path):
    """Check format of a Python file. """
    run_str = f"pycodestyle -v {py_file_path}"
    run(split(run_str))
    print("")


os.chdir(Path(__file__).parent)  # run from script directory

print("="*60)
print("Running pycodestyle")
print("="*60)
check_format("../aeroutils/atmosphere/atmosphere.py")
check_format("../aeroutils/flightcondition/flightcondition.py")
check_format("../aeroutils/units/units.py")
check_format("../aeroutils/constants/constants.py")
check_format("runtests.py")
check_format("checkusage.py")

if which('flake8'):
    print("="*60)
    print("Running flake8")
    print("="*60)
    os.chdir("../../src")
    run("flake8")
