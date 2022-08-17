"""
Aerospace utilities.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""
import importlib.metadata

from .atmosphere import Atmosphere
from .flightcondition import FlightCondition
from .boundarylayer import BoundaryLayer
from .isentropicflow import IsentropicFlow
from .normalshock import NormalShock
from .units import unit, dimless, printv

__version__ = importlib.metadata.version('flightcondition')
