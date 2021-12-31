*******************
Aerospace Utilities
*******************

About
=====

This module contains commonly-used aerospace utilities for problem solving.

    * **Flight Condition**: input altitude to compute common flight condition
      data.  Easily swap between true airspeed, calibrated airspeed, and
      equivalent airspeed.  Includes atmospheric data.
    * **Atmosphere**: input altitude to compute 1993 International Standard
      Atmosphere data.  Additional derived quantities also included.


Author
======

Matthew C. Jones <matt.c.jones.aoe@gmail.com>

Installation
============

Install Commands
----------------

Install using the *pip* package-management system.  The easiest method is to
open the terminal and run:

.. code-block:: bash

    pip install aeroutils

Alternatively, manually download the `source code
<https://github.com/MattCJones/aeroutils>`, unpack, and run:

.. code-block:: bash

    pip install <path/to/aeroutils>

Dependencies
------------

* `numpy <https://numpy.org>`_: widely-used package for scientific computing.
* `pint <https://pint.readthedocs.io>`_: package for dealing with units.

Usage
=====
In a Python script or an ipython notebook, import all utilities with:

.. code-block:: python

    from aeroutils import *

Flight Condition Package
------------------------

The :code:`flightcondition` package can be used to compute and interact with
common flight condition data.

Outputs include:

    * Mach number :code:`mach`
    * True airspeed :code:`TAS`
    * Calibrated airspeed :code:`CAS`
    * Equivalent airspeed :code:`EAS`
    * Dynamic pressure :code:`q_inf`
    * Reynolds number :code:`reynolds_number(ell)`
    * Reynolds number per unit length
      :code:`reynolds_number_by_unit_length(length_unit)`
    * Atmosphere data :code:`atm` (see :code:`atmosphere` below) 

Usage:

.. code-block:: python

    from aeroutils.flightcondition import FlightCondition
    from aeroutils.units import *

    # Compute flight conditions for a scalar or array of altitudes
    altitudes = [0, 10e3, 33.5e3] * unit('ft')
    fc = FlightCondition(altitudes, EAS=300*unit('knots'))
    print(f"\nAirspeed in multiple formats: {fc}")
    print(f"Even more data: {fc.tostring()}")
    print(f"Access atmospheric data (see Atmosphere class): {fc.atm}")

    # Or view fc formats individually:")
    print(f"\nThe Mach number is {fc.mach:.5g}")
    print(f"The true fc is {fc.TAS:.5g}")
    print(f"The calibrated fc is {fc.CAS:.5g}")
    print(f"The equivalent fc is {fc.EAS:.5g}")

    # Define flight condition with Mach number, TAS, CAS, or EAS:")
    fc = FlightCondition(altitudes, mach=0.4535*dimless)
    fc = FlightCondition(altitudes, TAS=300*unit('knots'))
    fc = FlightCondition(altitudes, CAS=300*unit('knots'))
    fc = FlightCondition(altitudes, EAS=300*unit('knots'))

    # Compute flight condition data based on input length scale
    ell = 5 * unit('ft')
    print(f"\nThe Reynolds number is {fc.reynolds_number(ell):.5g}")
    print(f"The Reynolds number per unit length is "
        f"{fc.reynolds_number_by_unit_length('in'):.5g}")

    # Use unit functionality to convert dimensions as desired:")
    print(f"\nThe dynamic pressure is {fc.q_inf.to('psi'):.5g}")

Atmosphere Package
------------------------

The code:`atmosphere` package can be used to compute and interact with common
standard atmosphere data and derived quantities.

Outputs include:

    * Pressure :code:`p`
    * Temperature :code:`T`
    * Density :code:`rho`
    * Sound speed :code:`a`
    * Dynamic viscosity :code:`mu`
    * Kinematic viscosity :code:`nu`
    * Thermal conductivity :code:`k`
    * Layer name :code:`layer.name`
    * Geometric altitude :code:`h`
    * Geopotential altitude :code:`H`
    * Acceleration due to gravity :code:`g`
    * Mean free path :code:`mean_free_path`

Usage:

.. code-block:: python

    from aeroutils.atmosphere import Atmosphere
    from aeroutils.units import *

    # Compute atmospheric data for a scalar or array of altitudes
    h = [0.0, 12.0, 33.5] * unit('km')
    atm = Atmosphere(h)
    print(f"Abbreviated output: {atm}")
    print(f"Extended output in Imperial units: "
        f"{atm.tostring(short_repr=False, imperial_units=False)}")

    # Access individual properties and convert to desired units: "
    print(f"\np={atm.p}\nT={atm.T.to('degR')}\nrho={atm.rho.to('kg/m^3')}")

    # Compute properties such as thermal conductivity, mean free path and
    # many more!
    print(f"\nthermal conductivity k={atm.k}"
        f"\nmean free path = {atm.mean_free_path} and many more!")

License
=======

aeroutils is licensed under the MIT LICENSE. See the LICENSE document.

Disclaimer
=======
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
