****************
Flight Condition
****************

About
=====

Airspeed conversions (true/calibrated/equivalent/Mach), atmospheric data, and
more with built-in unit checking.  Specific sub-modules include:

* :code:`flightcondition` : input altitude to compute common flight condition data.  Easily swap between Mach number, true airspeed, calibrated airspeed, and equivalent airspeed.  Includes atmospheric data.
* :code:`atmosphere` : input altitude to compute 1993 International Standard Atmosphere data.  Many relevant, derived quantities are included.
* :code:`units` : built-in unit-checking and conversion using `pint <https://pint.readthedocs.io>`_ package.

See examples below for usage!


Author
======

Matthew C. Jones <matt.c.jones.aoe@gmail.com>

Installation
============

Install Commands
----------------

Install using the :code:`pip` package-management system.  The easiest method is
to open the terminal and run:

.. code-block:: bash

    pip install flightcondition

Alternatively, manually download the `source code
<https://github.com/MattCJones/flightcondition>`_, unpack, and run:

.. code-block:: bash

    pip install <path/to/flightcondition>

Dependencies
------------

* `numpy <https://numpy.org>`_: package for scientific computing.

* `pint <https://pint.readthedocs.io>`_: package for dealing with units.

Usage
=====
In a Python script or an ipython notebook, import all utilities with,

.. code-block:: python

    from flightcondition import *

or more explicitly as shown in the following examples.


Flight Condition
----------------

The :code:`Flightcondition` class can be used to compute and interact with
common flight condition data.

Outputs include:

* Mach number :code:`mach`
* True airspeed :code:`TAS`
* Calibrated airspeed :code:`CAS`
* Equivalent airspeed :code:`EAS`
* Dynamic pressure :code:`q_inf`
* Reynolds number :code:`reynolds_number(ell)`
* Reynolds number per-unit-length :code:`reynolds_number_per_unit_length(length_unit)`
* Atmosphere data :code:`atm` (see :code:`atmosphere` below) 

Usage:

.. code-block:: python

    from flightcondition import FlightCondition, unit, dimless

    # Compute flight conditions for a scalar or array of altitudes
    altitudes = [0, 10, 33.5] * unit('kft')
    fc = FlightCondition(altitudes, EAS=300*unit('knots'))

    # Print flight condition data:
    print(f"{fc}")

    # Print extended output:
    print(f"\n{fc.tostring(full_output=True)}")

    # Access flight speed formats individually
    M_inf, U_inf, U_CAS, U_EAS = fc.mach, fc.TAS, fc.CAS, fc.EAS

    # Access atmospheric data alone (see Atmosphere class for options)
    atm = fc.atm  # access Atmosphere object 'atm'
    p, T, rho, nu, a = atm.p, atm.T, atm.rho, atm.nu, atm.a

    # Input true/calibrated/equivalent airspeed or Mach number
    fc_TAS = FlightCondition(altitudes, TAS=300*unit('knots'))
    fc_CAS = FlightCondition(altitudes, CAS=300*unit('knots'))
    fc_EAS = FlightCondition(altitudes, EAS=300*unit('knots'))
    fc_mach = FlightCondition(altitudes, mach=0.4535*dimless)

    # Specify desired units on input and output
    altitudes_in_km = [0, 3.048, 10.2108] * unit('km')
    fc_other_units = FlightCondition(altitudes, EAS=154.33*unit('m/s'))
    U_TAS = fc_other_units.TAS
    print(f"\nThe true airspeed in m/s is {U_TAS.to('m/s'):.5g}")
    print(f"The true airspeed in km/s is {U_TAS.to('km/s'):.5g}")

    # Compute additional derived quantities (see class for all options)
    print(f"\nThe dynamic pressure in psi is {fc.q_inf.to('psi'):.5g}")
    ell = 60 * unit('in')  # arbitrary length scale of interest
    print(f"The Reynolds number is {fc.reynolds_number(ell):.5g}")
    print(f"The Reynolds number per-unit-length [1/in] is "
        f"{fc.reynolds_number_per_unit_length('in'):.5g}")

Atmosphere
----------

The :code:`Atmosphere` class can be used to compute and interact with common
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

    from flightcondition import Atmosphere, unit

    # Compute atmospheric data for a scalar or array of altitudes
    h = [0.0, 12.7, 44.2, 81.0] * unit('km')
    atm = Atmosphere(h)

    # Print abbreviated output:
    print(f"\n{atm}")

    # Print extended output in US units:
    print(f"\n{atm.tostring(full_output=True, US_units=True)}")

    # See also the linspace() function from numpy, e.g.
    # h = linspace(0, 81.0, 82) * unit('km')

    # Access individual properties and convert to desired units: "
    p, T, rho, nu, a = atm.p, atm.T, atm.rho, atm.nu, atm.a
    print(f"\nThe pressure in psi is {p.to('psi'):.5g}")

    # Compute additional properties such as thermal conductivity,
    # mean free path, and more (see class for all options)
    print(f"\nThe thermal conductivity is {atm.k:.5g}"
        f"\nThe mean free path = {atm.mean_free_path:.5g}")

Units
-----

Conveniently input, output, and convert units using `pint
<https://pint.readthedocs.io>`_ units.

.. code-block:: python

    from flightcondition import unit, printv

    h = 33 * unit('km')
    print(h.to('kft'))
    # >>> 108.26771653543307 kft
    printv(h, to='kft')
    # >>> h = 108.27 kft

    U_inf = 20 * unit('knots')
    rho_inf = 1.225 * unit('kg/m^3')
    q_inf = 0.5*rho_inf*U_inf**2
    printv(q_inf, to='psi')
    # >>> q_inf = 0.0094042 psi

Note that `pint <https://pint.readthedocs.io>`_ does not support conflicting
unit registries so avoid interactions between :code:`flightcondition.unit` and
a separate :code:`pint.UnitRegistry`.


License
=======

:code:`flightcondition` is licensed under the MIT LICENSE. See the `LICENSE <https://github.com/MattCJones/flightcondition/blob/main/LICENSE>`_ document.

Disclaimer
==========
The software is provided "as is", without warranty of any kind, express or
implied, including but not limited to the warranties of merchantability,
fitness for a particular purpose and noninfringement. In no event shall the
authors or copyright holders be liable for any claim, damages or other
liability, whether in an action of contract, tort or otherwise, arising from,
out of or in connection with the software or the use or other dealings in the
software.
