****************
Flight Condition
****************

About
=====

Airspeed conversions (true/calibrated/equivalent/Mach), atmospheric data, and
more with built-in unit checking.  Specific sub-modules include:

* :code:`flightcondition` : input altitude to compute common flight condition
  data.  Easily swap between true airspeed, calibrated airspeed,
  equivalent airspeed, and Mach number.  Includes atmospheric data.
* :code:`atmosphere` : input altitude to compute `1993 International Standard
  Atmosphere
  <https://en.wikipedia.org/wiki/International_Standard_Atmosphere>`_ data.
  Many relevant, derived quantities are included. The upper limit is 80
  kilometers.
* :code:`units` : built-in unit-checking and conversion using `pint
  <https://pint.readthedocs.io>`_ package.

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
Import all utilities with,

.. code-block:: python

    from flightcondition import *

or more explicitly as shown in the following examples.


Flight Condition
----------------

The :code:`FlightCondition` class is used to compute and interact with
common flight condition data.  Inputs include altitude, airspeed in some
format, and an optional length scale.

**Input** arguments include:

#. :code:`h` *geometric altitude* - aliases are :code:`alt`, :code:`altitude`

#. Velocity (pick one):

   * :code:`TAS` *true airspeed* - aliases are :code:`tas`,
     :code:`true_airspeed`, :code:`U_inf`, :code:`V_inf`
   * :code:`CAS` *calibrated airspeed* - aliases are
     :code:`cas`, :code:`calibrated_airspeed`
   * :code:`EAS` *equivalent airspeed* - aliases are :code:`eas`,
     :code:`equivalent_airspeed`
   * :code:`M` *mach number* - aliases are :code:`mach`, :code:`Mach`,
     :code:`M_inf`, :code:`mach_number`

#. :code:`L` *length scale* (optional) - aliases are :code:`ell`, :code:`len`,
   :code:`length`, :code:`length_scale`, :code:`l`

See also the :code:`KTAS`, :code:`KCAS`, and :code:`KEAS` syntactic sugar.  For
example, :code:`KCAS=233` is equivalent to :code:`CAS=233*unit('knots')`.

**Outputs** include:

#. :code:`atm` *atmospheric* quantities - see :code:`Atmosphere` class below.
#. :code:`vel` *airspeed* quantities:

   * True airspeed :code:`TAS`
   * Calibrated airspeed :code:`CAS`
   * Equivalent airspeed :code:`EAS`
   * Mach number :code:`M`
   * Dynamic pressure :code:`q_inf`
   * Impact pressure :code:`q_c`
   * Stagnation pressure :code:`p_0`
   * Stagnation temperature :code:`T_0`
   * Reynolds number per unit length :code:`Re_by_L`

#. :code:`len` *length-scale* quantities:
   
   * Reynolds number :code:`Re`

Quantities may be accessed using the specific :code:`atm`, :code:`vel`, or
:code:`len` object.  Quantities may also be accessed using their full name with
the :code:`byname` object.  For example, Mach number can be accessed using
:code:`.vel.M` or by its full name using :code:`.byname.mach_number`

**Example usage**:

.. code-block:: python

    from flightcondition import FlightCondition, unit, dimless

    # Compute flight condition at 3 km, Mach 0.5
    fc = FlightCondition(3*unit('km'), M=0.5)

    # Uncomment to print summary of flight condition quantities:
    #print(f"{fc}")

    # Uncomment to print abbreviated output in US units:
    #print(f"\n{fc.tostring(full_output=False, unit_system="US")}")

    # Access true, calibrated, equivalent airspeeds
    KTAS = fc.vel.TAS.to('knots')
    KCAS = fc.vel.CAS.to('knots')
    KEAS = fc.vel.EAS.to('knots')
    print(f"Flying at {KTAS.magnitude:.4g} KTAS,"
        f" which is {KCAS.magnitude:.4g} KCAS,"
        f" or {KEAS.magnitude:.4g} KEAS")
    # >>> Flying at 319.4 KTAS, which is 277.7 KCAS, or 275.1 KEAS

    # Access atmospheric data (see Atmosphere class for more)
    atm = fc.atm  # access Atmosphere object
    h, p, T, rho, nu, a = atm.h, atm.p, atm.T, atm.rho, atm.nu, atm.a
    print(f"The ambient temperature at {h.to('km'):.4g} is {T:.4g}")
    # >>> The ambient temperature at 3 km is 268.7 K

    # Compute again instead using true airspeed and altitude in km
    fc = FlightCondition(3.048*unit('km'), TAS=401.7*unit('mph'))
    #print(f"{fc}")  # uncomment to print output

    # Compute for a range of altitudes at 275.14 knots-equivalent
    # airspeed with a characteristic length scale of 10 meters
    fc = FlightCondition([0, 9.8425, 20]*unit('kft'),
                        EAS=275.14*unit('kt'),
                        L=10*unit('m'))

    # Compute additional derived quantities
    # Explore the class data structure for all options
    print(f"\nThe dynamic pressure in psi is "
        f"{fc.vel.q_inf.to('psi'):.3g}")
    # >>> The dynamic pressure in psi is [1.78 1.78 1.78] psi
    print(f"The Reynolds number is {fc.len.Re:.3g}")
    # >>> The Reynolds number is [9.69e+07 8.82e+07 7.95e+07]

    # Alternatively access quantities by their full name
    print(fc.vel.TAS == fc.byname.true_airspeed)
    # >>> [ True  True  True]


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
* Mean free path :code:`MFP`

Usage:

.. code-block:: python

    from flightcondition import Atmosphere, unit

    # Compute atmospheric data for a scalar or array of altitudes
    h = [0.0, 44.2, 81.0] * unit('km')
    atm = Atmosphere(h)

    # Uncomment to print all atmospheric quantities:
    #print(f"\n{atm}")

    # Uncomment to print while specifying abbreviated output in US units:
    #print(f"\n{atm.tostring(full_output=False, unit_system="US")}")

    # See also the linspace() function from numpy, e.g.
    # h = linspace(0, 81.0, 82) * unit('km')

    # Access individual properties and convert to desired units: "
    p, T, rho, nu, a, k = atm.p, atm.T, atm.rho, atm.nu, atm.a, atm.k
    print(f"\nThe pressure in psi is {p.to('psi'):.3g}")
    # >>> The pressure in psi is [14.7 0.024 0.000129] psi

    # Compute additional properties such as mean free path
    # Explore the class data structure for all options
    print( f"\nThe mean free path = {atm.MFP:.3g}")
    # >>> The mean free path = [7.25e-08 4.04e-05 0.00564] yd

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

Command Line Interface
----------------------
A command line interface (CLI) is included for convenience but with limited
functionality.  Run :code:`flightcondition -h` for help.

An example call is given for the flight condition of 233
knots-equivalent-airspeed at 23 kilofeet with a length scale of 4 feet and
abbreviated output:

.. code-block:: bash

    flightcondition --alt 23 kft --EAS 233 knots --len 4 ft --abbreviated

.. code-block:: bash

    =========================================================
        Flight Condition (US units, abbreviated output)
    =========================================================
    -----------------  Altitude Quantities  -----------------
    geometric_altitude    h       = 23 kft
    pressure              p       = 857.25 lbf/ft²
    temperature           T       = 436.74 °R
    density               rho     = 1.1435×10⁻³ slug/ft³
    sound_speed           a       = 1024.5 ft/s
    kinematic_viscosity   nu      = 2.8509×10⁻⁴ ft²/s
    -----------------  Velocity Quantities  -----------------
    true_airspeed         TAS     = 566.99 ft/s
    calibrated_airspeed   CAS     = 401.94 ft/s
    equivalent_airspeed   EAS     = 393.26 ft/s
    mach_number           M       = 0.55344
    reynolds_per_length   Re_by_L = 1.6573×10⁵ 1/in
    -----------------   Length Quantities   -----------------
    length_scale          L       = 4 ft
    reynolds_number       Re      = 7.9551×10⁶

Alternatively use the :code:`--KEAS 233` syntactic sugar to omit the
:code:`knots` unit.  See also :code:`--KTAS` and  :code:`--KCAS`.

Assumptions
===========

* Atmospheric quantities follow the `1993 International Standard Atmosphere
  <https://en.wikipedia.org/wiki/International_Standard_Atmosphere>`_ model.
* Velocity computations include varying degrees of the following assumptions.
  If assumptions are broken for a particular quantity, that quantity returns
  :code:`nan`.

  - Continuum flow (mean free path is much smaller than the characteristic
    length scale)
  - Ideal gas
  - Thermally perfect gas
  - Calorically perfect gas
  - Adiabatic
  - Reversible (:code:`CAS`, :code:`q_c`, :code:`p_0`)

License
=======

:code:`flightcondition` is licensed under the MIT LICENSE. See the `LICENSE
<https://github.com/MattCJones/flightcondition/blob/main/LICENSE>`_ document.

Disclaimer
==========
The software is provided "as is", without warranty of any kind, express or
implied, including but not limited to the warranties of merchantability,
fitness for a particular purpose and noninfringement. In no event shall the
authors or copyright holders be liable for any claim, damages or other
liability, whether in an action of contract, tort or otherwise, arising from,
out of or in connection with the software or the use or other dealings in the
software.
