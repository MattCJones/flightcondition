# Flight Condition

## About

Airspeed conversions (true/calibrated/equivalent/Mach), atmospheric data, and
more with built-in unit checking.  Specific sub-modules include:

* `flightcondition`: input altitude to compute common flight condition data.
  Easily swap between true airspeed, calibrated airspeed, equivalent airspeed,
  and Mach number.  Includes atmospheric data.
* `atmosphere`: input altitude to compute
  [1993 International Standard Atmosphere](
  https://en.wikipedia.org/wiki/International_Standard_Atmosphere)
  data.  Many relevant, derived quantities are included. The upper limit is 80
  kilometers.
* `units`: built-in unit-checking and conversion using the
  [pint](https://pint.readthedocs.io) package.

![Flight Condition Demo](
https://github.com/MattCJones/videos/blob/main/flightcondition/flightcondition_demo.gif
)

## Author

Matthew C. Jones [matt.c.jones.aoe@gmail.com](matt.c.jones.aoe@gmail.com)

## Installation

### Install Commands

Install using the `pip` package-management system. The
easiest method is to open the terminal and run:

```python
pip install flightcondition
```

Alternatively, manually download the
[source code](https://github.com/MattCJones/flightcondition), unpack, and run:

```python
pip install <path/to/flightcondition>
```

### Dependencies

-   [numpy](https://numpy.org): package for scientific computing.
-   [pint](https://pint.readthedocs.io): package for dealing with units.

## Usage

Import all utilities with,

```python
from flightcondition import *
```

or more explicitly as shown in the following examples.

### Flight Condition

The `FlightCondition` class is used to compute and interact with common flight
condition data. Inputs include altitude, velocity in some format, and an
optional length scale.

#### Input Arguments
Input arguments include:

1.  **Altitude**: `h` (aliases include: `alt`, `altitude`)
2.  **Velocity** (pick one):
    -   *True airspeed*: `TAS` (aliases include: `tas`, `true_airspeed`,
        `U_inf`, `V_inf`)
    -   *Calibrated airspeed*: `CAS` (aliases include: `cas`,
        `calibrated_airspeed`)
    -   *Equivalent airspeed*: `EAS` (aliases include: `eas`,
        `equivalent_airspeed`)
    -   *Mach number*: `M` (aliases include: `mach`, `Mach`, `M_inf`,
        `mach_number`)
3.  **Length-scale** (optional): `L` (aliases include: `ell`, `bylen`,
    `length`, `length_scale`, `l`)

Input quantities must be dimensionalized - see the usage below.  Alternatively
use `KTAS`, `KCAS`, or `KEAS` for convenience.  For example, `KCAS=233`
is equivalent to `CAS=233*unit('knots')`.

#### Output Quantities

The following tables list the quantities and the variables used to access them.
Quantities may be accessed by either (a) their abbreviated variable, e.g.
`.TAS`, or (b) by their full names, e.g.  `byname.true_airspeed`. They
may also be accessed through their particular sub-category: `byalt`, `byvel`,
or `bylen`, e.g.  `.byvel.TAS` or `.byvel.byname.true_airspeed`.

| Altitude Quantity     | Variable | Full Name (via `.byname.`) |
|-----------------------|----------|----------------------------|
| Geometric altitude    | `h`      | `geometric_altitude`       |
| Geopotential altitude | `H`      | `geopotential_altitude`    |
| Pressure              | `p`      | `pressure`                 |
| Temperature           | `T`      | `temperature`              |
| Density               | `rho`    | `density`                  |
| Sound speed           | `a`      | `sounds_speed`             |
| Dynamic viscosity     | `mu`     | `dynamic_viscosity`        |
| Kinematic viscosity   | `nu`     | `kinematic_viscosity`      |
| Thermal conductivity  | `k`      | `thermal_conductivity`     |
| Gravity               | `g`      | `gravity`                  |
| Mean free path        | `MFP`    | `mean_free_path`           |

| Velocity Quantity                | Name      | Full Name (via `.byname.`)       |
|----------------------------------|-----------|----------------------------------|
| True airspeed (TAS)              | `TAS`     | `true_airspeed`                  |
| Calibrated airspeed (CAS)        | `CAS`     | `calibrated_airspeed`            |
| Equivalent airspeed (EAS)        | `EAS`     | `equivalent_airspeed`            |
| Mach number                      | `M`       | `mach_number`                    |
| Mach angle                       | `mu_M`    | `mach_angle`                     |
| Dynamic pressure                 | `q_inf`   | `dynamic_pressure`               |
| Impact pressure                  | `q_c`     | `impact_pressure`                |
| Stagnation pressure              | `q_0`     | `stagnation_pressure`            |
| Stagnation temperature           | `T_0`     | `stagnation_temperature`         |
| Recovery temperature (laminar)   | `Tr_lamr` | `recovery_temperature_laminar`   |
| Recovery temperature (turbulent) | `Tr_turb` | `recovery_temperature_turbulent` |
| Reynolds number per unit length  | `Re_by_L` | `recovery_temperature_turbulent` |

| Length-Scale Quantity                            | Name        | Full Name (via `.byname.`)       |
|--------------------------------------------------|-------------|----------------------------------|
| Length scale                                     | `L`         | `length_scale`                   |
| Reynolds number                                  | `Re`        | `reynolds_number`                |
| Boundary layer thickness (laminar)               | `h_BL_lamr` | `boundary_thickness_laminar`     |
| Boundary layer thickness (turbulent)             | `h_BL_turb` | `boundary_thickness_turbulent`   |
| Flat plate skin friction coefficient (laminar)   | `Cf_lamr`   | `friction_coefficient_laminar`   |
| Flat plate skin friction coefficient (turbulent) | `Cf_turb`   | `friction_coefficient_turbulent` |

#### Example Usage

```python
from flightcondition import FlightCondition, unit

# Compute flight condition at 3 km, Mach 0.5
fc = FlightCondition(h=3*unit('km'), M=0.5)

# Uncomment to print summary of flight condition quantities:
#print(f"{fc}")

# Uncomment to print abbreviated output in US units:
#print(f"\n{fc.tostring(full_output=False, units="US")}")

# Convert true, calibrated, equivalent airspeeds
KTAS = fc.TAS.to('knots')
KCAS = fc.CAS.to('knots')
KEAS = fc.EAS.to('knots')
print(f"Flying at {KTAS.magnitude:.4g} KTAS,"
      f" which is {KCAS.magnitude:.4g} KCAS,"
      f" or {KEAS.magnitude:.4g} KEAS")
# >>> Flying at 319.4 KTAS, which is 277.7 KCAS, or 275.1 KEAS

# Access atmospheric data (see Atmosphere class for more)
h, p, T, rho, nu, a = fc.h, fc.p, fc.T, fc.rho, fc.nu, fc.a
print(f"The ambient temperature at {h.to('km'):.4g} is {T:.4g}")
# >>> The ambient temperature at 3 km is 268.7 K

# Change airspeed to 300 KEAS and altitude to 12 kft
fc.EAS = 300 * unit('knots')
fc.h = 12 * unit('kft')
#print(f"{fc}")  # uncomment to print output

# Recompute for a range of altitudes at 275.14 knots-equivalent
# airspeed with a characteristic length scale of 10 meters
fc = FlightCondition(h=[0, 9.8425, 20]*unit('kft'),
                    EAS=275.14*unit('kt'),
                    L=10*unit('m'))

# Compute additional derived quantities - explore the class for more!
print(f"\nThe dynamic pressure in psi is {fc.q_inf.to('psi'):.3g}")
# >>> The dynamic pressure in psi is [1.78 1.78 1.78] psi
print(f"The Reynolds number is {fc.Re:.3g}")
# >>> The Reynolds number is [9.69e+07 8.82e+07 7.95e+07]
h_yplus100 = fc.wall_distance_from_yplus(100)
print(f"The wall distance where y+=100 is {h_yplus100.to('in'):.3g}")
# >>> The wall distance where y+=100 is [0.0126 0.0138 0.0153] in

# Alternatively access quantities by their full name
print(fc.TAS == fc.byname.true_airspeed)
# >>> [ True  True  True]

# Or by their sub-categories: `byalt`, `byvel`, or `bylen`
print(fc.byvel.TAS == fc.byvel.byname.true_airspeed)
# >>> [ True  True  True]
```

### Atmosphere

The `Atmosphere` class can be used to compute and interact with common standard
atmosphere data and derived quantities.  See the list of output quantities in
the `FlightCondition` documentation above.
See also `layer` for layer properties such as `layer.name` for the layer name.

Note that all `Atmosphere` properties can be accessed through the
`FlightCondition` class, however, this class stands on its own if the
additional velocity and length-scale quantities are not desired.

#### Input Arguments
The input argument is geometric altitude `h`.  Aliases include `alt` and
`altitude`.  Note that these geometric altitude must be input as dimensional
length quantities - see the usage below.  Alternatively input
un-dimensionalized numbers using `h_kft` or `h_km` for kilofeet and kilometers
respectively.

#### Example Usage

```python
from flightcondition import Atmosphere, unit

# Compute atmospheric data for a scalar or array of altitudes
h = [0.0, 44.2, 81.0] * unit('km')
atm = Atmosphere(h)

# Uncomment to print all atmospheric quantities:
#print(f"\n{atm}")

# Uncomment to print while specifying abbreviated output in US units:
#print(f"\n{atm.tostring(full_output=False, units="US")}")

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
```

### Units

Conveniently input, output, and convert units using
[pint](https://pint.readthedocs.io) units.

```python
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
```

Note that [pint](https://pint.readthedocs.io) does not support conflicting unit
registries so avoid interactions between `flightcondition.unit` and a separate
`pint.UnitRegistry`.

### Command Line Interface

A command line interface (CLI) is included for convenience but with limited
functionality. Run `flightcondition -h` for help.

An example call is given for the flight condition of 233
knots-equivalent-airspeed at 23 kilofeet with a length scale of 4 feet and
abbreviated output:

```bash
flightcondition --h 23 kft --EAS 233 knots --L 4 ft --no-full-output
```

```bash
===========================================================
   Flight Condition (units=US, full_output=False)
===========================================================
------------------  Altitude Quantities  ------------------
geometric_altitude  h       = 23 kft
pressure            p       = 857.25 lbf/ft²
temperature         T       = 436.74 °R
density             rho     = 1.1435×10⁻³ slug/ft³
sound_speed         a       = 1024.5 ft/s
kinematic_viscosity nu      = 2.8509×10⁻⁴ ft²/s
------------------  Velocity Quantities  ------------------
mach_number         M       = 0.55344
true_airspeed       TAS     = 335.93 kt
calibrated_airspeed CAS     = 238.14 kt
equivalent_airspeed EAS     = 233 kt
reynolds_per_length Re_by_L = 1.6573×10⁵ 1/in
------------------   Length Quantities   ------------------
length_scale        L       = 4 ft
reynolds_number     Re      = 7.9551×10⁶
```

Alternatively use the `--KEAS 233` syntactic sugar to omit the `knots` unit.
See also `--KTAS` and `--KCAS`.

## Assumptions

-   Atmospheric quantities follow the
    [1993 International Standard Atmosphere](
    https://en.wikipedia.org/wiki/International_Standard_Atmosphere)
    model.
-   Velocity computations include varying degrees of the following assumptions.
    Note that several assumptions break down for hypersonic flow.
    -   Continuum flow (mean free path is much smaller than the characteristic
        length scale)
    -   Ideal gas
    -   Thermally perfect gas
    -   Calorically perfect gas
    -   Adiabatic
    -   Reversible (`CAS`, `q_c`, `p_0`)

## Unit Testing
`flightcondition` is maintained using unit tests to maintain functionality and
accuracy.  Even so, see the Disclaimer below.

## License

`flightcondition` is licensed under the MIT LICENSE. See the
[LICENSE](https://github.com/MattCJones/flightcondition/blob/main/LICENSE)
document.

## Disclaimer

The software is provided "as is", without warranty of any kind, express or
implied, including but not limited to the warranties of merchantability,
fitness for a particular purpose and noninfringement. In no event shall the
authors or copyright holders be liable for any claim, damages or other
liability, whether in an action of contract, tort or otherwise, arising from,
out of or in connection with the software or the use or other dealings in the
software.
