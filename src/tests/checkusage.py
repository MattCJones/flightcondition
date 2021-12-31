#!/usr/bin/env python

import sys
# sys.path.append("../")
from pathlib import Path
package_root_directory = Path(__file__).resolve().parents[1]
sys.path.append(str(package_root_directory))


print("="*60)
print("Checking flightcondition usage:")
print("="*60)
############################################################

from aeroutils.flightcondition import FlightCondition
from aeroutils.units import *

# Compute flight conditions for a scalar or array of altitudes
altitudes = [0, 10e3, 33.5e3] * unit('ft')
fc = FlightCondition(altitudes, EAS=300*unit('knots'))
print(f"Airspeed in multiple formats: {fc}")
print(f"Even more data: {fc.tostring()}")
print(f"Access atmospheric data (see Atmosphere class): {fc.atm}")

# Or view fc formats individually:")
print(f"\nThe Mach number is {fc.mach:.5g}")
print(f"The true airspeed is {fc.TAS:.5g}")
print(f"The calibrated airspeed is {fc.CAS:.5g}")
print(f"The equivalent airspeed is {fc.EAS:.5g}")

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

############################################################
print("\n")
print("="*60)
print("Checking atmosphere usage:")
print("="*60)
############################################################
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
