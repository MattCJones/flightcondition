#!/usr/bin/env python
# Check that "Usage" commands run properly.

# flake8: noqa F401

############################################################
print("="*60)
print("Checking flightcondition usage:")
print("="*60)
############################################################

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

############################################################
print("")
print("="*60)
print("Checking atmosphere usage:")
print("="*60)
############################################################

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
      f"\nThe mean free path = {atm.MFP:.5g}")
