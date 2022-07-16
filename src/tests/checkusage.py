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

# Print all flight condition data:
print(f"{fc}")

# Print while specifying abbreviated output in SI units:
print(f"\n{fc.tostring(full_output=False, US_units=False)}")

# Access flight speed formats individually
M = fc.speed.M
U_inf, U_CAS, U_EAS = fc.speed.TAS, fc.speed.CAS, fc.speed.EAS

# Access atmospheric data alone (see Atmosphere class for options)
alt = fc.altitude  # access Atmosphere object 'altitude'
p, T, rho, nu, a = alt.p, alt.T, alt.rho, alt.nu, alt.a

# Input true/calibrated/equivalent airspeed or Mach number
fc_TAS = FlightCondition(altitudes, TAS=300*unit('knots'))
fc_CAS = FlightCondition(altitudes, CAS=300*unit('knots'))
fc_EAS = FlightCondition(altitudes, EAS=300*unit('knots'))
fc_M = FlightCondition(altitudes, M=0.4535*dimless)

# Specify desired units on input and output
altitudes_in_km = [0, 3.048, 10.2108] * unit('km')
lengthscale = 60 * unit('in')  # arbitrary length scale of interest
fc_other_units = FlightCondition(altitudes, EAS=154.33*unit('m/s'),
                                 ell=lengthscale)
U_TAS = fc_other_units.speed.TAS
print(f"\nThe true airspeed in m/s is {U_TAS.to('m/s'):.5g}")
print(f"The true airspeed in km/s is {U_TAS.to('km/s'):.5g}")

# Compute additional derived quantities (see class for all options)
print(f"\nThe dynamic pressure in psi is "
      f"{fc.speed.q_inf.to('psi'):.5g}")
print(f"The Reynolds number is {fc.length.Re:.5g}")
print(f"The Reynolds number per-unit-length [1/in] is "
      f"{fc.length.Re_by_ell.to('1/in'):.5g}")

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

# Print all atmospheric quantities:
print(f"\n{atm}")

# Print while specifying abbreviated output in US units:
print(f"\n{atm.tostring(full_output=False, US_units=True)}")

# See also the linspace() function from numpy, e.g.
# h = linspace(0, 81.0, 82) * unit('km')

# Access individual properties and convert to desired units: "
p, T, rho, nu, a = atm.p, atm.T, atm.rho, atm.nu, atm.a
print(f"\nThe pressure in psi is {p.to('psi'):.5g}")

# Compute additional properties such as thermal conductivity,
# mean free path, and more (see class for all options)
print(f"\nThe thermal conductivity is {atm.k:.5g}"
    f"\nThe mean free path = {atm.MFP:.5g}")
