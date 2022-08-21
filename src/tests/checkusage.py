#!/usr/bin/env python
# Check that "Usage" commands run properly.

# flake8: noqa F401

#######################################################################
maxchars = 71  # 80 - 1 - 4*2 = 71 columns available
print("="*maxchars)
print("Checking flightcondition usage:")
print("="*maxchars)
#######################################################################

from flightcondition import FlightCondition, unit, dimless

# Compute flight condition at 3 km, Mach 0.5
fc = FlightCondition(3*unit('km'), M=0.5)

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
fc = FlightCondition([0, 9.8425, 20]*unit('kft'),
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

#######################################################################
print("")
print("="*maxchars)
print("Checking atmosphere usage:")
print("="*maxchars)
#######################################################################

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
