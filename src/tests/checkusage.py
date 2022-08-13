#!/usr/bin/env python
# Check that "Usage" commands run properly.

# flake8: noqa F401

############################################################
print("="*60)
print("Checking flightcondition usage:")
print("="*60)
############################################################

from flightcondition import FlightCondition, unit, dimless

# Compute flight condition at 3 km, Mach 0.5
fc = FlightCondition(3*unit('km'), M=0.5)

# Uncomment to print summary of flight condition quantities:
#print(f"{fc}")

# Uncomment to print abbreviated output in US units:
#print(f"\n{fc.tostring(full_output=False, unit_system="US")}")

# Access true, calibrated, equivalent airspeeds
KTAS = fc.byvel.TAS.to('knots')
KCAS = fc.byvel.CAS.to('knots')
KEAS = fc.byvel.EAS.to('knots')
print(f"Flying at {KTAS.magnitude:.4g} KTAS,"
      f" which is {KCAS.magnitude:.4g} KCAS,"
      f" or {KEAS.magnitude:.4g} KEAS")
# >>> Flying at 319.4 KTAS, which is 277.7 KCAS, or 275.1 KEAS

# Access atmospheric data (see Atmosphere class for more)
atm = fc.byalt  # access Atmosphere object
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
      f"{fc.byvel.q_inf.to('psi'):.3g}")
# >>> The dynamic pressure in psi is [1.78 1.78 1.78] psi
print(f"The Reynolds number is {fc.bylen.Re:.3g}")
# >>> The Reynolds number is [9.69e+07 8.82e+07 7.95e+07]

# Alternatively access quantities by their full name
print(fc.byvel.TAS == fc.byvel.byname.true_airspeed)
# >>> [ True  True  True]


############################################################
print("")
print("="*60)
print("Checking atmosphere usage:")
print("="*60)
############################################################

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
