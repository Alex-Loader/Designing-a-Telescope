# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 17:15:12 2019

@author: Alex
"""
import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


#Defining the black body curve of a star
Lambda=(np.linspace(0,15000, 151))*1e-9 #using nanometres for the wavelengths
c=3e8
h=6.626e-34
kb=1.38e-23
T=6000
R=7e8 #radius of star
G=6.67e-11
Radiance=( (2*h*c**2) / (Lambda**5) ) * ( (1) / (np.exp( (h*c) / (Lambda*kb*T) ) - 1) ) #Using Plancks law and a wavelength dependence




Lower_wavelength=300e-9
Upper_wavelength=1000e-9
Lower_bound=Upper_bound=0


Waveband=Lambda[15:50]
#This isn't power! Need to account for solid angles/area
Luminosity=Lambda*4*np.pi*R**2
#This has given the spectral density for this waveband - now it has to reach Earth


sigma=5.67e-8
L0=Radiance * 4*np.pi*R**2
#L=4*np.pi*R**2*sigma*T**4 #Luminosity emitted by star

r=6e6 #radius of earth like planet
m_planet=5.97e24 #mass of earth like planet
g=G*m_planet/r**2 #calculating gravitational acceleration needed for scale height of the atmosphere
mu=4.65e-26 #Mean mass of molecules in atmosphere - assuming predominantly nitrogen similar to Earth
T_atm=250 #atmospheric temperature

H=kb*T_atm/(g*mu)
A=(10*r*H/R**2)
L1=L0*(1-(r**2/R**2))*1-((10*r*H/R**2)) #Proportion of luminosity blocked by transiting planet

E_atm=0.8 #atmospheric efficiency
L2=L1*E_atm #some light absorbed by atmosphere

D=4.73e18 #Distance from star to Earth
d=4 #radius of mirror
E_inst=0.2 #instrument efficiency

P1=(d**2*L2*E_inst)/(4*D**2)

E_det=0.5 #quantum efficieny of detector
E_ap=0.5 #apetrure efficieny due to diffraction patterns

P2=P1*E_det*E_ap

#assume a wavelength for now of 500nm
Photon_energy=h*c/Lambda

Photon_count=P2/Photon_energy

plt.figure()
plt.plot(Lambda, Radiance, label="Radiance")
plt.plot(Lambda, P2, label="Power Measured")
plt.plot(Lambda, Photon_count, label="Photon count from detector")
plt.yscale("log")
plt.legend(loc="best")

plt.figure()
plt.loglog(Lambda, Radiance, label="Radiance")
plt.loglog(Lambda, P2, label="Power Measured")
plt.loglog(Lambda, Photon_count, label="Photon count from detector")
plt.legend(loc="best")

plt.figure()
plt.plot(Lambda, Radiance)
plt.plot(Lambda, P2)
plt.plot(Lambda, Photon_count)
