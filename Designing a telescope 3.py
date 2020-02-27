
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 12:05:31 2020

@author: Alex
"""

import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.close("all")

"""Star and Instrument function"""
"""Needs to calculate the power of the star reaching our detector, and then the
amount of electrons created from the photon count"""

#Defining important constants:
c=3e8 #speed of light
h=6.626e-34 #Planck constant
kb=1.38e-23 #Boltzmann constant
G=6.67e-11 #Gravitational Constant
sigma=5.67e-8 #Stefan-Boltzmann constant
ly=9.46e15 #1 Light Year in metres
Rs=6.95e8 #Radius of the sun
Ms=1.989e30 #Mass of the sun
Re=6.371e6 #Radius of the Earth
Me=5.97e24 #Mass of the Earth
QE = 0.6 #Quantum Efficiency
Opt_eff = 0.6 #Optical Efficiency

def Star_and_Detector (Lambda, T, R, rad, D, R_power): 
    #Lambda = Wavelength, T = Temperature of star, R=Radius of the star, 
    #rad = radius of mirror, D = distance to star (ly), R_power = Reolving power

    D = D * ly
    R = R * Ms
    
    Photon_Energy = h * c / Lambda

    B_lambda = ( (2*h*c**2) / (Lambda**5) ) * ( 1 / (np.exp( (h*c) / (Lambda*kb*T) ) - 1) ) #Planck function for given wavelength

    Power_lambda = B_lambda * np.pi * (R**2 / D**2) * 1e-9 *np.pi * rad**2 #Calculates the power per wavelength

    Photon_Count_lambda = Power_lambda / Photon_Energy #Photon count for given wavelength with delta lambda of 1nm
    
    Bin_size = Lambda * 1e9 / R_power #Finding delta lambda given R power
    
    Photon_Count = Photon_Count_lambda * Bin_size #Photon count given true delta lambda
    
    E_Count = Photon_Count * QE * Opt_eff #Electron count given efficiencies

    return (Photon_Count, E_Count)

Photon_Count, Signal = Star_and_Detector(1000e-9, 3318, 0.42, 4, 31.8, 30)

"""Transit function"""
"""Calculates the period of rotation and the transit time of a planet, as well
as how much light this planet will block"""

def Transit (M, R, r):
    M = M*Ms
    R = R*Rs
    r = r*Re
    
    Lum_sun = 3.828e26 #Solar luminosity in Watts
    M_sun = 1.989e30 #Mass of sun in kg
    Lum_star = Lum_sun * (M / M_sun)**3.5
    
    a = np.sqrt(Lum_star / Lum_sun) #semi-major axis in AU 
    Period = np.sqrt(a**3) #Gives the period of the planet in years
    
    a_m = a * 1.5e11 #semi-major axis in metres
    Period_s = Period * 365.25 * 24 * 60 * 60 #Period is now in seconds
    
    v = 2 * np.pi * a_m / Period_s
    T14 = 2 * (R + r) / v #Time of transit (seconds)
    
    Blocked = (r**2 / R**2)
    
    return (T14, Period, Period_s, Blocked, a_m)

#This function returns the time of the transit, the total time of one orbit, the proportion of light blocked during the transit and the semi-major axis

T14, Period, Period_s, Blocked, a_m = Transit(0.41, 0.42, 4.327) 



"""Atmosphere function"""
"""Finds the scale height of the atmosphere and reduces the signal appropriately,
should result in our final signal"""


def Atmosphere(T, R, r, m, a_m, mu):    
    R = R*Rs
    r = r*Re
    m = m*Me
    
    g=G*m/r**2 #calculating gravitational acceleration needed for scale height of the atmosphere

    T_atm = T * np.sqrt(R / (2 * a_m))
    
    H = kb * T_atm / (g * mu) 
    A = (10 * r * H / R**2)
    
    return(T_atm, A)

T_atm, A = Atmosphere(3318, 0.42, 4.327, 21.36, a_m, 4.65e-26)

#Returns the atmospheric effective temperature of the planet, and the impact of the atmosphere on the readings

"""Noise function"""
"""Assuming only photon noise, calculates the relevant SNR of the observation"""

def Noise(T14, T_ex, Signal, Blocked, A):
    N = 20 #Assuming a small set number of exposures
    #N = T14 / T_ex #Maximising the number of exposures
    True_Signal = Signal * T_ex
    Atm_Signal = A #* True_Signal - need it as a proportion of the true signal
    
    Sigma_S = np.sqrt(True_Signal)
 
    Sigma_P = (np.sqrt(2) * (Sigma_S / np.sqrt(N/2))) / True_Signal #Is this right?
    
    Sigma_Ap = np.sqrt(Sigma_P) * np.sqrt(2)
    
    SNR = A / Sigma_Ap #Not true_signal, need to focus on smaller signal! Signal is not total signal but is P
    
    
    return (N, True_Signal, Atm_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR)


N, True_Signal, Atm_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR = Noise(T14, 1, Signal, Blocked, A)
