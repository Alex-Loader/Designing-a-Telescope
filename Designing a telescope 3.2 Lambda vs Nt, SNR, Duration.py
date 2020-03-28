# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 10:36:05 2020

@author: Alex
"""

import numpy as np 
import matplotlib.pyplot as plt

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
pc=3.26 #1 parsec in light years
Atm_Mass = 1.66e-27 #Atomic Mass Unit

Rs=696.34e6 #Radius of the sun
Ms=1.989e30 #Mass of the sun
Re=6.371e6 #Radius of the Earth
Me=5.97e24 #Mass of the Earth


def Star_and_Detector (Lambda, T, R, rad, D, R_power): 
    #Lambda = Wavelength, T = Temperature of star, R=Radius of the star, 
    #rad = radius of mirror, D = distance to star (ly), R_power = Reolving power

    D = D * ly
    R = R * Rs
    
    Photon_Energy = h * c / Lambda

    B_lambda = ( (2*h*c**2) / (Lambda**5) ) * ( 1 / (np.exp( (h*c) / (Lambda*kb*T) ) - 1) ) #Planck function for given wavelength

    Power_lambda = B_lambda * np.pi * (R**2 / D**2) * 1e-9 *np.pi * rad**2 #Calculates the power per wavelength

    Photon_Count_lambda = Power_lambda / Photon_Energy #Photon count for given wavelength with delta lambda of 1nm
    
    Bin_size = Lambda * 1e9 / R_power #Finding delta lambda given R power
    
    Photon_Count = Photon_Count_lambda * Bin_size #Photon count given true delta lambda
    
    E_Count = Photon_Count * QE * Opt_eff * atm_eff #Electron count given efficiencies

    return (Photon_Count, E_Count, Bin_size)

"""Transit function"""
"""Calculates the period of rotation and the transit time of a planet, as well
as how much light this planet will block"""

def Transit (M, R, r, a):
    M = M*Ms
    R = R*Rs
    r = r*Re
    
    a_m = a * 1.5e11 #semi-major axis in metres 
    
    Period_s = 2*np.pi * np.sqrt(a_m**3 / (G * M))
    Period = Period_s / (365.25 * 24 * 60 * 60)
    Period_d = Period * 365.25   
    
    v = 2 * np.pi * a_m / Period_s
    v2 = np.sqrt(G * (M*Ms) / a_m)
    
    T14 = 2 * (R + r) / v #Time of transit (seconds)
    
    Blocked = (r**2 / R**2)
    
    return (T14, Period, Period_d, Period_s, Blocked, a_m, v, v2)

#This function returns the time of the transit, the total time of one orbit, the proportion of light blocked during the transit and the semi-major axis


"""Atmosphere function"""
"""Finds the scale height of the atmosphere and reduces the signal appropriately,
should result in our final signal"""


def Atmosphere(T, R, r, m, a_m, mu):    
    R = R*Rs
    r = r*Re
    m = m*Me
    mu = mu*Atm_Mass
    
    g=G*m/r**2 #calculating gravitational acceleration needed for scale height of the atmosphere

    T_atm =  T * np.sqrt(R / (2 * a_m)) 
    
    H = kb * T_atm / (g * mu)  #In metres
    A = (10 * r * H / R**2) #ppm
    
    
    return(T_atm, g, H,  A)


#Returns the atmospheric effective temperature of the planet, and the impact of the atmosphere on the readings

"""Noise function"""
"""Assuming only photon noise, calculates the relevant SNR of the observation"""

def Noise(T14, T_ex, Signal, Blocked, A, Target):
    #N = 1 #Assuming a small set number of exposures
    N = 2 * T14 / T_ex #Maximising the number of exposures
    
    True_Signal = Signal * T_ex
    
    Sigma_S = np.sqrt(True_Signal)
 
    Sigma_P = (np.sqrt(2) * (Sigma_S / np.sqrt(N/2))) / (True_Signal)

    Sigma_Ap = Sigma_P * np.sqrt(2)
    
    SNR =  A / (Sigma_Ap)
    
    n = (Target/SNR)**2
    
    return (N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n)




Period_Array = np.zeros(4)
A_Array = np.zeros(4)
T14_Array = np.zeros(4)
Ap_Array = np.zeros(4)


#Hypothetical varying space telescope
Lambda = np.linspace(500e-9,12000e-9,12000)
rad = 2    
R_power_100 = 40
QE = 0.8
Opt_eff = 0.5
atm_eff = 1

#Earth and Trappist (M8V)
#Star
T = 2559
M = 0.0802
R = 0.117
D = 10
#Planet
m = 1
r = 1
a = (R*Rs) / (2 * (280 / T)**2) / 1.5e11
mu = 28.97 

Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_100)
T14, Period, Period_d, Period_s, Blocked, a_m, v, v2 = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
plt.loglog(Lambda, n, label="M-class system")

Period_Array[0] = Period_d
A_Array[0] = a_m
T14_Array[0] = T14
Ap_Array[0] = A



#Earth and Kepler 62 (K2V)
#Star
T = 4925
M = 0.69
R = 0.64
D = 10
#Planet
m = 1
r = 1
a = (R*Rs) / (2 * (280 / T)**2) / 1.5e11
mu = 28.97 

Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_100)
T14, Period, Period_d, Period_s, Blocked, a_m, v, v2 = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
plt.loglog(Lambda, n, label="K-class system")

Period_Array[1] = Period_d
A_Array[1] = a_m
T14_Array[1] = T14
Ap_Array[1] = A


#Earth and Sol (G2V)
#Star
T = 5778
M = 1
R = 1
D = 10
#Planet
m = 1
r = 1
a = (R*Rs) / (2 * (280 / T)**2) / 1.5e11
mu = 28.97 

Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_100)
T14, Period, Period_d, Period_s, Blocked, a_m, v, v2 = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
plt.loglog(Lambda, n, label="G-class system")

Period_Array[2] = Period_d
A_Array[2] = a_m
T14_Array[2] = T14
Ap_Array[2] = A


#Earth and Tau Bootis (F7IV)
#Star
T = 6400
M = 1.34
R = 1.46
D = 10
#Planet
m = 1
r = 1
a = (R*Rs) / (2 * (280 / T)**2) / 1.5e11
mu = 28.97 

Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_100)
T14, Period, Period_d, Period_s, Blocked, a_m, v, v2 = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
plt.loglog(Lambda, n, label="F-class system")

Period_Array[3] = Period_d
A_Array[3] = a_m
T14_Array[3] = T14
Ap_Array[3] = A


plt.ylabel("Number of Transits")
plt.xlabel("Wavelength Observed (m)")

plt.xlim(500e-9,12000e-9)
plt.title("Number of transits required to reach SNR=5 on transit spectroscopy measurements (logged axes)")
plt.legend(loc="best")
txt="This simulation used an R power of 100, a mirror of radius 2m with QE=0.8, Optical Efficiency=0.5"
plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='center', fontsize=12)



"""~~~~~~~~~~~~~~~~~~~ Non-Logged Plot ~~~~~~~~~~~~~~~~~"""
plt.figure()
#Earth and Trappist (M8V)
#Star
T = 2559
M = 0.0802
R = 0.117
D = 10
#Planet
m = 1
r = 1
a = (R*Rs) / (2 * (280 / T)**2) / 1.5e11
mu = 28.97 

Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_100)
T14, Period, Period_d, Period_s, Blocked, a_m, v, v2 = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
plt.plot(Lambda, n, label="M-class system")


#Earth and Kepler 62 (K2V)
#Star
T = 4925
M = 0.69
R = 0.64
D = 10
#Planet
m = 1
r = 1
a = (R*Rs) / (2 * (280 / T)**2) / 1.5e11
mu = 28.97 

Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_100)
T14, Period, Period_d, Period_s, Blocked, a_m, v, v2 = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
plt.plot(Lambda, n, label="K-class system")




#Earth and Sol (G2V)
#Star
T = 5778
M = 1
R = 1
D = 10
#Planet
m = 1
r = 1
a = (R*Rs) / (2 * (280 / T)**2) / 1.5e11
mu = 28.97 

Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_100)
T14, Period, Period_d, Period_s, Blocked, a_m, v, v2 = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
plt.plot(Lambda, n, label="G-class system")



#Earth and Tau Bootis (F7IV)
#Star
T = 6400
M = 1.34
R = 1.46
D = 10
#Planet
m = 1
r = 1
a = (R*Rs) / (2 * (280 / T)**2) / 1.5e11
mu = 28.97 

Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_100)
T14, Period, Period_d, Period_s, Blocked, a_m, v, v2 = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
plt.plot(Lambda, n, label="F-class system")

plt.ylabel("Number of Transits")
plt.xlabel("Wavelength Observed (m)")

plt.xlim(500e-9,12000e-9)
plt.title("Number of transits required to reach SNR=5 on transit spectroscopy measurements")
plt.legend(loc="best")
txt="This simulation used an R power of 100, a mirror of radius 2m with QE=0.8, Optical Efficiency=0.5"
plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='center', fontsize=12)

"""~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sigma P Plot for each star ~~~~~~~~~~~~~~~~~~~~~~~~~~"""
plt.figure()

#Earth and Trappist (M8V)
#Star
T = 2559
M = 0.0802
R = 0.117
D = 10
#Planet
m = 1
r = 1
a = (R*Rs) / (2 * (280 / T)**2) / 1.5e11
mu = 28.97 

Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_100)
T14, Period, Period_d, Period_s, Blocked, a_m, v, v2 = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
plt.loglog(Lambda, Sigma_P, label="M-class system")

Period_Array[0] = Period_d
A_Array[0] = a_m
T14_Array[0] = T14
Ap_Array[0] = A



#Earth and Kepler 62 (K2V)
#Star
T = 4925
M = 0.69
R = 0.64
D = 10
#Planet
m = 1
r = 1
a = (R*Rs) / (2 * (280 / T)**2) / 1.5e11
mu = 28.97 

Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_100)
T14, Period, Period_d, Period_s, Blocked, a_m, v, v2 = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
plt.loglog(Lambda, Sigma_P, label="K-class system")

Period_Array[1] = Period_d
A_Array[1] = a_m
T14_Array[1] = T14
Ap_Array[1] = A


#Earth and Sol (G2V)
#Star
T = 5778
M = 1
R = 1
D = 10
#Planet
m = 1
r = 1
a = (R*Rs) / (2 * (280 / T)**2) / 1.5e11
mu = 28.97 

Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_100)
T14, Period, Period_d, Period_s, Blocked, a_m, v, v2 = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
plt.loglog(Lambda, Sigma_P, label="G-class system")

Period_Array[2] = Period_d
A_Array[2] = a_m
T14_Array[2] = T14
Ap_Array[2] = A


#Earth and Tau Bootis (F7IV)
#Star
T = 6400
M = 1.34
R = 1.46
D = 10
#Planet
m = 1
r = 1
a = (R*Rs) / (2 * (280 / T)**2) / 1.5e11
mu = 28.97 

Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_100)
T14, Period, Period_d, Period_s, Blocked, a_m, v, v2 = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
plt.loglog(Lambda, Sigma_P, label="F-class system")

Period_Array[3] = Period_d
A_Array[3] = a_m
T14_Array[3] = T14
Ap_Array[3] = A


plt.ylabel("Sigma P")
plt.xlabel("Wavelength Observed (m)")

plt.xlim(500e-9,12000e-9)
plt.title("Sigma P of Exoplanet Transit")
plt.legend(loc="best")
txt="This simulation used an R power of 100, a mirror of radius 2m with QE=0.8, Optical Efficiency=0.5"
plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='center', fontsize=12)


"""~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SNR for each star ~~~~~~~~~~~~~~~~~~~~~~~~~~"""
plt.figure()

#Earth and Trappist (M8V)
#Star
T = 2559
M = 0.0802
R = 0.117
D = 10
#Planet
m = 1
r = 1
a = (R*Rs) / (2 * (280 / T)**2) / 1.5e11
mu = 28.97 

Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_100)
T14, Period, Period_d, Period_s, Blocked, a_m, v, v2 = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
plt.loglog(Lambda, SNR, label="M-class system")

Period_Array[0] = Period_d
A_Array[0] = a_m
T14_Array[0] = T14
Ap_Array[0] = A



#Earth and Kepler 62 (K2V)
#Star
T = 4925
M = 0.69
R = 0.64
D = 10
#Planet
m = 1
r = 1
a = (R*Rs) / (2 * (280 / T)**2) / 1.5e11
mu = 28.97 

Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_100)
T14, Period, Period_d, Period_s, Blocked, a_m, v, v2 = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
plt.loglog(Lambda, SNR, label="K-class system")

Period_Array[1] = Period_d
A_Array[1] = a_m
T14_Array[1] = T14
Ap_Array[1] = A


#Earth and Sol (G2V)
#Star
T = 5778
M = 1
R = 1
D = 10
#Planet
m = 1
r = 1
a = (R*Rs) / (2 * (280 / T)**2) / 1.5e11
mu = 28.97 

Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_100)
T14, Period, Period_d, Period_s, Blocked, a_m, v, v2 = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
plt.loglog(Lambda, SNR, label="G-class system")

Period_Array[2] = Period_d
A_Array[2] = a_m
T14_Array[2] = T14
Ap_Array[2] = A


#Earth and Tau Bootis (F7IV)
#Star
T = 6400
M = 1.34
R = 1.46
D = 10
#Planet
m = 1
r = 1
a = (R*Rs) / (2 * (280 / T)**2) / 1.5e11
mu = 28.97 

Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_100)
T14, Period, Period_d, Period_s, Blocked, a_m, v, v2 = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
plt.loglog(Lambda, SNR, label="F-class system")

Period_Array[3] = Period_d
A_Array[3] = a_m
T14_Array[3] = T14
Ap_Array[3] = A


plt.ylabel("SNR")
plt.xlabel("Wavelength Observed (m)")

plt.xlim(500e-9,12000e-9)
plt.title("SNR of Atmpsphere Spectroscopy")
plt.legend(loc="best")
txt="This simulation used an R power of 100, a mirror of radius 2m with QE=0.8, Optical Efficiency=0.5"
plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='center', fontsize=12)


"""~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Mission Duration for each star ~~~~~~~~~~~~~~~~~~~~~~~~~~"""
plt.figure()

#Earth and Trappist (M8V)
#Star
T = 2559
M = 0.0802
R = 0.117
D = 10
#Planet
m = 1
r = 1
a = (R*Rs) / (2 * (280 / T)**2) / 1.5e11
mu = 28.97 

Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_100)
T14, Period, Period_d, Period_s, Blocked, a_m, v, v2 = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
plt.loglog(Lambda, Time_Required, label="M-class system")

Period_Array[0] = Period_d
A_Array[0] = a_m
T14_Array[0] = T14
Ap_Array[0] = A



#Earth and Kepler 62 (K2V)
#Star
T = 4925
M = 0.69
R = 0.64
D = 10
#Planet
m = 1
r = 1
a = (R*Rs) / (2 * (280 / T)**2) / 1.5e11
mu = 28.97 

Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_100)
T14, Period, Period_d, Period_s, Blocked, a_m, v, v2 = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
plt.loglog(Lambda, Time_Required, label="K-class system")

Period_Array[1] = Period_d
A_Array[1] = a_m
T14_Array[1] = T14
Ap_Array[1] = A


#Earth and Sol (G2V)
#Star
T = 5778
M = 1
R = 1
D = 10
#Planet
m = 1
r = 1
a = (R*Rs) / (2 * (280 / T)**2) / 1.5e11
mu = 28.97 

Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_100)
T14, Period, Period_d, Period_s, Blocked, a_m, v, v2 = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
plt.loglog(Lambda, Time_Required, label="G-class system")

Period_Array[2] = Period_d
A_Array[2] = a_m
T14_Array[2] = T14
Ap_Array[2] = A


#Earth and Tau Bootis (F7IV)
#Star
T = 6400
M = 1.34
R = 1.46
D = 10
#Planet
m = 1
r = 1
a = (R*Rs) / (2 * (280 / T)**2) / 1.5e11
mu = 28.97 

Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_100)
T14, Period, Period_d, Period_s, Blocked, a_m, v, v2 = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
plt.loglog(Lambda, Time_Required, label="F-class system")

Period_Array[3] = Period_d
A_Array[3] = a_m
T14_Array[3] = T14
Ap_Array[3] = A


plt.ylabel("Mission Duration (years)")
plt.xlabel("Wavelength Observed (m)")

plt.xlim(500e-9,12000e-9)
plt.title("Mission Duration for Atmosphere Spectroscopy of Exoplanet Systems")
plt.legend(loc="best")
txt="This simulation used an R power of 100, a mirror of radius 2m with QE=0.8, Optical Efficiency=0.5"
plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='center', fontsize=12)
