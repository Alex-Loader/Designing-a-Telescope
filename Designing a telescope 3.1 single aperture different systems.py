
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 12:05:31 2020

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





"""Defining the Star and Planet"""
"""
#WASP-19b
#Star
T = 5568
M = 0.904
R = 1.004
D = 881.4
#Planet
m = 1.139*317.8
r = 1.41*11.209
a = 0.016 #In AU
mu = 2.22 #Value based on Jupiter


#WASP-19b
#Star
T = 5568
M = 0.904
R = 1.004
D = 881.4
#Planet
m = 1.139*317.8
r = 1.41*11.209
a = 0.016 #In AU
mu = 2.22 #Value based on Jupiter

#WASP-96b
#Star
T = 5540
M = 1.06
R = 1.05
D = 1161.12
#Planet
m = 0.48*317.8
r = 1.2*11.209
a = 0.0453 #In AU
mu = 2.22 #Value based on Jupiter


#GJ 1214b
#Star
T = 3026
M = 0.15
R = 0.216
D = 47.8
#Planet
m = 6.43
r = 2.6
a = 0.014 #In AU
mu = 18 


#HD 219134b
#Star
T = 4699
M = 0.81
R = 0.778
D = 20.35
#Planet
m = 4.74
r = 1.602
a = 0.03876 #In AU
mu = 18

"""
#Earth-like w Trappist
#Star
T = 2511
M = 0.089
R = 0.121
D = np.linspace(0,50,100)
#Planet
m = 1
r = 1
a = (R*Rs) / (2 * (280 / T)**2) / 1.5e11
mu = 28.97 



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
    
     
    Period = np.sqrt(a**3) #Gives the period of the planet in years
    
    a_m = a * 1.5e11 #semi-major axis in metres
    Period_s = Period * 365.25 * 24 * 60 * 60 #Period is now in seconds
    Period_d = Period * 365.25
    
    v = 2 * np.pi * a_m / Period_s
    
    T14 = 2 * (R + r) / v #Time of transit (seconds)
    
    Blocked = (r**2 / R**2)
    
    return (T14, Period, Period_d, Period_s, Blocked, a_m)

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
    N = T14 / T_ex #Maximising the number of exposures
    True_Signal = Signal * T_ex
    
    Sigma_S = np.sqrt(True_Signal)
 
    Sigma_P = (np.sqrt(2) * (Sigma_S / np.sqrt(N/2))) / (True_Signal) #Is this right?

    
    Sigma_Ap = np.sqrt(Sigma_P) * np.sqrt(2)
    
    SNR = np.sqrt(N) * A / Sigma_Ap #Not true_signal, need to focus on smaller signal! Signal is not total signal but is P
    
    n = (Target/SNR)**2
    
    return (N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n)






#Hypothetical varying space telescope
Lambda = 10000e-9
rad = np.linspace(0,20,200)    
R_power_100 = 100
R_power_40 = 40
QE = 0.8
Opt_eff = 0.5
atm_eff = 1

#Setting up multiple subplots
fig, axs = plt.subplots(4, sharex = True)
fig.suptitle("Time Required for a Detection of Atmosphere for Various Stellar Classifications")
fig.text(0.5, 0.04, "Aperture radius (m)", ha='center')
fig.text(0.04, 0.5, "Time required to reach SNR=5", va='center', rotation='vertical')


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
T14, Period, Period_d, Period_s, Blocked, a_m = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
axs[0].semilogy(rad, Time_Required, label="R=100")
Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_40)
T14, Period, Period_d, Period_s, Blocked, a_m = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
axs[0].semilogy(rad, Time_Required, label="R=40")
axs[0].legend(loc="best")

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
T14, Period, Period_d, Period_s, Blocked, a_m = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
axs[1].semilogy(rad, Time_Required, label="R=100")
Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_40)
T14, Period, Period_d, Period_s, Blocked, a_m = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
axs[1].semilogy(rad, Time_Required, label="R=40")
plt.ylim(0,30)



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
T14, Period, Period_d, Period_s, Blocked, a_m = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
axs[2].semilogy(rad, Time_Required, label="R=100")
Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_40)
T14, Period, Period_d, Period_s, Blocked, a_m = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
axs[2].semilogy(rad, Time_Required, label="R=40")
plt.ylim(0,20)
plt.legend(loc="best")



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
T14, Period, Period_d, Period_s, Blocked, a_m = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
axs[3].semilogy(rad, Time_Required, label="R=100")
Photon_Count, Signal, Bin_size = Star_and_Detector(Lambda, T, R, rad, D, R_power_40)
T14, Period, Period_d, Period_s, Blocked, a_m = Transit(M, R, r, a) 
T_atm, g, H, A = Atmosphere(T, R, r, m, a_m, mu)
N, True_Signal, Sigma_S, Sigma_P, Sigma_Ap, SNR, n = Noise(T14, 10, Signal, Blocked, A, 5)
Time_Required = n * Period
axs[3].semilogy(rad, Time_Required, label="R=40")
plt.ylim(0,10000)
