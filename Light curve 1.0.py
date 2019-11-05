# -*- coding: utf-8 -*
"""
Created on Thu Oct 31 13:42:29 2019

@author: Alex
"""
import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


z=np.linspace(-2,2,100) #centre-to-centre distance as a ratio of the stars radius
R=7e8 #radius of a star
r=2e8 #radius of a planet
p=r/R #ratio of radii
Z=abs(z)
light=np.ones_like(z)

for i in range(len(z)):
#    if 1+p<abs(z[i]): #The planet is not obstructing
#        light[i]=1
    if 1+p>=abs(z[i]) and 1-p<abs(z[i]): #The planet is moving in front of the star, but not in front yet
        light[i]=1-(1/np.pi)*((p**2*np.arccos((p**2+Z[i]**2-1)/(2*p*Z[i])))+np.arccos((1-p**2+Z[i]**2)/2*Z[i])-np.sqrt((4*Z[i]**2-(1-p**2+Z[i]**2)**2)/4))
    if abs(z[i])<=1-p: #The planet is directly infront of the star
        light[i]=1-p**2
    if abs(z[i])<=p-1: #The planet is larger than the star
        light[i]=0

plt.plot(z,light)