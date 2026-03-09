#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 15:32:03 2025

@author: Julien Denis
@mail:   julien.denis@univ-amu.fr

@author: Yosvany Silva-Solís
@mail:   yosvany.silva-solis@univ-amu.fr

Project: NetworKin
file:    networkin_frequency.py
object:  function to calculate Arrhenius-type frequency

"""

def networkin_frequency(WaZ,m,k0,E_a,Temp):
    
    # WaZ : logical to use Wert and Zener for pre-exponential factor
    # m   : mass of the specie
    # k0  : if WaZ = False -> pre-exponential factor [s-1]
    #       if WaZ = True  -> distance between two sites [Angström]
    # E_a : activation energy [eV]
    # Temp: temperature [K]
    
    import numpy as np
    
    from networkin_physical_constants import e_charge, k_B
            
    if WaZ == False:
        nu_0 = k0
    elif WaZ == True:
        nu_0 = np.sqrt(E_a*e_charge/(2.*m))/(k0*1.e-10)

    freq = nu_0*np.exp(-E_a/(k_B*Temp))
    
    return nu_0, freq