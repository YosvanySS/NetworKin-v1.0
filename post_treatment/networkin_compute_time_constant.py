#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 14:53:06 2025

@author: Julien Denis
@mail:   julien.denis@univ-amu.fr

@author: Yosvany Silva-Solís
@mail:   yosvany.silva-solis@univ-amu.fr

Project: NetworKin
file:    networkin_compute_time_constant.py
object:  script to compute the time constant from NetworKin simulation

"""

####  reset
from IPython import get_ipython;
get_ipython().run_line_magic('reset', '-sf')


#### Libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import h5py


#### Script inputs

    ########## INPUTS ##########
    
# NetworKin output filename
networkin_filename = './networkin_W_Cu_Dirichlet_300K' # without .h5 extension
name_simu = 'W/Cu full interface'

# site label used for the time constant calculation (string)
site_time_constant = 'c'

# steady-state time [s]
time_steady_state = 1e6

# percentage defining the time constant [%]
k_time_constant = 99.

# figure number
ifig = 1

    ############################


#### Read hdf5 result file from NetworKin

f = h5py.File(networkin_filename+'.h5','r')

network_label = [x.decode() for x in f["network_label"]]

time = np.array(f["time"])

Ci = np.array(f["Ci"])

stored_value = f['is_Ci_thermo_eq'][()]
is_Ci_thermo_eq = bool(stored_value)

if is_Ci_thermo_eq:
    Ci_thermo_eq = np.array(f['Ci_thermo_eq'])

f.close()


#### Time to reach k_time_constant% of steady-state concentration at site "site_time_constant"

    # We find the index of site "site_time_constant"
idx_site = network_label.index(site_time_constant)

    # We compute the steady-state concentration at site "site_time_constant"
spl_Ci = CubicSpline(time, Ci[idx_site,:])
Ci_steady_state = spl_Ci(time_steady_state)

    # We compute the time to reach k_time_constant% of steady-state
Ci_k_percent = k_time_constant/100.*Ci_steady_state
Ci_root = Ci[idx_site,:] - Ci_k_percent
spl_root = CubicSpline(time, Ci_root)

time_k_percent = spl_root.roots(extrapolate=False)[0]

print("Time to reach "+str(k_time_constant)+"% of steady-state:")
print("    - "+name_simu+" : ", time_k_percent)


#### We plot the result

ifig = ifig+1
fig = plt.figure(ifig)

p = plt.plot(time,Ci[idx_site],label=network_label[idx_site]+' - '+name_simu,linestyle='-')
        
plt.hlines(Ci_k_percent,time[0],time[-1],linestyles='dashed',color='k')
    
bottom, top = plt.ylim()
plt.vlines(time_k_percent,bottom,top,linestyles='dashed',color='k')
    
plt.grid()

plt.xlabel('time [s]')
plt.ylabel('$\mathrm{C_i}$ [$\mathrm{m^{-2}}$]')

plt.xscale("log")
plt.yscale("log")

plt.legend()