#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 14:53:06 2025

@author: Julien Denis
@mail:   julien.denis@univ-amu.fr

@author: Yosvany Silva-Solís
@mail:   yosvany.silva-solis@univ-amu.fr

Project: NetworKin
file:    networkin_plot.py
object:  script to plot results from NetworKin simulation

"""

####  reset
from IPython import get_ipython;
get_ipython().run_line_magic('reset', '-sf')


#### Libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py


#### INPUTS
    # NetworKin output filename 
#networkin_filename = 'networkin_simu'   # without .h5 extension
networkin_filename = 'networkin_W_Cu_scenario3_347K'
#networkin_filename = 'network_1_path' # without .h5 extension

    # figure number
ifig = 1


#### Read hdf5 result file
f = h5py.File(networkin_filename+'.h5','r')

network_label = [x.decode() for x in f['network_label']]

time = np.array(f['time'])

Ci = np.array(f['Ci'])

stored_value = f['is_Ci_thermo_eq'][()]
is_Ci_thermo_eq = bool(stored_value)

if is_Ci_thermo_eq:
    Ci_thermo_eq = np.array(f['Ci_thermo_eq'])

f.close()


#### We plot the results

fig1 = plt.figure(ifig)

for i in range(len(network_label)):
    p = plt.plot(time,Ci[i],label=network_label[i])
    color_plot = p[0].get_color()
    
    if is_Ci_thermo_eq:
        plt.hlines(Ci_thermo_eq[i],time[0],time[-1], \
                   label=network_label[i]+' thermo eq', \
                   linestyles='dashed',color=color_plot)
    
#    plt.plot(time[1:-1],Ci[:,i],label=network_label[i])

plt.grid()

plt.xlabel('time [s]')
plt.ylabel('$\mathrm{C_i}$ [$\mathrm{m^{-2}}$]')

plt.xscale('log')
plt.yscale('log')

plt.legend()