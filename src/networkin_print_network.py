#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 09:41:34 2025

@author: Julien Denis
@mail:   julien.denis@univ-amu.fr

@author: Yosvany Silva-Solís
@mail:   yosvany.silva-solis@univ-amu.fr

Project: NetworKin
file:    networkin_print_network.py
object:  function to print the network created by NetworKin

"""

def networkin_print_network(network,t0):
    
    # network: network object
    # t0     : normalisation time [s]
    
    import numpy as np

    for i in range(0,np.size(network)):
        
        if i > 0:
            print('----')

        print('Site label: ' , network[i].site_label)
        
        # loop on jumps
        for j in range(network[i].Njump):
            print('    Jump number '+str(j)+':')
            
            print('        Neighbor site     : ', \
                  network[i].jump[j].neighbor_site_label)
            print('        Forward  frequency: ', \
                  "{:e}".format(network[i].jump[j].forward_freq/t0))
            print('        Backward frequency: ', \
                  "{:e}".format(network[i].jump[j].backward_freq/t0))

            print('    ')