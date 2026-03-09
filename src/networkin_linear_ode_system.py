#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 17:39:51 2025

@author: Julien Denis
@mail:   julien.denis@univ-amu.fr

@author: Yosvany Silva-Solís
@mail:   yosvany.silva-solis@univ-amu.fr

Project: NetworKin
file:    networkin_linear_ode_system.py
object:  function defining the network linear ODE system to solve
    
"""
def networkin_linear_ode_system(t,y,network,network_BCs):
    
    # t          : time [s]
    # y, Ci      : concentration network site [m-3]
    # network    : network object
    # network_BCs: network boundary conditions
    
    import sys
    import numpy as np
    
    
    Ci = y
    
    dCidt = np.zeros(Ci.shape)
    
    
    # Network input boundary
    if network_BCs[0].BC_type == 'Dirichlet':
        y[0]  = network_BCs[0].BC_val
        Ci[0] = network_BCs[0].BC_val
        
        dCidt[0] = 0.
        
    elif network_BCs[0].BC_type == 'Neumann 1':
        dCidt[0] = network_BCs[0].BC_val
        
        # loop on jumps
        for j in range(network[0].Njump):
            dCidt[0] =   dCidt[0] \
                       - network[0].jump[j].forward_freq \
                         *Ci[0] \
                       + network[0].jump[j].backward_freq \
                         *Ci[network[0].jump[j].neighbor_site_idx]
    
    elif network_BCs[0].BC_type == 'Neumann 2':
        dCidt[0] = - network_BCs[0].diff_coeff*(Ci[0]-network_BCs[0].BC_val) \
                                              /network_BCs[0].diff_length
        
        # loop on jumps
        for j in range(network[0].Njump):
            dCidt[0] =   dCidt[0] \
                       - network[0].jump[j].forward_freq \
                         *Ci[0] \
                       + network[0].jump[j].backward_freq \
                         *Ci[network[0].jump[j].neighbor_site_idx]
    
    else:
        print('Boundary condition input_BC_type set to incorrect value: ', \
              network_BCs[0].BC_type)
        print('Please use "Dirichlet", "Neumann 1" or , \
              "Neumann 2" as boundary condition type.')
        print('Exiting...')
        sys.exit()
    
    
    # Network core
    for i in range(1,Ci.size-1,1):
        
        # loop on jumps
        for j in range(network[i].Njump):
            dCidt[i] =  dCidt[i] \
                      - network[i].jump[j].forward_freq \
                        *Ci[i] \
                      + network[i].jump[j].backward_freq \
                        *Ci[network[i].jump[j].neighbor_site_idx]
    
    
    # Network output boundary
    if network_BCs[1].BC_type == 'Dirichlet':
        y[-1]  = network_BCs[1].BC_val
        Ci[-1] = network_BCs[1].BC_val
        
        dCidt[-1] = 0.
        
    elif network_BCs[1].BC_type == 'Neumann 1':
        dCidt[-1] = network_BCs[1].BC_val
        
        # loop on jumps
        for j in range(network[-1].Njump):
            dCidt[-1] =  dCidt[-1] \
                       - network[-1].jump[j].forward_freq \
                         *Ci[-1] \
                       + network[-1].jump[j].backward_freq \
                         *Ci[network[-1].jump[j].neighbor_site_idx]
    
    elif network_BCs[1].BC_type == 'Neumann 2':
        dCidt[-1] = network_BCs[1].diff_coeff*(network_BCs[1].BC_val-Ci[-1]) \
                                             /network_BCs[1].diff_length
        
        # loop on jumps
        for j in range(network[-1].Njump):
            dCidt[-1] =  dCidt[-1] \
                       - network[-1].jump[j].forward_freq \
                         *Ci[-1] \
                       + network[-1].jump[j].backward_freq \
                         *Ci[network[-1].jump[j].neighbor_site_idx]
    
    else:
        print('Boundary condition output_BC_type set to incorrect value: ', \
              network_BCs[1].BC_type)
        print('Please use "Dirichlet", "Neumann 1" or , ' \
              '"Neumann 2" as boundary condition type.')
        print('Exiting...')
        sys.exit()
    
    return dCidt