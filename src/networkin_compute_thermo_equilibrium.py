#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 09:41:34 2025

@author: Julien Denis
@mail:   julien.denis@univ-amu.fr

@author: Yosvany Silva-Solís
@mail:   yosvany.silva-solis@univ-amu.fr

Project: NetworKin
file:    networkin_compute_thermo_equilibrium.py
object:  function to compute the thermodynamic equilibrium of the system

"""

def networkin_compute_thermo_equilibrium(network,network_BCs):
    
    # network    : network object
    # network_BCs: network boundary conditions    
    
    import numpy as np
    
    
    Ci_thermo_eq = np.zeros(len(network))
    
    
    if     (network_BCs[0].BC_type == 'Dirichlet' and
            network_BCs[1].BC_type == 'Neumann 1' and
            network_BCs[1].BC_val == 0.) \
        or (network_BCs[1].BC_type == 'Dirichlet' and
            network_BCs[0].BC_type == 'Neumann 1' and
            network_BCs[0].BC_val == 0.):
        
        
        if network_BCs[0].BC_type == 'Dirichlet':
            Ci_thermo_eq[0] = network_BCs[0].BC_val
        elif network_BCs[1].BC_type == 'Dirichlet':
            Ci_thermo_eq[1] = network_BCs[1].BC_val
        
        idx_zeros = np.where(Ci_thermo_eq==0)
        
        while idx_zeros[0].size != 0:
            
            for i in idx_zeros[0]:
                
                # loop on jumps
                for j in range(network[i].Njump):
                    
                    if (Ci_thermo_eq[network[i].jump[j].neighbor_site_idx]
                        != 0.):
                        
                        Ci_thermo_eq[i] = network[i].jump[j].backward_freq \
                                         /network[i].jump[j].forward_freq \
                                         * Ci_thermo_eq[network[i].jump[j].\
                                                        neighbor_site_idx]
                        
            idx_zeros = np.where(Ci_thermo_eq==0)
            
    
    else:
        print('Thermodynamic equilibrium conditions are not fulfilled ' \
              'for the boundary conditions.')

    return Ci_thermo_eq