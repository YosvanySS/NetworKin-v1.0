#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 17:39:51 2025

@author: Julien Denis
@mail:   julien.denis@univ-amu.fr

@author: Yosvany Silva-Solís
@mail:   yosvany.silva-solis@univ-amu.fr

Project: NetworKin
file:    networkin_nonlinear_ode_system_jacobian.py
object:  function defining the jacobian of the nonlinear ODE system to solve
    
"""
def networkin_nonlinear_ode_system_jacobian(t,Ci,network,network_BCs):
    
    # t          : time [s]
    # Ci         : concentration network site [m-3]
    # network    : network object
    # network_BCs: network boundary conditions
    
    import sys
    import numpy as np
    
    
    dCidt_jac = np.zeros((Ci.shape[0],Ci.shape[0]))
    
    # Network input boundary
    if network_BCs[0].BC_type == 'Dirichlet':
        dCidt_jac[0,:] = 0.
    
    elif network_BCs[0].BC_type == 'Neumann 1':
        
        # loop on jumps    
        for j in range(network[0].Njump):
            dCidt_jac[0,0] =  dCidt_jac[0,0] \
                            - network[0].jump[j].forward_freq \
                              *(1.-Ci[network[0].jump[j].neighbor_site_idx] \
                                   /network[network[0].jump[j].neighbor_site_idx]\
                                           .Ci_max) \
                            - network[0].jump[j].backward_freq \
                              *Ci[network[0].jump[j].neighbor_site_idx] \
                              /network[0].Ci_max

            dCidt_jac[0,network[0].jump[j].neighbor_site_idx] = \
                + network[0].jump[j].forward_freq \
                  *Ci[0] \
                  /network[network[0].jump[j].neighbor_site_idx]\
                          .Ci_max \
                + network[0].jump[j].backward_freq \
                  *(1.-Ci[0]/network[0].Ci_max)
                  
                  
                  
                  
                  
    elif network_BCs[0].BC_type == 'Neumann 2':
        
        dCidt_jac[0,0] = - network_BCs[0].diff_coeff \
                          /network_BCs[0].diff_length
        
        # loop on jumps    
        for j in range(network[0].Njump):
            dCidt_jac[0,0] =  dCidt_jac[0,0] \
                            - network[0].jump[j].forward_freq \
                              *(1.-Ci[network[0].jump[j].neighbor_site_idx] \
                                   /network[network[0].jump[j].neighbor_site_idx]\
                                           .Ci_max) \
                            - network[0].jump[j].backward_freq \
                              *Ci[network[0].jump[j].neighbor_site_idx] \
                              /network[0].Ci_max

            dCidt_jac[0,network[0].jump[j].neighbor_site_idx] = \
                + network[0].jump[j].forward_freq \
                  *Ci[0] \
                  /network[network[0].jump[j].neighbor_site_idx]\
                          .Ci_max \
                + network[0].jump[j].backward_freq \
                  *(1.-Ci[0]/network[0].Ci_max)
    
    else:
        print('Boundary condition input_BC_type set to incorrect value: ', \
              network_BCs[0].BC_type)
        print('Please use "Dirichlet", "Neumann 1" or , ' \
              '"Neumann 2" as boundary condition type')
        print('Exiting...')
        sys.exit()
    
    
    # Network core
    for i in range(1,Ci.size-1,1):
        
        # loop on jumps
        for j in range(network[i].Njump):
            dCidt_jac[i,i] =  dCidt_jac[i,i] \
                            - network[i].jump[j].forward_freq \
                              *(1.-Ci[network[i].jump[j].neighbor_site_idx] \
                                   /network[network[i].jump[j].neighbor_site_idx]\
                                           .Ci_max) \
                            - network[i].jump[j].backward_freq \
                              *Ci[network[i].jump[j].neighbor_site_idx] \
                              /network[i].Ci_max

            dCidt_jac[i,network[i].jump[j].neighbor_site_idx] = \
                + network[i].jump[j].forward_freq \
                  *Ci[i] \
                  /network[network[i].jump[j].neighbor_site_idx]\
                          .Ci_max \
                + network[i].jump[j].backward_freq \
                  *(1.-Ci[i]/network[i].Ci_max)
    
    
    # Network output boundary
    if network_BCs[1].BC_type == 'Dirichlet':
        dCidt_jac[-1,:] = 0
    
    elif network_BCs[1].BC_type == 'Neumann 1':
        
        # loop on jumps
        for j in range(network[-1].Njump):
            dCidt_jac[-1,-1] =  dCidt_jac[-1,-1] \
                              - network[-1].jump[j].forward_freq \
                                *(1.-Ci[network[-1].jump[j].neighbor_site_idx]\
                               /network[network[-1].jump[j].neighbor_site_idx]\
                                       .Ci_max) \
                              - network[-1].jump[j].backward_freq \
                                *Ci[network[-1].jump[j].neighbor_site_idx] \
                                /network[-1].Ci_max

            dCidt_jac[-1,network[-1].jump[j].neighbor_site_idx] = \
                + network[-1].jump[j].forward_freq \
                  *Ci[-1] \
                  /network[network[-1].jump[j].neighbor_site_idx]\
                          .Ci_max \
                + network[-1].jump[j].backward_freq \
                  *(1.-Ci[-1]/network[-1].Ci_max)
    
    elif network_BCs[1].BC_type == 'Neumann 2':
        dCidt_jac[-1,-1] = - network_BCs[1].diff_coeff \
                            /network_BCs[1].diff_length
        
        # loop on jumps
        for j in range(network[-1].Njump):
            dCidt_jac[-1,-1] =  dCidt_jac[-1,-1] \
                              - network[-1].jump[j].forward_freq \
                                *(1.-Ci[network[-1].jump[j].neighbor_site_idx]\
                               /network[network[-1].jump[j].neighbor_site_idx]\
                                       .Ci_max) \
                              - network[-1].jump[j].backward_freq \
                                *Ci[network[-1].jump[j].neighbor_site_idx] \
                                /network[-1].Ci_max

            dCidt_jac[-1,network[-1].jump[j].neighbor_site_idx] = \
                + network[-1].jump[j].forward_freq \
                  *Ci[-1] \
                  /network[network[-1].jump[j].neighbor_site_idx]\
                          .Ci_max \
                + network[-1].jump[j].backward_freq \
                  *(1.-Ci[-1]/network[-1].Ci_max)
    
    else:
        print('Boundary condition output_BC_type set to incorrect value: ', \
              network_BCs[1].BC_type)
        print('Please use "Dirichlet", "Neumann 1" or , ' \
              '"Neumann 2" as boundary condition type')
        print('Exiting...')
        sys.exit()
    
    
    return dCidt_jac