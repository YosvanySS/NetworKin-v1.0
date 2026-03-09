#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 23:20:25 2025

@author: Julien Denis
@mail:   julien.denis@univ-amu.fr

@author: Yosvany Silva-Solís
@mail:   yosvany.silva-solis@univ-amu.fr

Project: NetworKin
file:    networkin_classes.py
object:  definition of classes used in NetworKin
    
"""

# Interstitial site
class interstitial_site:
    def __init__(self, site_label, Ci_max=None):
        
        # site_label: label of the different sites
        
        self.site_label = site_label
        
        if Ci_max:
            self.Ci_max = Ci_max
        
        self.Njump = 0
        self.jump  = []


# Jump from interstitial site to neighbouring site
class jump:
    def __init__(self, neighbor_site_label, neighbor_site_idx, \
                 forward_freq, backward_freq):
        
        # neighbor_site_label: neighbor site label
        # neighbor_site_idx  : neighbor site index
        # forward_freq       : forward frequency  [s]
        # backward_freq      : backward frequency [s]
        
        self.neighbor_site_label = neighbor_site_label
        self.neighbor_site_idx   = neighbor_site_idx
        self.forward_freq        = forward_freq
        self.backward_freq       = backward_freq


# Boundary conditions
class BoundCond:
    def __init__(self, BC_type, BC_val, diff_coeff=None, diff_length=None):
        
        # BC_type: type of BC, "Dirichlet", "Neumann 1" or "Neumann 2"
        # BC_val : value of BC
        #
        # only if BC_type == "Neumann 2"
        #   diff_coeff : diffusion coefficient
        #   diff_length: diffusion length 
        
        self.BC_type = BC_type
        self.BC_val  = BC_val
        
        if BC_type == 'Neumann 2':
            self.diff_coeff  = diff_coeff
            self.diff_length = diff_length