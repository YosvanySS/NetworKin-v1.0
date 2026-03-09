#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 15:28:50 2025

@author: Julien Denis
@mail:   julien.denis@univ-amu.fr

@author: Yosvany Silva-Solís
@mail:   yosvany.silva-solis@univ-amu.fr

Project: NetworKin
file:    networkin_*network*_inputs.py
object:  input file for the network

"""

#### some necessary libraries and/or modules and/or functions
import numpy as np


    ################## NETWORK INPUTS #################

    # network_label: tuple with the label of the different sites in the network
    #                the order sets the order of the ODE system
    #                Restrictions : first and last sites are the network boundaries (in/out)

network_label = ('d', '31', '86', 'c')


    # network_Ci_max: maximum site density [m-2]
        # required if nonlinear simulation
network_Ci_max = np.zeros(np.shape(network_label))


    # Jumps in the network: jump = tuple (site 1 label, site 2 label,
    #                                     jump 1->2:
    #                                           - logical to use Wert and Zener for frequency pre-exp factor,
    #                                           - False: frequency pre-exp factor [s-1]; True: distance between 1 and 2 [Angström],
    #                                           - activation energy [eV],
    #                                     jump 2->1:
    #                                           - logical to use Wert and Zener for frequency pre-exp factor,
    #                                           - False: frequency pre-exp factor [s-1]; True: distance between 2 and 1 [Angström],
    #                                           - activation energy [eV]

network_jumps = (( 'd', '31', False, 1.e10, 0.3, False, 1.e10, 0.4),
                 ('31', '86', False, 1.e10, 0.3, False, 1.e10, 0.8),
                 ('86',  'c', False, 1.e10, 0.3, False, 1.e10, 0.4))


    ###################################################