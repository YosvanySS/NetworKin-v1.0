#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 15:28:50 2025

@author: Julien Denis
@mail:   julien.denis@univ-amu.fr

@author: Yosvany Silva-Solís
@mail:   yosvany.silva-solis@univ-amu.fr

Project: NetworKin
file:    networkin_*simu*_inputs.py
object:  input file for the simulation

"""

#### some necessary libraries and/or modules and/or functions
import numpy as np

from networkin_physical_constants import Da



    ################# NETWORK INPUTS ##################

from networkin_network_inputs import * # here specify the name of the network input file



    ############## SIMULATION PARAMETERS ##############

# temperature [K]
Temp = 300.

# mass of the specie [kg]
sp_mass = 1.00794*Da # Hydrogen

# initial conditions for surface density [m-2]
Ci_0 = np.zeros(np.shape(network_label)) # initial conditions for surface density [m-2]

# time at which the concentration must be saved [s] 
time = np.logspace(-12, 8, num=int(1e4))

# Boundary conditions:
#                     input_BC_type : type of boundary conditions for network input (first site in network_label)
#                     output_BC_type: type of boundary conditions for network output (last site in network_label)
#
#                                     Accepted boundary conditions:
#                                                                   - Dirichlet: imposed surface concentration => BC_val required
#                                                                   - Neumann 1: imposed flux density          => BC_val required
#                                                                   - Neumann 2: imposed flux density through steady state diffusion => BC_val, BC_lambda, BC_L and BC_diff_coeff required
#
#                     input_BC_val  : value of boundary conditions for network input  (first site in network_label)
#                     output_BC_val : value of boundary conditions for network output ( last site in network_label)
#                                     - Dirichlet: imposed surface density at the network location [m-2]
#                                     - Neumann 1: imposed particle flux density at the network location [m-2.s-1]
#                                     - Neumann 2: imposed density at the material boundary [m-3]

# input boundary condition
input_BC_type = 'Dirichlet'
input_BC_val  = 2.976635079807204723e+09 * 6.37e-10 #1.10e-10

    # required if input_BC_type = 'Neumann 2'
input_BC_lambda     = 0. # distance between two interstitial sites [m]
        # quantities related to the material where steady-state diffusion is assumed
input_BC_L          = 0. # material length [m]
input_BC_diff_coeff = 0. # diffusion coefficient of the specie in the material [m2.s-1]

# output boundary condition
output_BC_type = 'Neumann 1'
output_BC_val  = 0.

    # required if input_BC_type = 'Neumann 2'
output_BC_lambda     = 0. # distance between two interstitial sites [m]
        # quantities related to the material where steady-state diffusion is assumed
output_BC_L          = 0. # material length [m]
output_BC_diff_coeff = 0. # diffusion coefficient of the specie in the material [m2.s-1]

# Type of ODE system to solve:
    # 'linear'   : linear ODE system
    # 'nonlinear': nonlinear ODE system
ODE_system_type = 'linear'



    ################## SOLVER INPUTS ##################

# Networkin uses solve_ivp solver from scipy library
# https://docs.scipy.org/doc/scipy-1.13.1/reference/generated/scipy.integrate.solve_ivp.html

method_ivp = 'LSODA'   # Implicit solver used
#method_ivp = 'BDF'    # Implicit solver used
RTOL_ivp = 1.e-6       # Relative tolerance for implicit solver
ATOL_ivp = 1.e-14      # Absolute tolerance for implicit solver [m-2]
max_step_ivp = 100.    # Maximum time step [s]
#max_step_ivp = np.Inf  # Maximum time step [s]



    ###### ANALYTICAL MODEL OF THE LINEAR SYSTEM ######

compute_thermo_equ = True # True or False



    ################### OUTPUT FILE ###################

output_file_name = 'networkin_simu'



    ###################################################