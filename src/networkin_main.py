#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 17:08:50 2025

@author: Julien Denis
@mail:   julien.denis@univ-amu.fr

@author: Yosvany Silva-Solís
@mail:   yosvany.silva-solis@univ-amu.fr

Project: NetworKin
file:    networkin_main.py
object:  main file of NetworKin

"""

####  reset
#from IPython import get_ipython;
#get_ipython().run_line_magic('reset', '-sf')



#### Libraries, classes and functions
import sys
import numpy as np
import scipy.integrate as sp_int
import h5py

from networkin_logo                       import networkin_logo
from networkin_classes                    import interstitial_site, \
                                                 jump, \
                                                 BoundCond
from networkin_frequency                  import networkin_frequency
from networkin_print_network              import networkin_print_network
from networkin_linear_ode_system          import networkin_linear_ode_system
from networkin_linear_ode_system_jacobian import \
    networkin_linear_ode_system_jacobian
from networkin_nonlinear_ode_system       import networkin_nonlinear_ode_system
from networkin_nonlinear_ode_system_jacobian import \
    networkin_nonlinear_ode_system_jacobian
from networkin_compute_thermo_equilibrium import \
    networkin_compute_thermo_equilibrium



#### We print the NetworKin logo
networkin_logo()



#### Step 1: definition of inputs

print('---------------------   Reading inputs...   ---------------------')
print('    ')

from networkin_simu_inputs import * # here provide the name of the simulation input file

print('The simulation inputs are read.')
print('    ')
print('-----------------------------------------------------------------')
print('    ')
print('    ')



#### Step 2: preparing some variables

print('-----------------   Initialising variables...   -----------------')
print('    ')

    # We check if the selected boundary conditions are correct
        # Input BC 
if input_BC_type not in ['Dirichlet', 'Neumann 1', 'Neumann 2']:
    print('Boundary condition input_BC_type set to incorrect value: ', \
          input_BC_type)
    print('Please use "Dirichlet", "Neumann 1" or , ' \
          '"Neumann 2" as boundary condition type.')
    print('Exiting...')
    sys.exit()

        # Output BC 
if output_BC_type not in ['Dirichlet', 'Neumann 1', 'Neumann 2']:
    print('Boundary condition output_BC_type set to incorrect value: ', \
          output_BC_type)
    print('Please use "Dirichlet", "Neumann 1" or , ' \
          '"Neumann 2" as boundary condition type.')
    print('Exiting...')
    sys.exit()


    # We calculate the frequencies for each jump
jump_forward_freq  = np.zeros(np.shape(network_jumps)[0])
jump_backward_freq = np.zeros(np.shape(network_jumps)[0])

jump_forward_nu0  = np.zeros(np.shape(network_jumps)[0])
jump_backward_nu0 = np.zeros(np.shape(network_jumps)[0])

for i in range(len(network_jumps)):
    jump_forward_nu0[i] , jump_forward_freq[i]  \
        = networkin_frequency(network_jumps[i][2], sp_mass, \
                              network_jumps[i][3], network_jumps[i][4], Temp)
    jump_backward_nu0[i], jump_backward_freq[i] \
        = networkin_frequency(network_jumps[i][5], sp_mass, \
                              network_jumps[i][6], network_jumps[i][7] ,Temp)


    # We set the normalisation parameters of the ODE system
        # time
maxfreq = np.max(np.maximum(jump_forward_freq,jump_backward_freq))
minfreq = np.min(np.minimum(jump_forward_freq,jump_backward_freq))

t0 = 1./10.**((np.log10(maxfreq)+np.log10(minfreq))/2.)

        # surface density
if input_BC_type == 'Dirichlet' and output_BC_type != 'Dirichlet':
    C0 = input_BC_val
elif output_BC_type == 'Dirichlet' and input_BC_type != 'Dirichlet':
    C0 = output_BC_val
elif input_BC_type == 'Dirichlet' and output_BC_type == 'Dirichlet':
    C0 = np.max([input_BC_val,output_BC_val])
elif input_BC_type == 'Neumann 1' and output_BC_type == 'Neumann 2':
    C0 = (input_BC_val*output_BC_L/output_BC_diff_coeff+output_BC_val)\
        *output_BC_lambda
elif input_BC_type == 'Neumann 2' and output_BC_type == 'Neumann 1':
    C0 = (output_BC_val*input_BC_L/input_BC_diff_coeff+input_BC_val)\
        *input_BC_lambda
else:
    C0 = 1.
    
        # Flux density
flux_dens_0 = C0/t0

    # Normalisation of the different terms
        # Frequencies
jump_forward_freq  = t0*jump_forward_freq
jump_backward_freq = t0*jump_backward_freq

        # Time
time = time/t0

        # Concentrations
            # initial conditions
Ci_0 = Ci_0/C0

            # BCs
if input_BC_type == 'Dirichlet':
    input_BC_val        = input_BC_val/C0
    
elif input_BC_type == 'Neumann 1':
    input_BC_val        = input_BC_val/flux_dens_0
    
elif input_BC_type == 'Neumann 2':
    # normalisation parameters
    L0 = input_BC_L         # length
    n0 = C0/input_BC_lambda # density
    D0 = C0/t0*L0/n0        # diffusion coefficient
    
    # normalisation of the terms
    input_BC_val        = input_BC_val/n0
    input_BC_L          = input_BC_L/L0
    input_BC_diff_coeff = input_BC_diff_coeff/D0


if output_BC_type == 'Dirichlet':
    output_BC_val        = output_BC_val/C0
    
elif output_BC_type == 'Neumann 1':
    output_BC_val        = output_BC_val/flux_dens_0
    
elif output_BC_type == 'Neumann 2':
    # normalisation parameters
    L0 = output_BC_L         # length
    n0 = C0/output_BC_lambda # density
    D0 = C0/t0*L0/n0         # diffusion coefficient
    
    # normalisation of the terms
    output_BC_val        = output_BC_val/n0
    output_BC_L          = output_BC_L/L0
    output_BC_diff_coeff = output_BC_diff_coeff/D0


            # maximum concentrations
if ODE_system_type == 'nonlinear':
    
    # test if network_Ci_max variable exists
    try:
        network_Ci_max
    except:
        print('"network_Ci_max" variable is not defined.')
        print('It is required for nonlinear NetworKin system.')
        print('Exiting...')
        sys.exit()
    else:
        network_Ci_max = network_Ci_max/C0

            # solver variable
ATOL_ivp = ATOL_ivp/C0

    # Correction of initial conditions if Dirichlet BCs are imposed
if input_BC_type == 'Dirichlet':
    Ci_0[0] = input_BC_val

if output_BC_type == 'Dirichlet':
    Ci_0[-1] = output_BC_val

print('Variables initialised.')
print('    ')
print('-----------------------------------------------------------------')
print('    ')
print('    ')



#### Step 3: setting up the network

print('------------------   Network configuration...  ------------------')
print('    ')

dummy_network = []


    ## network initialisation
if ODE_system_type == 'linear':
    for i in range(len(network_label)):
        dummy_network.append(interstitial_site(network_label[i]))
    
elif ODE_system_type == 'nonlinear':
    
    # test if network_Ci_max variable exists
    try:
        network_Ci_max
    except:
        print('"network_Ci_max" variable is not defined.')
        print('It is required for nonlinear NetworKin system.')
        print('Exiting...')
        sys.exit()    
    
    # we init the interstitial site structure
    for i in range(len(network_label)):
        dummy_network.append(interstitial_site(network_label[i], \
                                               network_Ci_max[i]))
        
else:
    print('ODE system type ODE_system_type set to incorrect value: ', \
          ODE_system_type)
    print('Please use "linear" or "nonlinear" as ODE system type.')
    print('Exiting...')
    sys.exit()


    ## filling jumps

for i in range(len(network_jumps)):
    
        # for site 1
    idx_site          = network_label.index(network_jumps[i][0])
    idx_site_neighbor = network_label.index(network_jumps[i][1])
    
    dummy_network[idx_site].Njump = dummy_network[idx_site].Njump + 1
    
    forward_freq  = jump_forward_freq[i]
    backward_freq = jump_backward_freq[i]
    
    dummy_network[idx_site].jump.append(jump(network_jumps[i][1], \
                                             idx_site_neighbor, forward_freq, \
                                             backward_freq))
    
        # for site 2
    idx_site          = network_label.index(network_jumps[i][1])
    idx_site_neighbor = network_label.index(network_jumps[i][0])
    
    dummy_network[idx_site].Njump = dummy_network[idx_site].Njump + 1
    
    forward_freq  = jump_backward_freq[i]
    backward_freq = jump_forward_freq[i]
    
    dummy_network[idx_site].jump.append(jump(network_jumps[i][0], \
                                             idx_site_neighbor, forward_freq, \
                                             backward_freq))


network = tuple(dummy_network)
del dummy_network

print('Network configured:')
print('    ')

networkin_print_network(network,t0)

print('-----------------------------------------------------------------')
print('    ')
print('    ')



#### Step 4: setting up the network boundary conditions

dummy_network_BCs = []

    # input BC
if input_BC_type in ['Dirichlet','Neumann 1']:
    dummy_network_BCs.append(BoundCond(input_BC_type, input_BC_val))
elif input_BC_type == 'Neumann 2':
    dummy_network_BCs.append(BoundCond(input_BC_type, input_BC_val, \
                                       input_BC_diff_coeff, input_BC_L))

    # output BC
if output_BC_type in ['Dirichlet','Neumann 1']:
    dummy_network_BCs.append(BoundCond(output_BC_type, output_BC_val))
elif output_BC_type == 'Neumann 2':
    dummy_network_BCs.append(BoundCond(output_BC_type, output_BC_val, \
                                       output_BC_diff_coeff, output_BC_L))

network_BCs = tuple(dummy_network_BCs)
del dummy_network_BCs

#sys.exit()

#### Step 5: time integration

print('--------------   NetworKin integration pending...  --------------')
print('    ')

if ODE_system_type == 'linear':
    sol = sp_int.solve_ivp(networkin_linear_ode_system, \
                           [0, time[-1]], Ci_0, \
                           method=method_ivp, t_eval=time, \
                           args=(network,network_BCs), \
                           max_step=max_step_ivp/t0, \
                           rtol = RTOL_ivp, atol = ATOL_ivp, \
                           jac=networkin_linear_ode_system_jacobian)
    
elif ODE_system_type == 'nonlinear':
    sol = sp_int.solve_ivp(networkin_nonlinear_ode_system, \
                           [0, time[-1]], Ci_0, \
                           method=method_ivp, t_eval=time, \
                           args=(network,network_BCs), \
                           max_step=max_step_ivp/t0, \
                           rtol = RTOL_ivp, atol = ATOL_ivp, \
                           jac=networkin_nonlinear_ode_system_jacobian) 
else:
    print('ODE system type ODE_system_type set to incorrect value: ', \
          ODE_system_type)
    print('Please use "linear" or "nonlinear" as ODE system type.')
    print('Exiting...')
    sys.exit()


if sol.status == 0:
    print(sol.message)
elif sol.status == -1:
    print('Integration step failed.')
    print(sol.message)
    print('    ')
    print('Exiting...')
    sys.exit()

print('    ')
print('-----------------------------------------------------------------')
print('    ')
print('    ')



#### Step 6: we compute the thermodynamic equilibrium of the system is required

if compute_thermo_equ:
    print('-----   Computing the network thermodynamic equilibrium...  -----')
    print('    ')
    
    Ci_thermo_eq = networkin_compute_thermo_equilibrium(network,network_BCs)
    
    print('The thermodynamic equilibrium of the system is calculated.')
    print('    ')
    print('-----------------------------------------------------------------')
    print('    ')
    print('    ')



#### Step 7: saving results

print('---------------------   Saving results...   ---------------------')
print('    ')

f = h5py.File(output_file_name+'.h5','w')

network_label_dset = f.create_dataset('network_label', \
                                      shape = len(network_label), \
                                      data=network_label, \
                                      dtype=h5py.string_dtype())

    # solv_ivp
time_dset = f.create_dataset('time', shape=sol.t.shape, data=t0*sol.t, \
                             dtype='float64')
Ci_dset   = f.create_dataset('Ci'  , shape=sol.y.shape, data=C0*sol.y, \
                             dtype='float64')

    # thermodynamic equilibrium
is_Ci_thermo_eq_dset = f.create_dataset('is_Ci_thermo_eq', \
                                        data=int(compute_thermo_equ))

if compute_thermo_equ:    
    Ci_thermo_eq_dset = f.create_dataset('Ci_thermo_eq', \
                                         shape=Ci_thermo_eq.shape, \
                                         data=C0*Ci_thermo_eq, dtype='float64')

f.close()

print('The simulation results are saved.')
print('    ')
print('-----------------------------------------------------------------')
print('    ')
