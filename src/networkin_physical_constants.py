#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 15:04:40 2025

@author: Julien Denis
@mail:   julien.denis@univ-amu.fr

@author: Yosvany Silva-Solís
@mail:   yosvany.silva-solis@univ-amu.fr

Project:   NetworKin
file:      networkin_physical_constants.py
object:    module defining the physical constants used in NetworKin
reference: Bureau International des Poids et Mesures. (2025). Le Système international d'unités/The International System of Units [Brochure]. 9th edition. https://doi.org/10.59161/AUEZ1291

"""

# Planck constant [J.s]
h_P = 6.62607015e-34

# elementary charge [C]
e_charge = 1.602176634e-19

# Boltzmann constant
k_B_SI = 1.380649e-23 # [J.K-1]
k_B = k_B_SI/e_charge # [eV.K-1]

# Avogadro constant [mol-1]
N_A = 6.02214076e-23

# Dalton / unified atomic mass unit [kg]
Da = 1.66053906892e-27