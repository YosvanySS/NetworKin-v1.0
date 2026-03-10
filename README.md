# NetworKin-v1.0

**NetworKin** is a scientific simulation tool for modeling diffusion–reaction
processes on discrete interstitial networks by solving systems of ordinary
differential equations (ODEs).

It supports **linear** and **nonlinear** kinetic models, several **boundary
condition types**, and optional **thermodynamic equilibrium** calculations.

This code is designed for materials science applications involving atomic
transport and kinetic network modeling.

---

## Authors

**Julien Denis**  
Aix-Marseille Université

**Yosvany Silva-Solís**  
Aix-Marseille Université

---

## Features

- Discrete interstitial network modeling
- Linear and nonlinear kinetic ODE systems
- Jump-frequency calculations
- Dirichlet and Neumann boundary conditions
- Automatic normalization of physical quantities
- Thermodynamic equilibrium solver
- HDF5 output format for large datasets
- Modular architecture for extensibility

---

## Repository Structure


---

## Requirements

- Python ≥ 3.9
- numpy
- scipy
- h5py

Install dependencies with:

```bash
pip install numpy scipy h5py
