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
Aix Marseille Univ, CNRS, PIIM, Marseille, France

**Yosvany Silva-Solís**  
Aix Marseille Univ, PIIM, Marseille, France

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
```

---

## Simulation Workflow

### Step 1 — Input parsing
Loads simulation parameters from the user input file.

### Step 2 — Variable initialization
- Validates boundary conditions
- Computes jump frequencies
- Normalizes time and concentrations
- Prepares ODE solver parameters

### Step 3 — Network construction
Builds the interstitial network and allowed atomic jumps.

### Step 4 — Boundary conditions
Applies physical boundary conditions to the network.

### Step 5 — Time integration
Solves the kinetic ODE system using `scipy.integrate.solve_ivp`.

### Step 6 — Thermodynamic equilibrium (optional)
Computes equilibrium concentrations if requested.

### Step 7 — Output storage
Saves results in HDF5 format.

---

## Datasets

| Dataset           | Description                     |
| ----------------- | ------------------------------- |
| `time`            | Simulation time grid            |
| `Ci`              | Site concentrations vs time     |
| `network_label`   | Names of network sites          |
| `is_Ci_thermo_eq` | Equilibrium computed flag       |
| `Ci_thermo_eq`    | Thermodynamic equilibrium state |

---

## Citation

If you use **NetworKin** in your scientific work, please cite the article where the method and software are presented:

**Hydrogen transport at metallic interfaces: modeling the tungsten-copper system.**
*Silva-Solís, Y., Denis, J., Hodille, E. A., Ferro, Y.*
Physical Review Materials **Volume_xxx**, Pages_xxx (Year_xxx)  
DOI: https://doi.org/[DOI_xxx]

### BibTeX

```bibtex
@article{networkin_paper,
  author  = {Silva-Solís, Y., Denis, J., Hodille, E. A., Ferro, Y.},
  title   = {Hydrogen transport at metallic interfaces: modeling the tungsten-copper system.},
  journal = {Physical Review Materials},
  year    = {2026},
  volume  = {xxxxx},
  pages   = {xxxxx},
  doi     = {xxxxx}
}
```
