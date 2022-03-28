# Lithium Sulfur Battery Thermodynamic Equilibria Model Manuscript
The thermodynamic equilibrium voltage of a lithium sulfur cell can be found using a Galvanostatic Intermittent Titration Technique (GITT) relaxation experiment. We model this equilibrium voltage using both a 1D physics-based model (thermo + transport + kinetics) and a thermodynamic model. The thermodynamic model is based on solving 8-9 algebraic equations in different regions of depth-of-discharge (DOD). This repository contains all the experimental and simulated data, the python code used to make figures in the corresponding manuscript, and the Julia code for the thermodynamic equilibrium model.

For more information, or if this code is used, please go to or cite the following paper:
- Paper citation will be available upon acceptance


### Folder structure
*data*:
- *thermo_initial_guesses.csv*  for thermo equilibria model
- all other data used to make figures


*jupyter_notebooks*:

- *jlThermoModel.ipynb* notebook that walks through how to use LiS_thermo.jl to solve thermodynamic equilibria

- notebook containing code that generates figures for the manuscript

*src*:
- *LiS_thermo.jl* contains Julia code to solve thermodynamic model.
- *gitt.py* python file used to analyze data used for figures


### Software dependencies
 - PYTHON v3.6.5 (all other packages for making figures can be found in and installed using environment.yml)
 - Git Bash for Windows


 - JULIA v1.1.1
 - CSV v0.8.5
 - DataFrames v1.3.1
 - DifferentialEquations v6.4.0
 - FileIO v1.12.0
 - LineSearches v7.1.1
 - NLsolve v4.5.1
 - PyPlot v2.10.0
 - Sundials v3.6.0
 - DelimitedFiles
 - LinearAlgebra


### Getting started with Julia
- Instructions for downloading Julia can be found here: https://julialang.org/downloads/platform/
- Julia getting started documentation: https://docs.julialang.org/en/v1/manual/getting-started/
- How to install a Julia package and use it in a jupyter notebook: https://datatofish.com/install-package-julia/
