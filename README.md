# Lithium Sulfur Battery Thermodynamic Equilibria Model Manuscript
The thermodynamic equilibrium voltage of a lithium sulfur cell can be found using a Galvanostatic Intermittent Titration Technique (GITT) relaxation experiment. We model this equilibrium voltage using both a 1D full physics model and a thermodynamic model. The thermodynamic model is based off solving 8-9 algebraic equations in different regions of depth-of-discharge (DOD). This repository contains all the experimental and simulated data, the python code used to make figures in the corresponding manuscript, and the Julia code for the thermodynamic equilibrium model.

For more information, or if this code is used, please go to or cite the following paper:
- Paper citation will be available upon acceptance


### Folder structure
*data*:
- *thermo_initial_guesses.csv*


*jupyter_notebooks*:

- *jlThermoModel.ipynb* notebook that walks through how to use LiS_thermo.jl to solve thermodynamic equilibria

- notebook containing code that generates figures for the manuscript

*src*:
- *LiS_thermo.jl* contains Julia code to solve thermodynamic model.


### Software dependencies



### Getting started with Julia
- Instructions for downloading Julia can be found here: https://julialang.org/downloads/platform/
- Julia getting started documentation: https://docs.julialang.org/en/v1/manual/getting-started/
- How to install a Julia package and use it in a jupyter notebook: https://datatofish.com/install-package-julia/
