# Bacterial Growth with Chemical Nutrient & Hindering Byproduct

**Goal**: This report applies an optimal control approach to identify the ideal amount of chemical nutrient to add to a bacteria population to maximize population growth and minimize chemical usage. The chemical produces a byproduct that hinders bacterial growth, and this byproduct has a stronger effect on smaller populations.

**Methods**: 
- _Mathematical:_ I cast an optimal control problem by defining an ODE for bacterial growth, setting appropriate constraints on chemical concentration via the control, assuming initial conditions and parameter constraints, and finally constructing the Hamiltonian. Using the Hamiltonian, I calculated the adjoint equation, transversality condition, and the optimality condition.
- _Computational:_ I used a forward-backward sweep algorithm in MATLAB to solve the adjoint equation and differential equation (for bacterial growth) subject to the constraints of the transversality and optimality conditions with initial conditions and parameter ranges. 

## Repo Guide

In the "Full Report" file, you'll find my entire typeset report, including images and equations.

All .png files with the "m1" title are figures generated for this project.

In the "MATLAB" section, you'll find the MATLAB code I used to generate figures for the report. I used the Runge-Kutta 4 forward-backward sweep method to numerically estimate each graph because the differential equations we produced were impossible to solve analytically.
