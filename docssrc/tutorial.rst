Tutorials
*********

WIP

Bulk BSS H\ :sub:`2` Simulation (NVT)
=====================================

Basic LJ/ES bulk sorption modeling for H2 (using BSS model).

(repulsion-dispersion + electrostatics only; no polarizaion)

BSSP H\ :sub:`2` + MOF-5 Simulation (uVT)
=========================================

Basic Polarized bulk sorption modeling for H2 (using BSSP model).

(U = repulsion-dispersion + electrostatics + polarization)

PHAHST H\ :sub:`2` + HKUST-1 Simulation (uVT)
=============================================

Simulated Annealing
===================

This is a sample simulated annealing (using NVT) script for the MOF NOTT-112.
Two files are included.

input.pqr (1) contains the MOF with a single H2 molecule about 2.5A from
the CuC open metal site.

sa.inp (2) contains instructions for MPMC to run the annealing, starting from
40K and decreasing the temperature by a factor of 0.99999 each MC step. Pressure
is kept constant at 0.1atm.

S.A. generally takes less time than uptake simulations for several reasons:
Being an NVT model, the number of sorbate molecules is fixed (e.g. at 1)
The simluation need not 'equilibrate'; only run until a desired final 
temperature is reached.

Replaying a Trajectory
======================

Single Point Energy and Quantum Rotation Calculations
=====================================================

Included is a TGZ file containing the PDB files with the BSS model                                  
located at all four H2 sorption sites in MOF-5 (alpha. beta, gamma, and delta).                               
An input file to run the quantum rotation calculations is also included in the                                
TGZ file. Note that running quantum rotations calculations is predicated upon                                 
turning the quantum rotations part of the code on when compiling MPMC. By                                     
executing quantum rotation calculations on the four H2 sorption sites, you                                    
should obtain rotational levels that are very close to those shown in Table 1                                 
of Ivana's JCP 2012 paper (see included paper) for the respective sites. 

Paper: http://scitation.aip.org/content/aip/journal/jcp/137/1/10.1063/1.4730906

Potential Energy Surface Fitting (surf_fit)
===========================================

These are sample inputs for fitting noble gases in MPMC.

Use, for example,

/path/to/mpmc helium_fitting.input

to run.

All 3 files are required to run PES
\*.input (MPMC input)
\*.pqr (initial parameters/coordinates)
\*.dat (PES input data)

Potential Energy Surface Fitting (surf_fit and surf_multi_fit)
==============================================================

Multisorbate uVT
================

These input files will run a multi-sorbate simulation in an rht-MOF (NOTT-112).
The two sorbates are H2 (BSSP) and CO2 (PHAST). I have done 3 sorbates in this
system with success (e.g. N2). Pretty sure code is generalized to any number of
sorbates.

Files
1. insert.pqr
Only needed for multi-sorbate simulations (for gas sorption)
Includes the sorbate models to be included
2. input.pqr
Includes the MOF coordinates and initial sorbate (only 1 is needed)
3. \*.inp
The actual MPMC input file, which cites the above 2 sub-inputs.

This is one of the largest (in phase-space and atom count) systems we've ever
worked on, so it takes a while. Especially for large corrtime in \*.inp.

Coupled-Dipole van der Waals
============================

BSSP H\ :sub:`2` + MOF5 Simulation (uVT) with NVIDIA CUDA
=========================================================

1D Chain Replay
===============

Here we use 'ensemble replay' to read in a series of increasingly
larger 1D-chain samples (starting from 2-atoms and up to 512).

Each component of the energy is re-calcualte for each sample, for
various potential options (one sec of potential options per input
file). One may take the energy output files and process them (via scale.sh) to 
check for the small-size scalability and accuracy of our calculations.

3D Chain Replay
===============

Here we use 'ensemble replay' to read in a series of increasingly
larger crystal samples (starting from 2-atoms and up to 1024).

Each component of the energy is re-calcualte for each sample, for
various potential options (one sec of potential options per input
file). One may take the energy output files and process them (via scale.sh) to 
check for the small-size scalability and accuracy of our calculations.
