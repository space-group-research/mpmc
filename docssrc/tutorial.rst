Tutorials
*********

A series of representative example input and .pqr files exist in the tutorials_and_examples folder in the MPMC root directory.

BSS H\ :sub:`2` + MOF-5 Simulation (uVT)
========================================

This is a basic electrostatic/Lennard-Jones-based simulation of H\ :sub:`2` sorption in MOF-5 (using the BSS H\ :sub:`2` model and UFF parameters on MOF-5). The simulation may be started by simply invoking :code:`mpmc MOF5+BSS.inp`. Every corrtime steps averages are output to stdout, including the average potential energy, number of molecules sorbed, Qst, timing information, etc.


Several different output files are created in the simulation directory by default:


* \*.traj.pqr - this is the entire trajectory of the system at every corrtime
* \*.restart.pqr - this is the most recent snapshot of the system
* \*.final.pqr - output upon completion of the simulation, this is the last snapshot of the system
* histogram.dat - this is a histogram of sorbate positions in the .dx format
* frozen.dx - this is the frozen (MOF-5) atoms in the .dx format

BSSP H\ :sub:`2` + MOF-5 Simulation (uVT)
=========================================

This is a very similar simulation to the first tutorial with the main difference arising from the use of the BSSP polarizable H\ :sub:`2` model. The polarizability options in this input file should be robust for most situations, it is advised to double check the :code:`polar_max_iter` produces a converged polarization energy when simulating in a novel system however.

PHAHST H\ :sub:`2` + HKUST-1 Simulation (uVT)
=============================================

This is a similar simulation to the second tutorial with the Lennard-Jones repulsion/dispersion potential replaced by physically ground exponential repulsion/damped dispersion potential:

.. math::

    U_{rd} &= \sum_{i \neq j} \frac{F_0}{\beta_{ij}}e^{\beta_{ij}(r_{ij}-\rho_{ij})}+\sum_{n=3}^5 f_{2n}(\beta r_{ij} ) \frac{C_{2n,ij}}{r_{ij}^{2n}} \\
    f_{2n}( \beta r_{ij} ) &= 1 - e^{-\beta r_{ij}} \sum_{k=0}^{2n} \frac{(\beta r_{ij})^k}{k!} \\
    \rho_{ij} &= \frac{1}{2}(\rho_{ii} + \rho_{jj}) \\
    \beta_{ij} &= 2 \frac{\beta_{ii} \beta_{jj}}{\beta_{ii}+\beta_{jj}}\\
    C_{2n,ij} &= \sqrt{C_{2n,ii} C_{2n,jj}}

The use of the :code:`cavity_autoreject_repulsion` and :code:`polar_precision` options ensure that nuclear fusion doesn't happen due to a polarization catastrophe and finite exponential repulsion potential energy. For more information on the PHAHST force field please see: https://pubs.acs.org/doi/abs/10.1021/acs.jctc.0c00837\ .

Simulated Annealing
===================

This is a sample simulated annealing script for H\ :sub:`2` in the MOF NOTT-112. Simulated annealing is typically used to identify minimums in the potential energy surface, i.e. binding sites in material simulations. The input pqr contains the MOF with a single H\ :sub:`2` molecule about 2.5 angstrom from the CuC open metal site. The temperature starts at 40 K and is decreased by a factor of 0.99999 each MC step. At the end of the simulation the H\ :sub:`2` molecule will be in its lowest energy configuration around the open metal site. Small values for :code:`move_factor` and :code:`rot_factor` prevent the H\ :sub:`2` from straying away too far from its initial position.

Replaying a Trajectory
======================

The "replay" ensemble takes in a \*.traj.pqr file, either produced from a previous run or created artificially, and recalculates averages without running a simulation. This is useful in case some quantity needs to be recalculated without running a full simulation. The provided MOF5+BSS.traj-00000.pqr file was produced from the first tutorial.

Single Point Energy and Quantum Rotation Calculations
=====================================================

This folder contains the PQR files with the BSS model                                  
located at all four H\ :sub:`2` sorption sites in MOF-5 (alpha. beta, gamma, and delta).                               
An input file to run the quantum rotation calculations is also included in the                                
folder. Note that running quantum rotations calculations is predicated upon                                 
turning the quantum rotations part of the code on when compiling MPMC. By                                     
executing quantum rotation calculations on the four H\ :sub:`2` sorption sites, you                                    
should obtain rotational levels that are very close to those shown in Table 1                                 
of Ivana's JCP 2012 paper (see link) for the respective sites. Note that the :code:`total_energy`
ensemble here calculates the potential energy (and rotational eigenspectrum) but does not perform any
Monte Carlo steps.

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

Potential Energy Surface Fitting (surf_multi_fit)
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

BSSP H\ :sub:`2` + MOF5 Simulation (uVT) with NVIDIA CUDA
=========================================================

This is the second tutorial with the sole exception of the switch :code:`cuda on`. MPMC needs to be compiled with cuda and a cuda capable GPU must be present. Performance may be compared with the second tutorial.

1D Chain Replay
===============

Here we use :code:`ensemble replay` to read in a series of increasingly
larger 1D-chain samples (starting from 2-atoms and up to 512). Each component of the energy is re-calculated for each sample, for
various potential options (the different input files). One may take the energy output files and process them (via scale.sh) to 
check for the small-size scalability and accuracy of our calculations.

3D Chain Replay
===============

Here we use 'ensemble replay' to read in a series of increasingly
larger crystal samples (starting from 2-atoms and up to 1024).

Each component of the energy is re-calculated for each sample, for
various potential options (the different input files). One may take the energy output files and process them (via scale.sh) to 
check for the small-size scalability and accuracy of our calculations.

