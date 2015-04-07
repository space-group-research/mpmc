## MPMC

MPMC (Massively Parallel Monte Carlo) is an open-source Monte Carlo package primarily designed for the simulation of liquids, molecular interfaces and functionalized nanoscale materials. It was originally developed by Jon Belof and is now maintained by a group of researchers (Keith McLaughlin, Brant Tudor, Christian Cioce, Adam Hogan and Brian Space) in the Deparment of Chemistry and SMMARTT Materials Research Center at the University of South Florida. MPMC has been applied to the scientific research challenges of nanomaterials for clean energy, environmental sequestration and molecular detection. Developed to run efficiently on the most powerful supercomputing platforms, MPMC can scale to very large numbers of CPUs or GPUs (with support provided for NVidia's CUDA and the OpenCL architecture).

Optimized for the study of nanoscale interfaces, MPMC supports simulation of Coulomb and Lennard-Jones systems, many-body polarization, coupled-dipole van der Waals, quantum rotational statistics, semi-classical quantum effects, advanced importance sampling methods relevant to fluids, and numerous tools for the development of intermolecular potentials.

## Getting Started

#### Libraries

Currently a basic build of MPMC only requires a modern C compiler. Optionally MPMC may be configured to use either [CUDA](https://developer.nvidia.com/cuda-zone) or [OpenCL](https://www.khronos.org/opencl/). Configuring with Coupled-Dipole VDW requires [LAPACK](http://www.netlib.org/lapack/).

#### Downloading MPMC

MPMC can be downloaded with the following command.

```
git clone http://github.com/mpmccode/mpmc
```

#### Compiling MPMC

Once MPMC has been downloaded it may be compiled with the following commands.

```
cd build
./configure generic gcc no no no no
make
```

Make sure to add MPMC to your path after compiling. The different configuration options can be seen by running `./configure`

#### Updating MPMC

MPMC can be updated with the following command.

```
git pull
```

## Publications

* P. Nugent, T. Pham, K. McLaughlin, P.A. Georgiev, W. Lohstroh, J.P. Embs, M.J. Zaworotko, B. Space, and J. Eckert, "Dramatic effect of pore size reduction on the dynamics of hydrogen adsorbed in metal-organic materials.", Journal of Materials Chemistry A, 2, 13884-13891 (2014). 

* K.A. Forrest, T. Pham, K. McLaughlin, A. Hogan, and B. Space, "Insights into an intriguing gas sorption mechanism in a polar metal-organic framework with open-metal sites and narrow channels.", Chemical Communications, 50, 7283-7286 (2014). 

* T. Pham, K.A. Forrest, B. Tudor, S.K. Elsaidi, M.H. Mohamed, K. McLaughlin, C.R. Cioce, M.J. Zaworotko, and B. Space, "Theoretical Investigations of CO2 and CH4 Sorption in an Interpenetrated Diamondoid Metal-Organic Material.", Langmuir, 30(22), 6454-6462 (2014). 

* S.K. Elsaidi, M.H. Mohamed, L. Wojtas, A. Chanthapally, T. Pham, B. Space, J.J. Vittal, and M.J. Zaworotko, "Putting the Squeeze on CH4 and CO2 through Control over Interpenetration in Diamondoid Nets.", Journal of the American Chemical Society, 136(13), 5072-5077 (2014). 

* T. Pham, K.A. Forrest, A. Hogan, K. McLaughlin, J.L. Belof, J. Eckert, B. Space, "Simulations of Hydrogen Sorption in rht-MOF-1: Identifying the Binding Sites Through Explicit Polarization and Quantum Rotation Calculations.", Journal of Materials Chemistry A, 2:2088 (2014). 

* T. Pham, K.A. Forrest, J. Eckert, P.A. Georgiev, A. Mullen, R. Luebke, A.J. Cairns, Y. Belmabkhout, J.F. Eubank, K. McLaughlin, W. Lohstroh, M. Eddaoudi, and B. Space, "Investigating the Gas Sorption Mechanism in an rht-Metal-Organic Framework though Computational Studies.", Journal of Physical Chemistry C, 118(1), 439-456 (2014). 

* Cioce, Christian R., Keith McLaughlin, Jonathan L. Belof, and Brian Space. "A Polarizable and Transferable PHAST N2 Potential for use in Materials Simulation." Journal of Chemical Theory and Computation (2013). 

* Mullen, Ashley L., Tony Pham, Katherine A. Forrest, Christian R. Cioce, Keith McLaughlin, and Brian Space. "A Polarizable and Transferable PHAST CO2 Potential for Materials Simulation." Journal of Chemical Theory and Computation (2013). 

* Nugent, Patrick, Youssef Belmabkhout, Stephen D. Burd, Amy J. Cairns, Ryan Luebke, Katherine Forrest, Tony Pham et al. "Porous materials with optimal adsorption thermodynamics and kinetics for CO2 separation." Nature (2013). 

* McLaughlin, Keith, Christian R. Cioce, Tony Pham, Jonathan L. Belof, and Brian Space. "Efficient calculation of many-body induced electrostatics in molecular systems." The Journal of Chemical Physics 136 (2013): 184112 

* Nugent, Patrick, Vanessah Rhodus, Tony Pham, Brant Tudor, Katherine Forrest, Lukasz Wojtas, Brian Space, and Michael Zaworotko. "Enhancement of CO2 selectivity in a pillared pcu MOM platform through pillar substitution." Chemical Communications 49, no. 16 (2013): 1606-1608. 

* Pham, Tony, Katherine A. Forrest, Keith McLaughlin, Brant Tudor, Patrick Nugent, Adam Hogan, Ashley Mullen, Christian R. Cioce, Michael J. Zaworotko, and Brian Space. "Theoretical Investigations of CO2 and H2 Sorption in an Interpenetrated Square-Pillared Metal-Organic Material." The Journal of Physical Chemistry C (2013). 

* Pham, Tony, Katherine A. Forrest, Patrick Nugent, Youssef Belmabkhout, Ryan Luebke, Mohamed Eddaoudi, Michael J. Zaworotko, and Brian Space. "Understanding Hydrogen Sorption in a Metal–Organic Framework with Open-Metal Sites and Amide Functional Groups." The Journal of Physical Chemistry C 117, no. 18 (2013): 9340-9354. 

* Forrest, Katherine A., Tony Pham, Patrick Nugent, Stephen D. Burd, Ashley Mullen, Lukasz Wojtas, Michael J. Zaworotko, and Brian Space. "Examining the Effects of Different Ring Configurations and Equatorial Fluorine Atom Positions on CO2 Sorption in (bpy) 2SiF6?." Crystal Growth & Design 13, no. 10 (2013): 4542-4548. 

* Forrest, Katherine A., Tony Pham, Adam Hogan, Keith McLaughlin, Brant Tudor, Patrick Nugent, Stephen D. Burd et al. "Computational Studies of CO2 Sorption and Separation in an Ultramicroporous Metal–Organic Material." The Journal of Physical Chemistry C 117, no. 34 (2013): 17687-17698. 

* Nugent, Patrick S., Vanessah Lou Rhodus, Tony Pham, Katherine Forrest, Lukasz Wojtas, Brian Space, and Michael J. Zaworotko. "A Robust Molecular Porous Material with High CO2 Uptake and Selectivity." Journal of the American Chemical Society 135, no. 30 (2013): 10950-10953. 

* Mohamed, Mona H., Sameh K. Elsaidi, Tony Pham, Katherine A. Forrest, Brant Tudor, Lukasz Wojtas, Brian Space, and Michael J. Zaworotko. "Pillar substitution modulates CO 2 affinity in “mmo” topology networks." Chemical Communications 49, no. 84 (2013): 9809-9811. 

* Mohamed, Mona H., Sameh K. Elsaidi, Lukasz Wojtas, Tony Pham, Katherine A. Forrest, Brant Tudor, Brian Space, and Michael J. Zaworotko. "Highly Selective CO2 Uptake in Uninodal 6-Connected “mmo” Nets Based upon MO42–(M= Cr, Mo) Pillars." Journal of the American Chemical Society 134, no. 48 (2012): 19556-19559. 

* Forrest, Katherine A., Tony Pham, Keith McLaughlin, Jonathan L. Belof, Abraham C. Stern, Michael J. Zaworotko, and Brian Space. "Simulation of the Mechanism of Gas Sorption in a Metal–Organic Framework with Open Metal Sites: Molecular Hydrogen in PCN-61." The Journal of Physical Chemistry C 116, no. 29 (2012): 15538-15549. 

* Stern, Abraham C., Jonathan L. Belof, Mohamed Eddaoudi, and Brian Space. "Understanding hydrogen sorption in a polar metal-organic framework with constricted channels." The Journal of Chemical Physics 136 (2012): 034705. 

* McLaughlin, Keith, Christian R. Cioce, Jonathan L. Belof, and Brian Space. "A molecular H2 potential for heterogeneous simulations including polarization and many-body van der Waals interactions." The Journal of Chemical Physics 136 (2012): 194302. 

* Belof, Jonathan L., Abraham C. Stern, and Brian Space. "A Predictive Model of Hydrogen Sorption for Metal− Organic Materials." The Journal of Physical Chemistry C 113, no. 21 (2009): 9316-9320. 

* Belof, Jonathan L., Abraham C. Stern, and Brian Space. "An Accurate and Transferable Intermolecular Diatomic Hydrogen Potential for Condensed Phase Simulation." Journal of Chemical Theory and Computation 4, no. 8 (2008): 1332-1337. 

* Belof, Jonathan L., Abraham C. Stern, Mohamed Eddaoudi, and Brian Space. "On the mechanism of hydrogen storage in a metal-organic framework material." Journal of the American Chemical Society 129, no. 49 (2007): 15202-15210. 
