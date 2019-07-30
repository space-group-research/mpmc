## MPMC

MPMC (Massively Parallel Monte Carlo) is an open-source Monte Carlo package primarily designed for the simulation of liquids, molecular interfaces and functionalized nanoscale materials. It was originally developed by Jon Belof and is now maintained by a group of researchers (Keith McLaughlin, Brant Tudor, Christian Cioce, Adam Hogan, Douglas Franz and Brian Space) in the Department of Chemistry and SMMARTT Materials Research Center at the University of South Florida. MPMC has been applied to the scientific research challenges of nanomaterials for clean energy, environmental sequestration and molecular detection. Developed to run efficiently on the most powerful supercomputing platforms, MPMC can scale to very large numbers of CPUs or GPUs (with support provided for NVidia's CUDA and the OpenCL architecture).

Optimized for the study of nanoscale interfaces, MPMC supports simulation of Coulomb and Lennard-Jones systems, many-body polarization, coupled-dipole van der Waals, quantum rotational statistics, semi-classical quantum effects, advanced importance sampling methods relevant to fluids, and numerous tools for the development of intermolecular potentials.

## Getting Started

#### Libraries

Currently a basic build of MPMC only requires a modern C compiler and CMake. Optionally MPMC may be configured to use either [CUDA](https://developer.nvidia.com/cuda-zone) or [OpenCL](https://www.khronos.org/opencl/). Configuring with Coupled-Dipole VDW or QM Rotation requires [LAPACK](http://www.netlib.org/lapack/).

#### Supported Platforms
MPMC supports compilation on Linux, macOS, and Windows; however MPMC is primarily tested on Linux and support is not guaranteed on other platforms.

#### Downloading MPMC

MPMC can be downloaded with the following command.

```
git clone https://github.com/mpmccode/mpmc
```

#### Compiling MPMC

Once MPMC has been downloaded it may be compiled with the following commands.

```
bash compile.sh
```

or

```
mkdir build
cd build
cmake -DQM_ROTATION=OFF -DVDW=OFF -DMPI=OFF -DOPENCL=OFF -DCUDA=OFF -DCMAKE_BUILD_TYPE=Release -Wno-dev ../
make
```

Make sure to add MPMC to your path after compiling.

#### Updating MPMC

MPMC can be updated with the following command.

```
git pull
```
#### MPMC Test Suite

An end-to-end test suite for MPMC is currently under development.
If cloning MPMC anew, use the following to include the tests:
```
git clone https://github.com/mpmccode/mpmc --recurse-submodules
```
To clone the submodule into an existing MPMC installation, use this instead:
```
cd ${local_mpmc_dir}
git submodule init
git submodule update
```
To run the tests, make sure you have Python 3 installed, compile MPMC normally, and then run:
```
cd mpmc_testing
python run_tests.py
```
More information about the test suite can be found in its [repository](https://github.com/LucianoLaratelli/mpmc_testing).

#### Wiki Link

https://github.com/mpmccode/mpmc/wiki

## Publications

* T. Pham, K.A. Forrest, A. Hogan, B. Tudor, K. McLaughlin, J.L. Belof, J. Eckert, and B. Space, "Understanding Hydrogen Sorption in In-soc-MOF: A Charged Metal-Organic Framework with Open-Metal Sites, Narrow Channels, and Counterions.", Crystal Growth & Design, 15(3), 1460-1471 (2015).

* Y. Zhang, B. Li, R. Krishna, Z. Wu, D. Ma, Z. Shi, T. Pham, K.A. Forrest, B. Space, and S. Ma, "Highly selective adsorption of ethylene over ethane in a MOF featuring the combination of open metal site and π-complexation.", Chemical Communications, 51, 2714-2717 (2015).

* T. Pham, K.A. Forrest, R. Banerjee, M.G. Orcajo, J. Eckert, and B. Space, "Understanding the H2 Sorption Trends in the M-MOF-74 Series (M= Mg, Ni, Co, Zn).", The Journal of Physical Chemistry C, 119(2), 1078-1090 (2014).

* T. Pham, K.A. Forrest, K. McDonald, and B. Space, "Modeling PCN-61 and PCN-66: Isostructural rht-Metal–Organic Frameworks with Distinct CO2 Sorption Mechanisms.", Crystal Growth & Design, 14(11), 5599-5607 (2014).

* T. Pham, K.A. Forrest, K. McLaughlin, J. Eckert, and B. Space, "Capturing the H2–Metal Interaction in Mg-MOF-74 Using Classical Polarization.", The Journal of Physical Chemistry C, 118(39), 22683-22690 (2014).

* T. Pham, K.A. Forrest, P. Georgiev, W. Lohstroh, D. Xue, A. Hogan, M. Eddaoudi, B. Space, and J. Eckert, "A high rotational barrier for physisorbed hydrogen in an fcu-metal–organic framework.", Chemical Communications, 50(91), 14109-14112 (2014).

* P. Nugent, T. Pham, K. McLaughlin, P.A. Georgiev, W. Lohstroh, J.P. Embs, M.J. Zaworotko, B. Space, and J. Eckert, "Dramatic effect of pore size reduction on the dynamics of hydrogen adsorbed in metal-organic materials.", Journal of Materials Chemistry A, 2, 13884-13891 (2014).

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
