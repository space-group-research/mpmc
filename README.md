## MPMC

![MPMC](docssrc/mpmc.png)

MPMC (Massively Parallel Monte Carlo) is an open-source Monte Carlo package primarily designed for the simulation of liquids, molecular interfaces and functionalized nanoscale materials. It was originally developed by [Jon Belof](http://people.llnl.gov/belof1), includes contributions from Kieth McLaughlin, Brant Tudor, Christian Cioce, Adam Hogan and Douglas Franz, and is currently maintained by the [Brian Space group](http://drbrian.space/) in the Department of Chemistry at North Carolina State University. MPMC has been applied to the scientific research challenges of nanomaterials for clean energy, environmental sequestration, and molecular detection. Developed to run efficiently on the most powerful supercomputing platforms, MPMC can scale to very large numbers of CPUs or GPUs (with support provided for NVidia's CUDA and the OpenCL architecture).

Optimized for the study of nanoscale interfaces, MPMC supports many common intermolecular potentials including Lennard-Jones and damped dispersion paired with exponential repulsion, many-body polarization, coupled-dipole van der Waals, quantum rotational statistics, semi-classical quantum effects, advanced importance sampling methods relevant to fluids, and numerous tools for the development of intermolecular potentials.

## Getting Started

#### Libraries

Currently a basic build of MPMC only requires a modern C compiler and CMake. Optionally MPMC may be configured to use [OpenMPI](https://www.open-mpi.org/), [CUDA](https://developer.nvidia.com/cuda-zone) or [OpenCL](https://www.khronos.org/opencl/). Configuring with Coupled-Dipole VDW or QM Rotation requires [LAPACK](http://www.netlib.org/lapack/).

#### Supported Platforms

MPMC supports compilation on Linux, macOS, and Windows; however MPMC is primarily tested on Linux and support is not guaranteed on other platforms.

#### Downloading MPMC

MPMC can be downloaded with the following command.

```
git clone https://github.com/mpmccode/mpmc
```

#### Compiling MPMC

Once MPMC has been downloaded it may be compiled with the following command.

```
bash compile.sh
```

#### Documentation Link

Additional documentation, including tutorials, detailed command lists, and PHAST/PHAHST models are available at the following link.

http://mpmccode.github.io/mpmc/


## Selected Publications



* [Next-Generation Accurate, Transferable, and Polarizable Potentials for Material Simulations.](https://pubs.acs.org/doi/10.1021/acs.jctc.0c00837)\
Hogan, A.; Space, B. *J. Chem. Theory Comput.* **2020**, 16 (12), 7632-7644.

* [On the Mechanism of Hydrogen Storage in a Metal-Organic Framework Material.](https://pubs.acs.org/doi/abs/10.1021/ja0737164)\
Belof, J. L.; Stern, A. C.; Eddaoudi, M.; Space, B.,\
*J. Am. Chem. Soc.* **2007**, 129 (49), 15202-15210.

* [Introduction of π-Complexation into Porous Aromatic Framework for Highly Selective Adsorption of Ethylene over Ethane.](https://pubs.acs.org/doi/abs/10.1021/ja502119z)\
Li, B.; Zhang, Y.; Krishna, R.; Yao, K.; Han, Y.; Wu, Z.; Ma, D.; Shi, Z.; Pham, T.; Space, B.; Liu, J.; Thallapally, P. K.; Liu, J.; Chrzanowski, M.; Ma, S.\
*J. Am. Chem. Soc.* **2014**, 136 (24), 8654-8660.

* [A Robust Molecular Porous Material with High CO2 Uptake and Selectivity](https://pubs.acs.org/doi/10.1021/ja4054948)\
Nugent, P.S.; Rhodus, V.L.; Pham, T.; Forrest, K.; Wojtas, L.; Space, B.; Zaworotko, M.J.\
*J. Am. Chem. Soc.* **2013**, 135 (30), 10950–10953. 

* [Capturing the H2–Metal Interaction in Mg-MOF-74 Using Classical Polarization.](https://pubs.acs.org/doi/10.1021/jp508249c)\
Pham, T.; Forrest, K. A.; McLaughlin, K.; Eckert, J.; Space, B.\
*J. Phys. Chem. C* **2014**, 118, 22683–22690.

* [Efficient calculation of many-body induced electrostatics in molecular systems.](https://aip.scitation.org/doi/full/10.1063/1.4829144)\
McLaughlin, K.; Cioce, C. R.; Pham, T.; Belof, J. L.; Space, B.\
*J. Chem. Phys.* **2013**, 139, 184112. 

* [A Polarizable and Transferable PHAST N2 Potential For Use in Materials Simulation.](https://pubs.acs.org/doi/abs/10.1021/ct400526a)\
Cioce, C. R.; McLaughlin, K.; Belof, J. L.; Space B.\
*J. Chem. Theory Comput.* **2013**, 9, 5550–5557. 

* [A Polarizable and Transferable PHAST CO2 Potential For Materials Simulation.](https://pubs.acs.org/doi/10.1021/ct400549q)\
Mullen, A. L.; Pham, T.; Forrest, K. A.; Cioce, C. R.; McLaughlin, K.; Space, B.\
*J. Chem. Theory Comput.* **2013**, 9, 5421–5429. 

* [A molecular H2 potential for heterogeneous simulations including polarization and many-body van der Waals interactions.](https://aip.scitation.org/doi/10.1063/1.4717705)\
McLaughlin, K.; Cioce, C. R.; Belof, J. L.; Space, B.,\
*J. Chem. Phys.* **2012**, 136 (19).

