Introduction
************

MPMC (Massively Parallel Monte Carlo) is an open-source Monte Carlo package primarily designed for the simulation of liquids, molecular interfaces and functionalized nanoscale materials. It was originally developed by Jon Belof and is now maintained by a group of researchers (Keith McLaughlin, Brant Tudor, Christian Cioce, Adam Hogan, Douglas Franz and Brian Space) in the Department of Chemistry and SMMARTT Materials Research Center at the University of South Florida. MPMC has been applied to the scientific research challenges of nanomaterials for clean energy, environmental sequestration and molecular detection. Developed to run efficiently on the most powerful supercomputing platforms, MPMC can scale to very large numbers of CPUs or GPUs (with support provided for NVidia's CUDA and the OpenCL architecture).

Optimized for the study of nanoscale interfaces, MPMC supports simulation of Coulomb and Lennard-Jones systems, many-body polarization, coupled-dipole van der Waals, quantum rotational statistics, semi-classical quantum effects, advanced importance sampling methods relevant to fluids, and numerous tools for the development of intermolecular potentials.

Citing MPMC
===========



Supported Systems
=================

Currently a basic build of MPMC only requires a modern C compiler and CMake. Optionally MPMC may be configured to use OpenMPI, CUDA or OpenCL. Configuring with Coupled-Dipole VDW or QM Rotation requires LAPACK.

MPMC supports compilation on Linux, macOS, and Windows; however MPMC is primarily tested on Linux and support is not guaranteed on other platforms.

Downloading MPMC
================



Building MPMC
=============



Running MPMC
============



MPMC Testing Suite
==================



License
=======

MPMC is liscensed under the GNU GPL v3 license, a copy is located in the root directory.
