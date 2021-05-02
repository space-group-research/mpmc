Introduction
************

.. image:: mpmc.png
  :width: 600
  :align: center

MPMC (Massively Parallel Monte Carlo) is an open-source Monte Carlo package primarily designed for the simulation of liquids, molecular interfaces and functionalized nanoscale materials. It was originally developed by `Jon Belof <http://people.llnl.gov/belof1>`_, includes contributions from Keith McLaughlin, Brant Tudor, Christian Cioce, Adam Hogan and Douglas Franz, and is currently maintained by the `Brian Space group <http://drbrian.space/>`_ in the Department of Chemistry at North Carolina State University. MPMC has been applied to the scientific research challenges of nanomaterials for clean energy, environmental sequestration, and molecular detection. Developed to run efficiently on the most powerful supercomputing platforms, MPMC can scale to very large numbers of CPUs or GPUs (with support provided for NVidia's CUDA).

Optimized for the study of nanoscale interfaces, MPMC supports many common intermolecular potentials including Lennard-Jones and damped dispersion paired with exponential repulsion, many-body polarization, coupled-dipole van der Waals, quantum rotational statistics, semi-classical quantum effects, advanced importance sampling methods relevant to fluids, and numerous tools for the development of intermolecular potentials.

Citing MPMC
===========

 | `MPMC and MCMD: Free High‐Performance Simulation Software for Atomistic Systems. <https://onlinelibrary.wiley.com/doi/full/10.1002/adts.201900113>`_
 | Franz, D. M.; Belof, J. L.; McLaughlin, K.; Cioce, C. R.; Tudor, B.; Hogan, A.; Laratelli, L.; Mulcair, M.; Mostrom, M.; Navas, A.; Stern, A. C.; Forrest, K. A.; Pham, T.; Space, B.
 | *Adv. Theory Simul.* **2019**, DOI: 10.1002/adts.201900113.

Supported Systems
=================

Currently a basic build of MPMC only requires a modern C compiler and CMake. Optionally MPMC may be configured to use `OpenMPI <https://www.open-mpi.org/>`_ or `CUDA <https://developer.nvidia.com/cuda-zone>`_. Configuring with Coupled-Dipole VDW or QM Rotation options requires `LAPACK <http://www.netlib.org/lapack/>`_.

MPMC supports compilation on Linux, macOS, and Windows; however MPMC is primarily tested on Linux and support is not guaranteed on other platforms.

Downloading MPMC
================

MPMC can be downloaded from the `releases page on GitHub <https://github.com/mpmccode/mpmc/releases>`_ or with the following command:

.. code-block:: none

    git clone https://github.com/mpmccode/mpmc

Building MPMC
=============

Once MPMC has been downloaded it may be quickly compiled with the following bash script:

.. code-block:: none

    bash compile.sh

Note that `CMake v2.8 <https://cmake.org/>`_ is required for compilation. More complete control may be exercised by invoking CMake directly:

.. code-block:: none

    mkdir build
    cd build
    cmake -DQM_ROTATION=OFF -DVDW=OFF -DMPI=OFF -DCUDA=OFF -DCMAKE_BUILD_TYPE=Release -Wno-dev ../
    make

Make sure to add MPMC to your path after compiling!

Running MPMC
============

MPMC accepts only one argument on the command line, the location of the MPMC input script:

.. code-block:: none

    mpmc mpmc.inp

Updating MPMC
=============

MPMC can be updated with the following command:

.. code-block:: none

    git pull
    
and then rebuilding as necessary.

MPMC Tutorials
==============

Example MPMC input scripts and PQRs are available in the tutorials_and_examples

MPMC Testing Suite
==================

An end-to-end test suite for MPMC is currently under development. If cloning MPMC anew, use the following to include the tests:

.. code-block:: none

    git clone https://github.com/mpmccode/mpmc --recurse-submodules

To clone the submodule into an existing MPMC installation, use this instead:

.. code-block:: none

    cd mpmc
    git submodule init
    git submodule update

To run the tests, make sure you have Python installed, compile MPMC normally, and then run:

.. code-block:: none

    cd mpmc_testing
    python run_tests.py

More information about the test suite can be found in its `repository <https://github.com/LucianoLaratelli/mpmc_testing>`_.

License
=======

MPMC is liscensed under the GNU GPL v3 license, a copy is located in the `root directory <https://github.com/mpmccode/mpmc/blob/master/LICENSE>`_.

