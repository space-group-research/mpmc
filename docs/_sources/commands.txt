Commands
********

Main Commands
=============

.. csv-table::
    :header: "Command","Description"
    :widths: 20,40

    "cuda [on|off]", "Turns on/off nVIDIA CUDA **(default = off)**"
    "preset_seeds [int] [int] [int] [int]", "Initialize the random number generator using the specified seed values"
    "numsteps [int]", "Total number of steps (moves) for Monte Carlo simulations"
    "corrtime [int]", "Sets the correlation time in steps. Each corrtime averages and trajectory outputs are produced"
    "temperature [double]", "Sets the (initial) temperature"
    "pressure [double]", "Sets the initial pressure for NPT and uVT simulations. In uVT the pressure is converted to a fugacity via a fugacity option or assuming the ideal gas limit (P==f)"

Monte Carlo Options
-------------------


Ensembles and Ensemble Specific Commands
----------------------------------------

NVE Options
-----------

uVT Options
-----------

Replay Options
--------------

Input / Output Commands
=======================

Potential Commands
==================

Mixing Rules
------------

Ewald/Wolf Options
------------------

Polarization Options
--------------------

PHAHST Options
--------------

Annealing / Tempering Commands
==============================

Quantum Rotation Commands
=========================

Surface Fitting Commands
========================

External Tools
==============


