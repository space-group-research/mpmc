Input File Specifications
*************************

The bulk of MPMC runs require two plain text files, the input script and the PQR data file. The input script contains all of the commands defining the type of run to be performed as well as details about the state point, control of the potentials, control of the input/output, etc. The PQR file contains the atomic positions of the system as well as force field parameters such as partial charges and van der Waals coefficients on each atom.

MPMC Input Script
=================

The MPMC input script contains a series of commands, usually of the form :code:`[command name] [on|off|value]`, and one per line. Comments may be included by beginning a line with :code:`!` or :code:`#`. Whitespace is ignored and the order of the commands is not important as the entire input script is read and then the simulation is started. A minimal MPMC input script contains the ensemble to simulate in, the temperature (and possibly pressure), the number of steps, and the output frequency. As an example, a minimal input script for a :math:`\mu VT` simulation of H\ :sub:`2` sorption in MOF-5 is provided below.

.. code-block:: none

    job_name MOF5+BSS

    ensemble uvt
    
    temperature 298.0
    pressure    1.0

    numsteps 100
    corrtime 4

    insert_probability 0.667

    pqr_input input.pqr
    abcbasis  25.669 25.669 25.669 90 90 90


The full list of commands is available in :doc:`commands`.

PQR File
=============

PQR files used by MPMC contain additional columns compared to standard .pqr or .pdb files to support inclusion of the force field parameters. 

.. code-block:: none

    1     2    3          4           5    6         7  8  9  10    11      12              13          14        15     16         17  18  19
    ATOM  ID#  Element    MolecLabel  M/F  MolecID   x  y  z  mass  charge  polarizability  LJ epsilon  LJ sigma  Omega  GWP alpha  C6  C8  C10


1: ATOM verbatim


2: Atom ID, starting from 1 to N\ :sub:`atoms`


3: Element, can include additional numbers/letters in the case of multiple atom types (e.g. "ZN", "C1", "O2", "H2G", etc)


4: Molecular label, doesn't have to be unique (e.g. "MOF" or "H2")


5: M = Movable, F = Frozen (a MOF, e.g., is normally frozen, while sorbates are movable.)


6: Molecule ID, starting from 1 to N\ :sub:`molecules`


7-9: X, Y, Z coordinates in Angstroms


10: Mass of atom in amu


11: Partial charge in e


12: Polarizability in Angstrom\ :sup:`3`


13: :math:`\epsilon` (in K) for Lennard-Jones simulations or :math:`\beta` (in Angstrom\ :sup:`-1`) for PHAHST simulations


14: :math:`\sigma` (in Angstrom) for Lennard-Jones simulations or :math:`\rho` (in Angstrom) for PHAHST simulations


15: :math:`\omega` (in a.u.) for many-body van der Waals interactions


16: :math:`\alpha` for gaussian wave packet Coulombic interactions (normally not needed)


17-19: Dispersion coefficients for PHAHST simulations


For typical Lennard-Jones simulations columns 15-19 are not needed and if omitted will default to 0. An excerpt of the PQR file from the second tutorial, BSSP H\ :sub:`2` sorption in MOF-5, is provided below as an example.

.. code-block:: none

    ATOM      1 ZN   MOF F   1       7.568   5.314  -7.516  65.3900   1.8530  0.16000 62.39930  2.46200
    ATOM      2 ZN   MOF F   1       5.335  -5.287  -5.283  65.3900   1.8530  0.16000 62.39930  2.46200
    ATOM      3 ZN   MOF F   1       5.335   7.547   7.551  65.3900   1.8530  0.16000 62.39930  2.46200
    ATOM      4 ZN   MOF F   1       5.335  -7.520   7.551  65.3900   1.8530  0.16000 62.39930  2.46200
    ATOM      5 ZN   MOF F   1      -5.266   7.547  -7.516  65.3900   1.8530  0.16000 62.39930  2.46200
    [...]
    ATOM    425 H2G  H2  M    2      0.000   0.000   0.000  0.00000 -0.74640  0.69380 12.76532  3.15528
    ATOM    426 H2E  H2  M    2      0.371   0.000   0.000  1.00800  0.37320  0.00044  0.00000  0.00000
    ATOM    427 H2E  H2  M    2     -0.371  -0.000   0.000  1.00800  0.37320  0.00044  0.00000  0.00000
    ATOM    428 H2N  H2  M    2      0.363   0.000   0.000  0.00000  0.00000  0.00000  2.16726  2.37031
    ATOM    429 H2N  H2  M    2     -0.363  -0.000   0.000  0.00000  0.00000  0.00000  2.16726  2.37031


Surface Fitting Files
=====================

Surf_multi_fit Files
--------------------
