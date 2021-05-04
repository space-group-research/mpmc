Commands
********

Main Commands
=============

.. csv-table::
    :header: "Command","Description"
    :widths: 20,40

    "job_name [string]", "Base job name for output files. **(required)**"
    "seed [int]", "Initialize the random number generator using the specified seed values. **(default = randomized on startup based on current time)**"
    "numsteps [int]", "Total number of steps (moves) for Monte Carlo simulations. **(required for Monte Carlo)**"
    "corrtime [int]", "Sets the correlation time in steps. Every corrtime steps averages and trajectory outputs are produced. **(required for Monte Carlo)**"
    "temperature [double]", "Sets the (initial) temperature. **(required for Monte Carlo)**"
    "pressure [double]", "Sets the initial pressure for NPT and uVT simulations. In uVT the pressure is converted to a fugacity via a fugacity option or assuming the ideal gas limit (P==f). **(required for NPT/uVT Monte Carlo)**"
    "cuda [on|off]", "Turns on/off nVIDIA CUDA for the calculation of the polarization equations. **(default = off)**"

General Monte Carlo Options
---------------------------

.. csv-table::
    :header: "Command","Description"
    :widths: 20,40

    "move_factor [double]", "Scale factor for translational moves in a fraction of the simulation box (1.0 = full box length). **(default = 1.0)**"
    "rot_factor [double]", "Scale factor for rotational moves in a fraction of a complete rotation (1.0 = full 360 degrees). **(default = 1.0)**"



Ensembles and Ensemble Specific Commands
----------------------------------------

.. csv-table::
    :header: "Command","Description"
    :widths: 20,40

    "ensemble [ENSEMBLE]", "Chooses the ensemble to simulate in. ENSEMBLE may be one of 'npt', 'nve', 'nvt', 'uvt', 'surf', 'surf_fit', 'total_energy' or 'replay'. **(required)**"

\

.. csv-table::
    :header: "ENSEMBLE","Description"
    :widths: 20,40

    "npt", "Monte Carlo simulation in the isobaric-isothermal ensemble."
    "nve", "Monte Carlo simulation in the micro-canonical ensemble."
    "nvt", "Monte Carlo simulation in the canonical ensemble."
    "uvt", "Monte Carlo simulation in the grand canonical ensemble."
    "surf", "Potential energy surface scan for some specified orientation."
    "surf_fit", "Perform surface fitting via simulated annealing over parameter space."
    "total_energy", "Calculate the single-point energy for a given system."
    "replay", "Recalculate observables over a set of configurations read from an input trajectory."

NVE Options
-----------

.. csv-table::
    :header: "Command","Description"
    :widths: 20,40

    "total_energy", "Sets total energy for NVE Monte Carlo. **(required for NVE Monte Carlo)**"

uVT Options
-----------

.. csv-table::
    :header: "Command","Description"
    :widths: 20,40

    "pressure [double]", "Sets the initial pressure for NPT and uVT simulations. In uVT the pressure is converted to a fugacity via a fugacity option or assuming the ideal gas limit (P==f). **(required for NPT/uVT Monte Carlo)**"
    "h2_fugacity [on|off]", "Used for converting pressure to fugacity of H2 uVT simulations. (Zhou/Shaw/BACK depending on state point) **(default = off)**"
    "co2_fugacity [on|off]", "Used for converting pressure to fugacity of CO2 uVT simulations. (Peng-Robinson) **(default = off)**"
    "ch4_fugacity [on|off]", "Used for converting pressure to fugacity of CH4 uVT simulations. (Peng-Robinson/BACK depending on state point) **(default = off)**"
    "n2_fugacity [on|off]", "Used for converting pressure to fugacity of N2 uVT simulations. (Zhou/Peng-Robinson/BACK depending on state point) **(default = off)**"
    "user_fugacities [double] ([double]) ([double]) (...)", "Specifies the fugacities for species in insert_input. Accepts up to eight arguments."
    "insert_probability [double]", "In a uVT simulation, the probability to randomly choose an insertion/deletion move over a translational/rotational move. **(default = 0.0)**"
    "free_volume [double]", "Used for statistics calculations in MOF sorption simulations."
    "cavity_bias [on|off]", "Use cavity bias for insertions in uVT Monte Carlo. **(default = off)**"
    "cavity_grid [int]", "Number of grid points in each dimension used in cavity_bias."
    "cavity_radius [double]", "Radius of sphere used to determine if a cavity grid point is unoccupied."

NPT Options
-----------

.. csv-table::
    :header: "Command","Description"
    :widths: 20,40

    "pressure [double]", "Sets the initial pressure for NPT and uVT simulations. In uVT the pressure is converted to a fugacity via a fugacity option or assuming the ideal gas limit (P==f). **(required for NPT/uVT Monte Carlo)**"
    "volume_probability [double]", "In a NPT simulation, the probability to randomly choose a volume change move over a translational/rotational move. **(default = 0.0)**"
    "volume_change_factor [double]", "In a NPT simulation, the scale factor for the change of cell parameters in a fraction of the current cell parameters. **(default = 0.25)**"

Replay Options
--------------

.. csv-table::
    :header: "Command","Description"
    :widths: 20,40

    "calc_pressure [on|off]", "Perform Frenkel pressure calculation over the uvt/nvt input trajectory. **(default = off)**"
    "calc_pressure_dv [double]", "Size of volume moves to perform in the pressure calculation."

Input / Output Commands
=======================

.. csv-table::
    :header: "Command","Description"
    :widths: 20,40

    "pqr_input [filename]", "Specifies input pqr file. **(default based on job_name)**"
    "pqr_output [filename]", "Specifies filename for writing final output pqr file(s). Clobbers existing file(s). **(default based on job_name)**"
    "pqr_restart [filename]", "Specifies filename for writing restart pqr file(s). Clobbers existing file(s). **(default based on job_name)**"
    "traj_input [filename]", "Specifies input trajectory file for 'ensemble replay'."
    "traj_output [filename]", "Specifies filename for writing trajectory pqr file(s). Clobbers existing file(s). **(default based on job_name)**"
    "energy_output [filename]", "Specifies filename for writing observables log. Clobbers existing file(s). **(default based on job_name)**"
    "energy_output_csv [filename]", "Specifies filename for writing observables log in csv format. Clobbers existing file(s)."
    "xyz_output [filename]", "Specifies filename for writing trajectory in xyz format. Clobbers existing file(s)."
    "pop_histogram [on|off]", "Turns on population histogram. **(default = off)**"
    "pop_histogram_output [filename]", "Specifies filename for writing popular histogram. Clobbers existing file(s). **(default = histogram.dx)**"
    "dipole_output [filename]", "Specifies filename for writing induced dipole data. Clobbers existing file(s)."
    "field_output [filename]", "Specifies filename for writing total electrostatic field for each molecule. Clobbers existing file(s)."
    "frozen_output [filename]", "Specifies filename for writing frozen atoms in dx format. Clobbers existing file(s)."
    "insert_input [filename]", "Specifies filename for reading molecules for performing insertions in uVT simulations."
    "parallel_restarts [on|off]", "Forces each MPMC thread to restart from its own pqr_restart file. **(default = off)**"
    "long_output [on|off]", "Prints additional sigfigs for atom xyz info in output pqr's. **(default = off, unless box has a dimension >= 100 Å)**"
    "read_pqr_box [on|off]", "Reads simulation box dimensions from pqr input file. **(default = off)**"
    "wrapall [on|off]", "Wraps atoms back into the simulation box on output. **(default = on)**"
    "basis1 [double] [double] [double]", "Specifies the basis vector's x-, y- and z- components."
    "basis2 [double] [double] [double]", "Specifies the basis vector's x-, y- and z- components."
    "basis3 [double] [double] [double]", "Specifies the basis vector's x-, y- and z- components."
    "[abcbasis|carbasis] [double] [double] [double] [double] [double] [double]", "Specifies the basis vectors of the unit-cell by a, b, c, alpha, beta, gamma."

Potential Commands
==================

.. csv-table::
    :header: "Command","Description"
    :widths: 20,40

    "cavity_autoreject_absolute [on|off]", "Automatically rejects any monte carlo move which would put two sites (not on the same molecule) that are too close. **(default = off)**"
    "cavity_autoreject_scale [double]", "Sets threshold (distance in Angstroms) for triggering cavity_autoreject_absolute."
    "cavity_autoreject_repulsion [double]", "Automatically rejects any monte carlo move where the repulsive energy is greater than the value input. Currently only implemented in combination with disp_expansion."
    "feynman_hibbs [on|off]", "Turns on Feynman-Hibbs quantum corrections. **(default = off)**"
    "feynman_hibbs_order [2|4]", "Specifies highest-order Feynman-Hibbs terms to use."
    "pbc_cutoff [double]", "Override the default cutoff distance for interactions. **(default = half the shortest simulation box dimension)**"
    "scale_charge [double]", "Scales the charges on all frozen atoms."
    "rd_lrc [on|off]", "Turns on long-range corrections to repulsion/dispersion energies via integration from r_cutoff to infinity. **(default = on)**"
    "rd_only [on|off]", "Only calculate repulsion/dispersion energies. (excludes coupled dipole vdW) **(default = off)**"
    "sg [on|off]", "Silvera-Goldman potential (hard coded, see src/energy/sg.c for details). **(default = off)**"
    "dreiding [on|off]", "Dreiding potential. (see src/energy/dreiding.c for details) **(default = off)**"

Lennard-Jones Mixing Rules
--------------------------

.. csv-table::
    :header: "Command","Description"
    :widths: 20,40

    "waldmanhagler [on|off]", "Use Waldman-Hagler mixing rules for Lennard-Jones RD rather than the default (Lorentz-Berthelot). **(default = off)**"
    "halgren_mixing [on|off]", "Use Halgren mixing rules for Lennard-Jones RD rather than the default (Lorentz-Berthelot). **(default = off)**"
    "c6_mixing [on|off]", "Use the known C6 mixing rule to calculate the Lennard-Jones epsilon. The Lennard-Jones sigma is calculated using arithmetic mean. **(default = off)**"

Ewald/Wolf Options
------------------

.. csv-table::
    :header: "Command","Description"
    :widths: 20,40

    "wolf [on|off]", "Calculates permanent electrostatics via wolf method. If off permanent electrostatics are handled via the Ewald summation method. **(default = off)**"
    "polar_ewald_full [on|off]", "Full ewald polarization (induced and static) for periodic systems. (used in conjunction to polarization/polarvdw). **(default = off)**"
    "polar_ewald [on|off]", "Partial ewald polarization (static-only) for periodic systems. (used in conjunction to polarization/polarvdw). **(default = off)**"
    "polar_ewald_alpha [int]", "Sets alpha/damping parameter for polar ewald calculation."
    "polar_wolf_full [on|off]", "Full wolf polarization (induced and static) for periodic systems. (used in conjunction to polarization/polarvdw). **(default = off)**"
    "polar_wolf [on|off]", "Partial wolf polarization (static-only) for periodic systems. (used in conjunction to polarization/polarvdw). **(default = off)**"
    "[polar_wolf_damp|polar_wolf_alpha] [int]", "Sets alpha/damping parameter for polar wolf calculation."
    "polar_wolf_lookup [on|off]", "Uses a lookup table for calculation of erfc's in wolf calculation. Grid size probably needs to be tweaked in the source. **(default = off)**"
    "ewald_alpha [double]", "Overrides default alpha for ewald and wolf permanent electrostatics and polar_ewald. **(default = 3.5/pbc_cutoff)**"
    "ewald_kmax [int]", "Sets the maximum k-vectors to include in ewald sums for permanent electrostatics and polarization. **(default = 7)**"

Polarization Options
--------------------

.. csv-table::
    :header: "Command","Description"
    :widths: 20,40

    "polarization [on|off]", "Turns on Thole-Applequist polarization. **(default = off)**"
    "polar_damp_type [off|none|linear|exponential]", "Type of polarization damping. (off=none)"
    "polar_damp [double]", "Polarization exponential damping constant (to help avoid polarization catastrophe). **(required if polar_damp_type != off)**"
    "polar_ewald [on|off]", "Calculate induced polarization via ewald summation. **(default = off)**"
    "polarizability_tensor [on|off]", "Prints the molecular polarizability tensor for the system. **(default = off)**"
    "polar_zodid [on|off]", "Calculates polarization energy via zeroth-order iteration. **(default = off)**"
    "polar_iterative [on|off]", "Full iterative method for calculation polarization energy. **(default = off)**"
    "polar_palmo [on|off]", "Iterative polar correction due to Kim Palmo. **(default = off)**"
    "polar_gs [on|off]", "Gauss-Seidel smoothing for iterative polarization. **(default = off)**"
    "polar_gs_ranked [on|off]", "Ranked Gauss-Seidel smoothing for iterative polarization. **(default = off)**"
    "polar_sor [on|off]", "(Linear??) polarization overrelaxation. **(default = off)**"
    "polar_esor [on|off]", "Exponential polarization overrelaxation. **(default = off)**"
    "polar_gamma [double]", "Polarization overrelaxation constant."
    "polar_precision [double]", "Terminate polarization iterative solver when all dipole fluctuations are within this tolerance. **(either polar_precision or polar_max_iter required if polarization = on)**"
    "polar_max_iter [int]", "Terminate polarization iterative solver after a fixed number of iterations. **(either polar_precision or polar_max_iter required if polarization = on)**"
    "polar_self [on|off]", "Include molecular self-induction. **(default = off)**"
    "polar_rrms [on|off]", "Calculate root-mean-square fluctuation in dipoles elements during iterative solution. **(default = off)**"

PHAHST Options
--------------

.. csv-table::
    :header: "Command","Description"
    :widths: 20,40

    "[phahst|disp_expansion] [on|off]", "Activates a RD potential similar to the Tang-Toennies potential. :math:`E_{rd} = -\frac{C6}{r^6}-\frac{C8}{r^8}-\frac{C10}{r^{10}}+596.725194095/ \epsilon * \mathrm{exp}(- \epsilon * ( r - \sigma))`. **(default = off)**"
    "damp_dispersion [on|off]", "Damps the PHAHST dispersion interaction according to Tang and Toennies's incomplete gamma functions. **(default = on)**"
    "extrapolate_disp_coeffs [on|off]", "Extrapolates C10 from C6 and C8. **(default = off)**"

Coupled-Dipole Van der Waals Options
------------------------------------

.. csv-table::
    :header: "Command","Description"
    :widths: 20,40

    "[cdvdw|polarvdw] [on|off|evects|comp]", "Turns on coupled-dipole method van der Waals. Evects prints eigenvectors and comp prints a comparision to a two-body decomposition. Also activates polarization. **(default = off)**"
    "vdw_fh_2be [on|off] ", "Uses two-body expansion for calculation of Feynman Hibbs in coupled-dipole vdW calculations. **(default = off)**"
    "cdvdw_9th_repulsion [on|off]", "Use 9th power mixing rule :math:`rep_{ij} = (rep_{ii}^{1/9}+rep_{jj}^{1/9})^9` for repulsion interactions (used in conjunction with coupled-dipole vdW) **(default = off)**"
    "cdvdw_sig_rep [on|off]", "Calculate repulsion using :math:`\frac{3}{2} \hbar w_i w_j \alpha_i \alpha_j / ( w_i + w_j ) * sig6` with WH mixing for sigma) (used in conjunction with coupled-dipole vdW) **(default = off)**"
    "cdvdw_exp_rep [on|off]", "Uses exponential repulsion :math:`\sigma * \mathrm{exp}(-\frac{r}{2 \epsilon})`, using some mixing rule I found somewhere -- see source code. **(default = off)**"

Miscellaneous Options
---------------------

None of these are guaranteed to work.\

.. csv-table::
    :header: "Command","Description"
    :widths: 20,40

    "spectre [on|off]", "??? **(default = off)**"
    "spectre_max_charge [double]", ""
    "spectre_max_target [double]", ""
    "rd_anharmonic [on|off]", "1-dimensional anharmonic spring potential. **(default = off)**"
    "rd_anharmonic_k [double]", "Harmonic term (order 2)."
    "rd_anharmonic_g [double]", "Anharmonic term (order 4)."
    "feynman_kleinert [on|off]", "Iterative Feynman-Kleinert correction for anharmonic bond potential. **(default = off)**"

Annealing / Tempering Commands
==============================

.. csv-table::
    :header: "Command","Description"
    :widths: 20,40

    "simulated_annealing [on|off]", "Turns on simulated annealing for MC simulations. **(default = off)**"
    "simulated_annealing_schedule [double]", "(Exponential) decay constant for the temperature in a simulated annealing MC simulation."
    "simulated_annealing_target [double]", "Target temperature in a simulated annealing MC simulation."
    "simulated_annealing_linear [on|off]", "Sets a linear ramp throughout the entire simulation instead of exponential decay. **(default = off)**"

Parallel Tempering Options
--------------------------

.. csv-table::
    :header: "Command","Description"
    :widths: 20,40

    "parallel_tempering [on|off]", "Turns on parallel tempering for Monte Carlo simulations. **(default = off)**"
    "ptemp_freq [int]", "How often to perform bath swaps when performing Monte Carlo with parallel tempering. **(default = 20)**"
    "max_temperature [double]", "Sets the temperature for the hottest bath in a MC simulation with parallel tempering."


Quantum Rotation Commands
=========================

.. csv-table::
    :header: "Command","Description"
    :widths: 20,40

    "quantum_rotation [on|off]", "Enables quantum rotational eigenspectrum calculation. **(default = off)**"
    "quantum_rotation_hindered [on|off]", "Calculates the rotational energy levels using the hindered potential :math:`\textrm{sin}^2 \theta`. **(default = off)**"
    "quantum_rotation_hindered_barrier [double]", ""
    "quantum_rotation_B [double]", "Sets the rotational constant. For H2, it is 85.35060622 Kelvin."
    "quantum_rotation_level_max [int]", "Number of rotational energy levels to solve for (Equal to (l_max + 1)2). **(default = 36)**"
    "quantum_rotation_l_max [int]", "Number of rotational energy levels to solve for (Equal to (l_max + 1)2). **(default = 36)**"
    "quantum_rotation_sum [int]", "Number of rotational energy levels to sum over. **(default = 10)**"
    "quantum_vibration [on|off]", "**(default = off)**"

Surface Fitting Commands
========================

.. csv-table::
    :header: "Command","Description"
    :widths: 20,40

    "fit_start_temp [double]", "Intial temperature for parameter annealing during surface fitting. **(default = 50000)**"
    "fit_schedule [double]", "Temperature (exponential) decay constant for parameter annealing during surface fitting. **(default = 0.999)**"
    "fit_max_energy [double]", "Maximum energy values to be considered during surin two weeksface fitting. **(default = 2000)**"
    "fit_input [file]", "Specifies fit input file. Call multiple times to specify multiple fit geometries."
    "surf_descent [on|off]", "Only accept parameter moves that lower the square error (rather than a Monte Carlo approach)."
    "surf_weight_constant [double]", "Exponential weighting factor used in surface fitting to prioritize fitting at lower potential energies. **(default = 0.5)**"
    "surf_scale_q [double]", "Magnitude of charge fluctuations during surface fitting **(default = 0)**"
    "surf_scale_r [double]", "Magnitude of position fluctuations of non-H2E/H2Q sites (when that site exists in multiples/doesn't apply to CoM site) **(default = 0.001)**"
    "surf_scale_epsilon [double]", "Magnitude of epsilon fluctations to sites with non-zero epsilon. **(default = 1.0)**"
    "surf_scale_sigma [double]", "Magntiude of sigma fluctuations to sites with non-zero sigma. **(default = 0.1)**"
    "surf_scale_omega [double]", "Magnitude of omega fluctuations to sites with non-zero omega. **(default = 0.001)**"
    "surf_scale_pol [double]", "Magnitude of alpha (polarizabilities) fluctuations to sites with non-zero alpha. **(default = 0)**"
    "surf_qshift [on|off]", "Adjusts position of H2Q sites while adjusting charges of H2Q and H2G to remain charge neutral and conserve quadrupole. **(default = off)**"
    "surf_global_axis [on|off]", "Use quaternions to rotate the molecules about the cartesian axes during surface fitting (as opposed to the local axes which is the default) **(default = off)**"
    "surf_scale_pol [double]", "Magnitude of polarizability fluctuations to sites with non-zero polarizability. **(default = 0)**"
    "surf_scale_c6 [double]", "Magnitude of C6 fluctuations to sites with non-zero C6. **(default = 0)**"
    "surf_scale_c8 [double]", "Magnitude of C8 fluctuations to sites with non-zero C8. **(default = 0)**"
    "surf_scale_c10 [double]", "Magnitude of C10 fluctuations to sites with non-zero C10. **(default = 0)**"
    "surf_multi_fit [on|off]", "Fit using a separate trajectory file (multi_fit_input) specifying *ab initio* energies and atomic positions. More flexible than the default algorithm but requires explicit preparation of atomic coordinates. **(default = off)**"
    "multi_fit_input [string]", "Filename for surf_multi_fit."
    "surf_do_not_first_list [string] ([string] ...)", "List of atom types to ignore during fitting. Currently only implemented with surf_multi_fit."

Surface Scan Options
--------------------

.. csv-table::
    :header: "Command","Description"
    :widths: 20,40

    "surf_decomp [on|off]", "Decompose the total energy into it's components (i.e. electrostatic, polarization, Lennard-Jones, CPM-vdW) **(default = off)**"
    "surf_min [double]", "Minimum center-of-mass to center-of-mass pair separation in the generation of the potential energy surface. **(default = 0.25)**"
    "surf_max [double]", "Maximum center-of-mass to center-of-mass pair separation in the generation of the potential energy surface. **(default = 25.0)**"
    "surf_inc [double]", "Increment for center-of-mass to center-of-mass pair separation in the generation of the potential energy surface. **(default = 0.25)**"
    "surf_ang [double]", "Angular increment (in rads) for generation of the isotropic potential energy surface."
    "surf_preserve [on|off]", "Preserve orientation while calculating surface curves. **(default = off)**"
    "surf_preserve_rotation [double] [double] [double] [double] [double] [double]", "Perform rotation to molecules prior to calculating surface curves. (a1, b1, c1, a2, b2, c2)"
    "surf_print_level [1-6]", "Verbosity of surface trajectory data written to surf_output. **(default = 3)**"
    "surf_output", "Output file for surface curves."

External Tools
==============


.. csv-table::
    :header: "Tool","Description"
    :widths: 20,40

    "traj_pqr2pdb", "Designed a workaround for viewing long_output PQR trajectories in VMD! This is useful for ultra-dense systems, or if you're just precision obsessed like me :D. The shell script (traj_pqr2pdb.sh) will convert your PQR to a PDB trajectory, which you can then read into VMD. You'll do so like this:
    
    
    [user@machine]$ ./traj_pqr2pdb.sh INPUT.pqr > OUTPUT.pdb [user@machine]$ vmd OUTPUT.pdb -e traj_pqr2pdb.vmd
    
    
    More specifically, the shell script will temporarily reduce the precision of your atomic coordinates (.6f) to .3f, and store that removed precision in another column of the file. The PDB that is created is a readable trajectory file by VMD. The secondary script that you load into VMD along with the PDB will reassign the atomic coordinates of each atom, in each frame, restoring the precision.
    
    
    Note: Although the file mode bits should be preserved, the shell script should be executable but the vmd script should not be."
