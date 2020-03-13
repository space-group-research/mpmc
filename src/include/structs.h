#ifndef STRUCTS_H
#define STRUCTS_H

#include "defines.h"

#ifdef OPENCL
#include "CL/cl.h"
#endif

/* complex value */
typedef struct _complex_t {
    double real;
    double imaginary;
} complex_t;

typedef struct _ptemp {
    int *index;
    double *templist;
} ptemp_t;

typedef struct _pair {
    int frozen;  //are they both MOF atoms, for instance
    int rd_excluded, es_excluded;
    int attractive_only;
    int recalculate_energy;
    double lrc;               //LJ long-range correction
    double last_volume;       //what was the volume when we last calculated LRC? needed for NPT
    double epsilon, sigma;    //LJ
    double r, rimg, dimg[3];  //separation and separation with nearest image
    double d_prev[3];         //last known position
    double rd_energy, es_real_energy, es_self_intra_energy;
    double sigrep;
    double c6, c8, c10;
    struct _atom *atom;
    struct _molecule *molecule;
    struct _pair *next;
} pair_t;

typedef struct _atom {
    int id, bond_id;
    char atomtype[MAXLINE];
    int frozen, adiabatic, spectre, target;
    double mass, charge, polarizability, epsilon, sigma, omega;
    double c6, c8, c10, c9;
    double es_self_point_energy;
    double pos[3], wrapped_pos[3];  //absolute and wrapped (into main unit cell) position
    double ef_static[3], ef_static_self[3], ef_induced[3], ef_induced_change[3];
    double mu[3], old_mu[3], new_mu[3];
    double dipole_rrms;
    double rank_metric;
    int gwp_spin;
    double gwp_alpha;
    int site_neighbor_id;  // dr fluctuations will be applied along the vector from this atom to the atom identified by this variable
    pair_t *pairs;
    double lrc_self, last_volume;  // currently only used in disp_expansion.c
    struct _atom *next;

} atom_t;

typedef struct _molecule {
    int id;
    char moleculetype[MAXLINE];
    double mass;
    int frozen, adiabatic, spectre, target;
    double com[3], wrapped_com[3];  //center of mass
    double iCOM[3];                 // initial Center of Mass
    int nuclear_spin;
    double rot_partfunc_g, rot_partfunc_u, rot_partfunc;
#ifdef QM_ROTATION
    double *quantum_rotational_energies;
    complex_t **quantum_rotational_eigenvectors;
    int *quantum_rotational_eigensymmetry;
    double quantum_rotational_potential_grid[QUANTUM_ROTATION_GRID][QUANTUM_ROTATION_GRID];
#endif /* QM_ROTATION */
#ifdef XXX
    /* XXX - vib work in progress */
    double *quantum_vibrational_energies;
    complex_t **quantum_vibrational_eigenvectors;
    int *quantum_vibrational_eigensymmetry;
#endif /* XXX */
    atom_t *atoms;
    struct _molecule *next;
} molecule_t;

//stores vdw energies for each molecule within the coupled dipole model
typedef struct _vdw {
    char mtype[MAXLINE];
    double energy;
    struct _vdw *next;
} vdw_t;

//constants for peng_robinson equation of state
typedef struct _peng_robinson_constants {
    double Tc;
    double Pc;
    double w;
} peng_robinson_constants;

//info about periodic boundary and unit cell geometry
typedef struct _pbc {
    double basis[3][3];            /* unit cell lattice (A) */
    double reciprocal_basis[3][3]; /* reciprocal space lattice (1/A) */
    double cutoff;                 /* radial cutoff (A) */
    double volume;                 /* unit cell volume (A^3) */
} pbc_t;

typedef struct _cavity {
    int occupancy;
    double pos[3];
} cavity_t;

typedef struct _histogram {
    int ***grid;
    int x_dim, y_dim, z_dim;
    double origin[3];
    double delta[3][3];
    int count[3];
    int n_data_points;
    int norm_total;
} histogram_t;

/* used in quantum rotation code */
typedef struct _spherical_harmonic {
    double legendre;  /* legendre polynomial part */
    double real;      /* real part of the spherical harmonic */
    double imaginary; /* imaginary part of the spherical harmonic */

} spherical_harmonic_t;

typedef struct _nodestats {
    int accept, reject;
    int accept_insert, accept_remove, accept_displace, accept_adiabatic, accept_spinflip, accept_volume, accept_ptemp;
    int reject_insert, reject_remove, reject_displace, reject_adiabatic, reject_spinflip, reject_volume, reject_ptemp;
    double boltzmann_factor;
    double acceptance_rate;
    double acceptance_rate_insert, acceptance_rate_remove, acceptance_rate_displace;
    double acceptance_rate_adiabatic, acceptance_rate_spinflip, acceptance_rate_volume, acceptance_rate_ptemp;
    double cavity_bias_probability;
    double polarization_iterations;
} nodestats_t;

typedef struct _avg_nodestats {
    int counter;
    double boltzmann_factor, boltzmann_factor_sq;
    double acceptance_rate;
    double acceptance_rate_insert, acceptance_rate_remove, acceptance_rate_displace, acceptance_rate_ptemp;
    double acceptance_rate_adiabatic, acceptance_rate_spinflip, acceptance_rate_volume, acceprance_rate_ptemp;
    double cavity_bias_probability, cavity_bias_probability_sq;
    double polarization_iterations, polarization_iterations_sq;
} avg_nodestats_t;

typedef struct _observables {
    double energy;
    double coulombic_energy, rd_energy, polarization_energy, vdw_energy, three_body_energy;
    double dipole_rrms;
    double kinetic_energy; /* for NVE */
    double temperature;    /* for NVE */
    double volume;         /* for NPT */
    double N, NU;
    double spin_ratio;               /* ortho:para spin ratio */
    double frozen_mass, total_mass;  //updated in average.c
} observables_t;

typedef struct _avg_observables {
    /* these are uncorrelated averages of the observables */
    double energy, energy_sq, energy_error;
    double N, N_sq, N_error;
    double coulombic_energy, coulombic_energy_sq, coulombic_energy_error;
    double rd_energy, rd_energy_sq, rd_energy_error;
    double polarization_energy, polarization_energy_sq, polarization_energy_error;
    double vdw_energy, vdw_energy_sq, vdw_energy_error;
    double three_body_energy, three_body_energy_sq, three_body_energy_error;
    double dipole_rrms, dipole_rrms_sq, dipole_rrms_error;
    double density, density_sq, density_error;
    double pore_density, pore_density_error;
    double percent_wt, percent_wt_error;
    double percent_wt_me, percent_wt_me_error;
    double excess_ratio, excess_ratio_error;

    /* needed for heat capacity error propagation */
    double energy_sq_sq, energy_sq_error;

    /* for NVE MC */
    double kinetic_energy, kinetic_energy_sq, kinetic_energy_error;
    double temperature, temperature_sq, temperature_error;

    /* for NPT */
    double volume, volume_sq, volume_error;

    /* needed for qst */
    double NU;

    /* ortho:para spin ratio */
    double spin_ratio, spin_ratio_sq, spin_ratio_error;

    /* these quantities are node stats, not observables */
    double boltzmann_factor, boltzmann_factor_sq, boltzmann_factor_error;
    double cavity_bias_probability, cavity_bias_probability_sq, cavity_bias_probability_error;
    double polarization_iterations, polarization_iterations_sq, polarization_iterations_error;

    double acceptance_rate;
    double acceptance_rate_insert, acceptance_rate_remove, acceptance_rate_displace;
    double acceptance_rate_adiabatic, acceptance_rate_spinflip, acceptance_rate_volume, acceptance_rate_ptemp;

    /* these quantities are based on averages; the error is not easily calculated */
    double qst;
    double qst_nvt;
    double heat_capacity;
    double heat_capacity_error;
    double compressibility;
    double compressibility_error;

} avg_observables_t;

// used for storing global sorbate averages
typedef struct _sorbateAverages {
    double avgN, avgN_sq, avgN_err;                             // average sorbate count
    double percent_wt, percent_wt_sq, percent_wt_err;           // weight percent for this sorbate (sorb_mass / total_mass)
    double percent_wt_me, percent_wt_me_sq, percent_wt_me_err;  // weight percent for this sorbate (sorb_mass / frozen_mass)
    double excess_ratio, excess_ratio_sq, excess_ratio_err;     // excess adsorption ratio
    double pore_density, pore_density_sq, pore_density_err;     // mass / volume of pore
    double density, density_sq, density_err;                    // mass of sorbate / volume
    double selectivity, selectivity_err;                        // sorbate's selectivity ratio relative to all other sorbates in the insert list.
    struct _sorbateAverages *next;
} sorbateAverages_t;

// local sorbate data array
typedef struct _sorbateInfo {
    char id[16];  // identifying tag for the sorbate, e.g. CH4, CO2 or H2
    double mass;  // mass of this sorbate.
    int currN;    // sorbate count for the current step
    double percent_wt;
    double percent_wt_me;
    double excess_ratio;
    double pore_density;
    double density;
} sorbateInfo_t;

typedef struct _grid {
    histogram_t *histogram;
    histogram_t *avg_histogram;
} grid_t;

/* unused --  kmclaugh 2012 APR 16
// begin mpi message struct 
typedef struct _message {
	double energy;
	double coulombic_energy, rd_energy, polarization_energy, vdw_energy;
	double kinetic_energy, temperature, volume, N, NU, spin_ratio;
	double boltzmann_factor, boltzmann_factor_sq;
	double acceptance_rate;
	double acceptance_rate_insert, acceptance_rate_remove, acceptance_rate_displace;
	double acceptance_rate_spinflip, acceptance_rate_volume;
	double cavity_bias_probability, cavity_bias_probability_sq;
	double polarization_iterations, polarization_iterations_sq;
	int ***grid;
} message_t;
*/

typedef struct _checkpoint {
    int movetype, biased_move;
    int thole_N_atom;  //used for keeping track of thole matrix size (allocated)
    molecule_t *molecule_backup, *molecule_altered;
    molecule_t *head, *tail;
    observables_t *observables;
} checkpoint_t;

// used in surface fitting
typedef struct _param_type {
    char atomtype[MAXLINE];
    int ntypes, axis;
    double charge, epsilon, sigma, omega, dr, pol, c6, c8, c10;
    double last_charge, last_epsilon, last_sigma, last_dr, last_omega, last_pol, last_c6, last_c8, last_c10;
    struct _param_type *next;
} param_t;

typedef struct _param_global {
    double alpha;
    double last_alpha;
    param_t *type_params;
} param_g;

typedef struct _file_pointers {
    FILE *fp_energy;
    FILE *fp_energy_csv;
    FILE *fp_xyz;
    FILE *fp_field;
    FILE *fp_histogram;
    FILE *fp_frozen;
    //	FILE *fp_traj; //unused
    FILE *fp_traj_replay;
    FILE *fp_surf;
} file_pointers_t;

// For accomodating an arbitrary number of fit_input files
typedef union {
    int count;
    char *filename;
} intechar;

typedef struct _fileNode {
    intechar data;  // data.count is for use in the list head only & indicates the # of list elements
                    // data.filename is for use in the nodes
    struct _fileNode *next;
} fileNode_t;

// Metadata and a pointer to (r,E) (i.e. (distance, Energy)) point data
typedef struct _curveData {
    char *id;                      // Id associated with the curve data.
    char *filename;                // Filename associated with curve data
    int weight;                    // Integer representation of relative curve-weight
    double normalized_weight;      // Weight as a percent of total
    double alpha1, beta1, gamma1,  // Orientation specs for two molecules
        alpha2, beta2, gamma2;
    int nPoints;  // Number of data points
    double *r,    // Array of molecule-molecule distances
        *input,   // Array of data points associated with each r-value
        *output,
        *global;
} curveData_t;

// used with qshift option to surf_fit
typedef struct _qshiftData {
    double qH2Q, qH2G, drH2Q;
} qshiftData_t;

typedef struct _surf_preserve_rotation {
    double alpha1, alpha2, beta1, beta2, gamma1, gamma2;
} surf_preserve_rotation;

#ifdef OPENCL
typedef struct _ocl {
    /* kernel that performs a single iteration through the dipole field equations */
    cl_context context;
    cl_kernel kernel;
    cl_kernel palmo_echg;
    cl_kernel thole_estat;
    cl_kernel potential_reduction;
    cl_program program;
    cl_command_queue queue;
    cl_device_id *device_id;
    cl_platform_id *platforms;

} ocl_t;
#endif

typedef struct _system {
    int ensemble;
    int gwp;
    int cuda;
    int opencl;
#ifdef OPENCL
    ocl_t *ocl;
#endif

    //for manual specification of random seeds
    int preset_seeds_on;
    uint32_t preset_seeds;
    int rng_initialized;

    // is this a restart of a parallel job?
    int parallel_restarts;

    //surface fitting options
    int surf_fit_multi_configs;
    char *multi_fit_input;
    int surf_fit_arbitrary_configs;
    int surf_qshift_on, surf_scale_epsilon_on, surf_scale_r_on, surf_scale_omega_on, surf_scale_sigma_on, surf_scale_q_on, surf_scale_pol_on;
    int surf_weight_constant_on, surf_global_axis_on, surf_descent, surf_scale_alpha_on, surf_scale_c6_on, surf_scale_c8_on, surf_scale_c10_on;
    fileNode_t fit_input_list;
    double surf_scale_epsilon, surf_scale_r, surf_scale_omega, surf_scale_sigma, surf_scale_q, surf_scale_alpha, surf_scale_pol, surf_scale_c6, surf_scale_c8, surf_scale_c10;
    double surf_quadrupole, surf_weight_constant;
    double fit_start_temp, fit_max_energy, fit_schedule;
    int fit_boltzmann_weight;
    double fit_best_square_error;
    char **surf_do_not_fit_list;

    //surf options
    int surf_preserve, surf_decomp;
    double surf_min, surf_max, surf_inc, surf_ang;
    surf_preserve_rotation *surf_preserve_rotation_on;
    int surf_virial;
    char *virial_output;
    double *virial_coef;
    double virial_tmin, virial_tmax, virial_dt;
    int virial_npts;
    int ee_local;
    double range_eps, range_sig, step_eps, step_sig;

    //monte carlo controls
    int numsteps, corrtime, step;
    int ptemp_freq;
    double move_factor, rot_factor, insert_probability;
    double adiabatic_probability, spinflip_probability, gwp_probability, volume_probability;
    double volume_change_factor, last_volume;  //NPT
    //auto-reject options
    //first is in terms of sigma and only applies to LJ; latter is in Angstroms and applies to all pairs
    int cavity_autoreject_absolute;
    int count_autorejects;
    //parallel tempering options
    int parallel_tempering;
    double max_temperature;
    ptemp_t *ptemp;

    //cavity stuff
    int cavity_bias, cavity_grid_size;
    cavity_t ***cavity_grid;
    int cavities_open;
    double cavity_radius, cavity_volume, cavity_autoreject_scale, cavity_autoreject_repulsion;

    //spectre
    int spectre;
    double spectre_max_charge, spectre_max_target;

    //simulated annealing
    int simulated_annealing, simulated_annealing_linear;
    double simulated_annealing_schedule;
    double simulated_annealing_target;

    //force-field options
    int rd_only, rd_anharmonic;
    double rd_anharmonic_k, rd_anharmonic_g;
    int sg, dreiding, waldmanhagler, lj_buffered_14_7, halgren_mixing, c6_mixing, disp_expansion;
    int extrapolate_disp_coeffs, damp_dispersion, disp_expansion_mbvdw;
    int axilrod_teller, midzuno_kihara_approx;
    //es_options
    int wolf;
    double ewald_alpha, polar_ewald_alpha;
    int ewald_alpha_set, polar_ewald_alpha_set;
    int ewald_kmax;
    //thole options
    int polarization, polarvdw, polarizability_tensor;
    int cdvdw_exp_repulsion, cdvdw_sig_repulsion, cdvdw_9th_repulsion;
    int iter_success;  //flag set when iterative solver fails to converge (when polar_precision is used)
    int polar_iterative, polar_ewald, polar_ewald_full, polar_zodid, polar_palmo, polar_rrms;
    int polar_gs, polar_gs_ranked, polar_sor, polar_esor, polar_max_iter, polar_wolf, polar_wolf_full, polar_wolf_alpha_lookup;
    double polar_wolf_alpha, polar_gamma, polar_damp, field_damp, polar_precision;
    int damp_type;
    double **A_matrix, **B_matrix, C_matrix[3][3]; /* A matrix, B matrix and polarizability tensor */
    vdw_t *vdw_eiso_info;                          //keeps track of molecule vdw self energies
    double *polar_wolf_alpha_table, polar_wolf_alpha_lookup_cutoff;
    int polar_wolf_alpha_table_max;  //stores the total size of the array

    // energy-corrections
    int feynman_hibbs, feynman_kleinert, feynman_hibbs_order;
    int vdw_fh_2be;  //2BE method for polarvdw
    int rd_lrc, rd_crystal, rd_crystal_order;

    // uvt fugacity functions
    int h2_fugacity, co2_fugacity, ch4_fugacity, n2_fugacity, user_fugacities;

    // i/o options
    int wrapall;
    char *job_name;  // (CRC)
    char *pqr_input, *pqr_output, *pqr_restart, *traj_input, *traj_output, *energy_output, *energy_output_csv, *surf_output, *xyz_output;
    int read_pqr_box_on;   //read box basis from pqr
    int long_output;       // prints extended (%11.6f) coordinates
    int surf_print_level;  // sets the amount of output (1-6) that correspond to the nested loops in surface.c
    char *dipole_output, *field_output, *histogram_output, *frozen_output;
    char *insert_input;
    double max_bondlength; /* threshold to bond (re:output files) */
    // insertions from a separate linked list
    int num_insertion_molecules;             // the number of elements found in both lists below:
    molecule_t *insertion_molecules;         // linked list of molecules to be randomly inserted
    molecule_t **insertion_molecules_array;  // array providing direct access to elements of above list

    // quantum rotation stuff
    int quantum_rotation, quantum_rotation_hindered;
    double quantum_rotation_B;
    double quantum_rotation_hindered_barrier;
    int quantum_rotation_level_max, quantum_rotation_l_max, quantum_rotation_theta_max, quantum_rotation_phi_max, quantum_rotation_sum;
    int quantum_vibration;

    // histogram stuff
    grid_t *grids;
    int calc_hist; /* flag to calculate a 3D histogram */
    double hist_resolution;
    int n_histogram_bins;

    //observables
    double temperature, pressure, free_volume, total_energy, N;
    double *fugacities;
    int fugacitiesCount;

    //atom array
    int natoms;
    atom_t **atom_array;
    molecule_t **molecule_array;

    //replay option
    int calc_pressure;
    double calc_pressure_dv;

    //misc
    double scale_charge;
    int independent_particle;

    pbc_t *pbc;
    molecule_t *molecules;

    nodestats_t *nodestats;
    avg_nodestats_t *avg_nodestats;
    observables_t *observables;
    avg_observables_t *avg_observables;

    // Linked list head that will keep track of separate average-observables
    // for each sorbate in the system.
    int sorbateCount;                  // Number of sorbates in the system.
    sorbateInfo_t *sorbateInfo;        //stores an array of sorbate Info
    int sorbateInsert;                 //which sorbate was last inserted
    sorbateAverages_t *sorbateGlobal;  //where the global average is stored

    checkpoint_t *checkpoint;
    file_pointers_t file_pointers;

} system_t;

#endif
