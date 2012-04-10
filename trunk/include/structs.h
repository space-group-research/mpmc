
/* complex value */
typedef struct _complex_t {

	double real;
	double imaginary;

} complex_t;

typedef struct _pair {
	int frozen;
	int rd_excluded, es_excluded;
	int attractive_only;
	double lrc;
//	double charge;
//	double polarizability;
	double last_volume; //what was the volume when we last calculated LRC? needed for NPT
	double epsilon;
	double sigma;
//	double gwp_alpha;
//	int gwp_spin;
	double r;
//	double d[3];
	double d_prev[3];
	int recalculate_energy;
	double rimg;
	double dimg[3];
//	double r2, r3, r5;		/* 2nd, 3rd and 5th powers for A matrix calc */
//	double r2img, r3img, r5img;
	double rd_energy, es_real_energy, es_self_intra_energy;
	struct _atom *atom;
	struct _molecule *molecule;
	struct _pair *next;
} pair_t;


typedef struct _atom {
	int    id;
        int    bond_id;
	char   atomtype[MAXLINE];
	int    frozen,
               adiabatic,
               spectre,
               target;
	double mass;
	double charge;
	double polarizability;
	double epsilon;
	double sigma;
	double omega; //for linear algebra VDW calculations
	double es_self_point_energy;
	double pos[3];
	double wrapped_pos[3];
	double ef_static[3];
	double ef_static_self[3];
	double ef_induced[3];
	double ef_induced_change[3];
	double mu[3];
	double old_mu[3];
	double new_mu[3];
	double dipole_rrms;
	double rank_metric;
	int gwp_spin;
	double gwp_alpha;
	int    site_neighbor_id; // dr fluctuations will be applied along the vector from this atom to the atom identified by this variable
	pair_t *pairs;
	struct _atom *next;
} atom_t;


typedef struct _molecule {
	int id;
	char moleculetype[MAXLINE];
	double mass;
	int frozen, adiabatic, spectre, target;
	double com[3];
	double wrapped_com[3];
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

typedef struct _vdw {
	char mtype[MAXLINE];
	double energy;
	struct _vdw * next;
} vdw_t;

typedef struct _pbc {
	double basis[3][3];		/* unit cell lattice (A) */
	double reciprocal_basis[3][3];	/* reciprocal space lattice (1/A) */
	double cutoff;			/* radial cutoff (A) */
	double volume;			/* unit cell volume (A^3) */
} pbc_t;

typedef struct _cavity {

	int occupancy;
	double pos[3];

} cavity_t;

typedef struct _histogram {
	int ***grid;
	int x_dim,y_dim,z_dim;
	double origin[3];
	double delta[3][3];
	int count[3];
	int n_data_points;
	int norm_total;
} histogram_t;

/* used in quantum rotation code */
typedef struct _spherical_harmonic {

	double legendre;        /* legendre polynomial part */
	double real;            /* real part of the spherical harmonic */
	double imaginary;       /* imaginary part of the spherical harmonic */

} spherical_harmonic_t;

typedef struct _nodestats {

	int accept, reject;
	int accept_insert, reject_insert;
	int accept_remove, reject_remove;
	int accept_displace, reject_displace;
	int accept_adiabatic, reject_adiabatic;
	int accept_spinflip, reject_spinflip;
	int accept_volume, reject_volume;

	double boltzmann_factor;
	double acceptance_rate;
	double acceptance_rate_insert;
	double acceptance_rate_remove;
	double acceptance_rate_displace;
	double acceptance_rate_adiabatic;
	double acceptance_rate_spinflip;
	double acceptance_rate_volume;
	double cavity_bias_probability;
	double polarization_iterations;

} nodestats_t;

typedef struct _avg_nodestats {

	double boltzmann_factor;
	double boltzmann_factor_sq;

	double acceptance_rate;
	double acceptance_rate_insert;
	double acceptance_rate_remove;
	double acceptance_rate_displace;
	double acceptance_rate_adiabatic;
	double acceptance_rate_spinflip;
	double acceptance_rate_volume;

	double cavity_bias_probability;
	double cavity_bias_probability_sq;

	double polarization_iterations;
	double polarization_iterations_sq;

} avg_nodestats_t;

typedef struct _observables {

	double energy;
	double coulombic_energy;
	double rd_energy;
	double polarization_energy;
	double vdw_energy; //for linear algebra VDW
	double dipole_rrms;
	double kinetic_energy;		/* for NVE */
	double temperature;		/* for NVE */
	double volume; /* for NPT */
	double N;
	double NU;
	double spin_ratio;		/* ortho:para spin ratio */

} observables_t;


typedef struct _avg_observables {

	/* these are uncorrelated averages of the observables */
	double energy;
	double energy_sq;
	double energy_error;

	/* needed for heat capacity error propagation */
	double energy_sq_sq;
	double energy_sq_error;

	double coulombic_energy;
	double coulombic_energy_sq;
	double coulombic_energy_error;

	double rd_energy;
	double rd_energy_sq;
	double rd_energy_error;

	double polarization_energy;
	double polarization_energy_sq;
	double polarization_energy_error;

	double vdw_energy; //for linear algebra VDW
	double vdw_energy_sq;
	double vdw_energy_error;

	double dipole_rrms;
	double dipole_rrms_sq;
	double dipole_rrms_error;

	/* for NVE MC */
	double kinetic_energy;
	double kinetic_energy_sq;
	double kinetic_energy_error;

	double temperature;
	double temperature_sq;
	double temperature_error;

	/* for NPT */
	double volume;
	double volume_sq;
	double volume_error;

	double N;
	double N_sq;
	double N_error;

	/* needed for qst */
	double NU;

	/* ortho:para spin ratio */
	double spin_ratio;
	double spin_ratio_sq;
	double spin_ratio_error;

	/* these quantities are node stats, not observables */
	double boltzmann_factor;
	double boltzmann_factor_sq;
	double boltzmann_factor_error;

	double acceptance_rate;
	double acceptance_rate_insert;
	double acceptance_rate_remove;
	double acceptance_rate_displace;
	double acceptance_rate_adiabatic;
	double acceptance_rate_spinflip;
	double acceptance_rate_volume;

	double cavity_bias_probability;
	double cavity_bias_probability_sq;
	double cavity_bias_probability_error;

	double polarization_iterations;
	double polarization_iterations_sq;
	double polarization_iterations_error;

	/* these quantities are based on averages; the error is not easily calculated */
	double qst;
	double heat_capacity;
	double compressibility;

	/* the error of these propogates */
	double density;
	double density_sq;
	double density_error;

	double pore_density;
	double pore_density_error;

	double percent_wt;
	double percent_wt_error;

	double percent_wt_me;
	double percent_wt_me_error;

	double excess_ratio;
	double excess_ratio_error;

} avg_observables_t;

typedef struct _grid {
	histogram_t *histogram;
	histogram_t *avg_histogram;
} grid_t;

/* begin mpi message struct */
typedef struct _message {
	/* observables */
	double energy;
	double coulombic_energy;
	double rd_energy;
	double polarization_energy;
	double kinetic_energy;
	double temperature;
	double volume;
	double N;
	double NU;
	double spin_ratio;

	/* avg_nodestats */
	double boltzmann_factor;
	double boltzmann_factor_sq;

	double acceptance_rate;
	double acceptance_rate_insert;
	double acceptance_rate_remove;
	double acceptance_rate_displace;
	double acceptance_rate_spinflip;
	double acceptance_rate_volume;

	double cavity_bias_probability;
	double cavity_bias_probability_sq;

	double polarization_iterations;
	double polarization_iterations_sq;

	/* histogram */
	int ***grid;
} message_t;

typedef struct _checkpoint {

	int movetype, biased_move;
	int N_atom, N_atom_prev;
	molecule_t *molecule_backup, *molecule_altered;
	molecule_t *head, *tail;
	observables_t *observables;

} checkpoint_t;

typedef struct _param {
	char atomtype[MAXLINE];
	int ntypes;
	double charge;
	double epsilon;
	double sigma;
	double omega;
	double dr;
	int axis;
	double last_charge;
	double last_epsilon;
	double last_sigma;
	double last_dr;
	double last_omega;
	struct _param *next;
} param_t;

typedef struct _file_pointers {

	FILE *fp_energy;
	FILE *fp_traj;
	FILE *fp_dipole;
	FILE *fp_field;
	FILE *fp_histogram;
	FILE *fp_frozen;
	
} file_pointers_t;



// For accomodating an arbitrary number of fit_input files
typedef union {
    int count;
    char *filename;
} intechar;

typedef struct _fileNode {
    intechar data;            // data.count is for use in the list head only & indicates the # of list elements
                              // data.filename is for use in the nodes
    struct _fileNode *next;
} fileNode_t;


// Metadata and a pointer to (r,E) (i.e. (distance, Energy)) point data
typedef struct _curveData {
    char   *id;                   // Id associated with the curve data.
    char   *filename;             // Filename associated with curve data
    int    weight;                // Integer representation of relative curve-weight
    double normalized_weight;     // Weight as a percent of total
    double alpha1, beta1, gamma1, // Orientation specs for two molecules
           alpha2, beta2, gamma2;
    int    nPoints;               // Number of data points
    double *r,                    // Array of molecule-molecule distances
           *input,                // Array of data points associated with each r-value
           *output,
           *global;
} curveData_t;

// used with qshift option to surf_fit
typedef struct _qshiftData {
	double qH2Q;
	double qH2G;
	double drH2Q;
} qshiftData_t;

typedef struct _surf_preserve_rotation {
	double alpha1, alpha2, beta1, beta2, gamma1, gamma2;
} surf_preserve_rotation;

typedef struct _system {

	////////////////////////////////////////////////////////////
	int tag; // Testing purposes only, erase when finished.  //
	//////////////////////////////////////////////////////////

	int ensemble;
	int gwp;
	int wolf;
	int surf_qshift_on;
	int surf_scale_epsilon_on, surf_scale_r_on, surf_scale_omega_on, surf_scale_sigma_on, surf_scale_q_on;
	int surf_weight_constant_on;
	int surf_preserve, surf_decomp;
	//some surf_preserve options
	surf_preserve_rotation * surf_preserve_rotation_on;
	long unsigned int seed;
	int numsteps, corrtime, step;
	double move_probability, rot_probability, insert_probability, adiabatic_probability, spinflip_probability, gwp_probability, volume_probability;
	double volume_change_factor, last_volume;
	int cavity_bias, cavity_grid_size;
	cavity_t ***cavity_grid;
	int cavities_open;
	double cavity_radius, cavity_volume;
	int cavity_autoreject, cavity_autoreject_absolute; //first is in terms of sigma and only applies to LJ; latter is in Angstroms and applies to all pairs
	double cavity_autoreject_scale;
	double temperature, pressure, fugacity, free_volume, total_energy, N;
	int spectre;
	double spectre_max_charge, spectre_max_target;
	int simulated_annealing;
	double simulated_annealing_schedule;
	int rd_only, rd_anharmonic, rd_lrc;
	int h2_fugacity, co2_fugacity, ch4_fugacity, n2_fugacity;
	double rd_anharmonic_k, rd_anharmonic_g;
	int feynman_hibbs, feynman_kleinert, feynman_hibbs_order;
	int sg, dreiding;
	int fvm, wpi, wpi_grid;
	int wrapall;
	double scale_charge, scale_rd;
	double ewald_alpha;
	int ewald_kmax;
	int polarization, polarizability_tensor;
	int polarvdw; //for linear algebra VDW
	int polar_iterative, polar_ewald, polar_zodid, polar_self, polar_palmo, polar_gs, polar_gs_ranked, polar_sor, polar_esor, polar_max_iter;
	double polar_gamma;
	double polar_damp, field_damp, polar_precision;
	int quantum_rotation, quantum_rotation_hindered;
	double quantum_rotation_B;
	double quantum_rotation_hindered_barrier;
	int quantum_rotation_level_max, quantum_rotation_l_max, quantum_rotation_theta_max, quantum_rotation_phi_max, quantum_rotation_sum;
	int quantum_vibration;
	int damp_type;
	int cuda;
	int independent_particle;

	double **A_matrix, **B_matrix, C_matrix[3][3];	/* A matrix, B matrix and polarizability tensor */
	char *pdb_input, *pdb_output, *pdb_restart, *traj_output, *energy_output;
	int read_pdb_box_on;
	char *dipole_output, *field_output, *histogram_output, *frozen_output;
	char *insert_input;

        // surface fitting parameters
        fileNode_t fit_input_list;
	double surf_min, surf_max, surf_inc, surf_ang;
	double surf_scale_epsilon, surf_scale_r, surf_scale_omega, surf_scale_sigma, surf_scale_q;
	double surf_quadrupole, surf_weight_constant;
        double fit_start_temp;
        double fit_max_energy;
        double fit_schedule;

	grid_t *grids;
	int calc_hist; /* flag to calculate a 3D histogram */
	double hist_resolution;
	int n_histogram_bins;

	double max_bondlength; /* threshold to bond (re:output files) */

	pbc_t *pbc;
	molecule_t *molecules;
	int        num_insertion_molecules;        // the number of elements found in both lists below:
	molecule_t    *insertion_molecules,        // linked list of molecules to be randomly inserted
	             **insertion_molecules_array;  // array providing direct access to elements of above list

	nodestats_t *nodestats;
	avg_nodestats_t *avg_nodestats;
	observables_t *observables;
	avg_observables_t *avg_observables;

	checkpoint_t *checkpoint;

	file_pointers_t file_pointers;

	vdw_t * vdw_eiso_info;

} system_t;

