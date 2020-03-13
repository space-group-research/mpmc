#ifndef FXN_PROTOTYPES_H
#define FXN_PROTOTYPES_H

/* energy */
double energy(system_t *);
double energy_no_observables(system_t *);
double cavity_absolute_check(system_t *);
double lj(system_t *);
double lj_nopbc(system_t *);
double exp_repulsion(system_t *);
double exp_repulsion_nopbc(system_t *);
double dreiding(system_t *);
double dreiding_nopbc(molecule_t *);
double lj_buffered_14_7(system_t *);
double lj_buffered_14_7_nopbc(system_t *);
double disp_expansion_lrc(const system_t *, pair_t *, const double);
double disp_expansion_lrc_self(const system_t *, atom_t *, const double);
double disp_expansion(system_t *);
double disp_expansion_nopbc(system_t *);
double axilrod_teller(system_t *system);
double factorial(int);
double tt_damping(int, double);
void countN(system_t *);
void update_com(molecule_t *);
void flag_all_pairs(system_t *);
void pair_exclusions(system_t *, molecule_t *, molecule_t *, atom_t *, atom_t *, pair_t *);
void minimum_image(system_t *, atom_t *, atom_t *, pair_t *);
void pairs(system_t *);
void setup_pairs(system_t *);
void update_pairs_insert(system_t *);
void update_pairs_remove(system_t *);
void unupdate_pairs_insert(system_t *);
void unupdate_pairs_remove(system_t *);
double pbc_cutoff(pbc_t *);
double pbc_volume(pbc_t *);
void pbc(system_t *);
double sg(system_t *);
double sg_nopbc(molecule_t *);
double coulombic(system_t *);
double coulombic_wolf(system_t *);
double coulombic_real(system_t *);
double coulombic_reciprocal(system_t *);
double coulombic_background(system_t *);
double coulombic_nopbc(molecule_t *);
double coulombic_real_gwp(system_t *);
double coulombic_reciprocal_gwp(system_t *);
double coulombic_background_gwp(system_t *);
double coulombic_nopbc_gwp(system_t *);
double coulombic_kinetic_gwp(system_t *);
double morse_energy(double, double, double, double);
double anharmonic(system_t *);
double anharmonic_energy(double, double, double);
double anharmonic_fk(double, double, double, double, double);
double anharmonic_fh_second_order(double, double, double, double, double);
double anharmonic_fh_fourth_order(double, double, double, double, double);
double h2_bond_energy(double r);
int bessik(float, float, float *, float *, float *, float *); /* NR function */
double besselK(double, double);
void rebuild_arrays(system_t *);  //builds atom and molecule arrays for the current config

/* io */
void write_observables_csv(FILE *, system_t *, observables_t *, double);
void write_molecules_xyz(system_t *, FILE *);  //L
void update_sorbate_info(system_t *);
int safe_atof(char *, double *);
int safe_atoi(char *a, int *i);
int safe_atou(char *a, uint32_t *i);
int safe_atol(char *a, long unsigned int *l);

int check_system(system_t *);
system_t *read_config(char *);
int setup_simulation_box(FILE *, system_t *);
int check_config(system_t *);
molecule_t *read_molecules(FILE *, system_t *);
int read_pqr_box(FILE *, system_t *);
system_t *setup_system(char *);
char *make_filename(char *, int);
void error(char *);
void output(char *);
void clear_nodestats(nodestats_t *);
void clear_node_averages(avg_nodestats_t *);
void clear_observables(observables_t *);
void clear_sorbate_averages(sorbateAverages_t *, int);
void clear_root_averages(avg_observables_t *);
void clear_avg_nodestats(system_t *);
void calc_system_mass(system_t *);
void track_ar(nodestats_t *);
void update_nodestats(nodestats_t *, avg_nodestats_t *);
void update_root_averages(system_t *, observables_t *, avg_observables_t *);
void update_root_sorb_averages(system_t *, sorbateInfo_t *);
void update_root_nodestats(system_t *, avg_nodestats_t *, avg_observables_t *);
int write_performance(int, system_t *);
int print_observables(system_t *);
int write_averages(system_t *);
int write_molecules(system_t *, FILE *);
int write_molecules_wrapper(system_t *, char *);
void write_observables(FILE *, system_t *, observables_t *, double);
void write_dipole(system_t *);
void write_field(system_t *);
void write_states(system_t *);
void write_surface_traj(FILE *, system_t *);
int wrapall(molecule_t *, pbc_t *);
void spectre_wrapall(system_t *);
int open_files(system_t *);
int open_surf_traj_file(system_t *);
void close_files(system_t *);
curveData_t *readFitInputFiles(system_t *, int);
molecule_t *read_insertion_molecules(system_t *);
void count_sorbates(system_t *);
void write_virial_output(system_t *, double, double, double);
#ifdef OPENCL
ocl_t *setup_ocl();
#endif

/* main */
void die(int);
void usage(char *);
#ifdef QM_ROTATION
void free_rotational(system_t *);
#endif /* QM_ROTATION */
void free_pairs(molecule_t *);
void free_atoms(molecule_t *);
void free_molecule(system_t *, molecule_t *);
void free_molecules(molecule_t *);
void free_averages(system_t *system);
void free_matrices(system_t *system);
void free_cavity_grid(system_t *system);
void cleanup(system_t *);
void terminate_handler(int, system_t *);
int memnullcheck(void *, int, int, char *);
int filecheck(void *, char *, int);
void free_all_molecules(system_t *, molecule_t *);
void free_all_pairs(system_t *);

/* mc */
void temper_system(system_t *, double);
void enumerate_particles(system_t *);
void boltzmann_factor(system_t *, double, double, double);
void register_accept(system_t *);
void register_reject(system_t *);
int mc(system_t *);
molecule_t *copy_molecule(system_t *, molecule_t *);
void translate(system_t *system, molecule_t *, pbc_t *, double);
void rotate(system_t *system, molecule_t *, pbc_t *, double);
void displace(system_t *system, molecule_t *, pbc_t *, double, double);
void displace_1D(system_t *, molecule_t *, double);
void displace_gwp(system_t *system, molecule_t *, double);
void spectre_displace(system_t *, molecule_t *, double, double, double);
void spectre_charge_renormalize(system_t *);
void make_move(system_t *);
void checkpoint(system_t *);
void restore(system_t *);
double surface_energy(system_t *, int);
void molecule_rotate_euler(molecule_t *, double, double, double, int);
void molecule_rotate_quaternion(molecule_t *, double, double, double, int);
int surface_dimer_geometry(system_t *, double, double, double, double, double, double, double, int);
int surface_dimer_parameters(system_t *, param_g *);
int surface_dimer_geometry_virial(system_t *, double, double, double, double, double, double, int);
void reset_molecule_position(molecule_t *);
void surface_curve(system_t *, double, double, double, double *);
int surface(system_t *);
int surface_fit(system_t *);
int surface_virial(system_t *);
void setup_cavity_grid(system_t *);
void cavity_volume(system_t *);
void cavity_probability(system_t *);
void cavity_update_grid(system_t *);
void qshift_do(system_t *, qshiftData_t *, double, double);
double calcquadrupole(system_t *);
void volume_change(system_t *);
void revert_volume_change(system_t *);

/*surface_fit*/
void free_all_mem(int, curveData_t *, param_g *, qshiftData_t *, double *);
void apply_new_parameters(param_g *);
void surface_curve(system_t *, double, double, double, double *);
double error_calc(system_t *, int, int, curveData_t *, double);
int alloc_curves(int, int, curveData_t *);
void output_pqrs(system_t *, int, curveData_t *);
void output_params(double, double, param_g *);
param_g *record_params(system_t *);
void surf_perturb(system_t *, double, qshiftData_t *, param_g *);
void output_fit(int, int, curveData_t *, double, double *);
void get_curves(system_t *, int, curveData_t *, double, double, double);
void revert_parameters(system_t *, param_g *);
void new_global_min(system_t *, int, int, curveData_t *);

/* other "ensembles" */
int replay_trajectory(system_t *);
int calculate_te(system_t *);

/* polarization */
double polar(system_t *);
void thole_amatrix(system_t *);
void thole_bmatrix(system_t *);
void thole_bmatrix_dipoles(system_t *);
void thole_polarizability_tensor(system_t *);
void thole_field(system_t *);
void thole_field_wolf(system_t *);
void thole_field_nopbc(system_t *);
void thole_field_real(system_t *);
void thole_field_recip(system_t *);
void thole_field_self(system_t *);
double *polar_wolf_alpha_lookup_init(system_t *);
double polar_wolf_alpha_getval(system_t *, double);
int thole_iterative(system_t *);
void invert_matrix(int, double **, double **);
int countNatoms(system_t *);
void thole_resize_matrices(system_t *);
void print_matrix(int N, double **matrix);
void ewald_estatic(system_t *);
void ewald_full(system_t *);
void calc_dipole_rrms(system_t *);
int are_we_done_yet(system_t *, int);

/* polarization - CUDA */
#ifdef CUDA
void *polar_cuda(void *);
#endif /* CUDA */

/* linear algebra - VDW */
double vdw(system_t *);
void free_vdw_eiso(vdw_t *);

/* pimc */
int pimc(system_t *);

/* histogram */
int histogram(system_t *);
int cart2frac(double *, double *, system_t *);
int frac2cart(double *, double *, system_t *);
void setup_histogram(system_t *);
void allocate_histogram_grid(system_t *);
void write_histogram(FILE *, int ***, system_t *);
void zero_grid(int ***, system_t *);
void population_histogram(system_t *);
void mpi_copy_histogram_to_sendbuffer(char *, int ***, system_t *);
void mpi_copy_rcv_histogram_to_data(char *, int ***, system_t *);
void update_root_histogram(system_t *);
void write_histogram(FILE *, int ***, system_t *);

/* dxwrite */
void write_frozen(FILE *fp_frozen, system_t *system);

#ifdef QM_ROTATION
/* quantum rotation */
void quantum_system_rotational_energies(system_t *);
void quantum_rotational_energies(system_t *, molecule_t *, int, int);
void quantum_rotational_grid(system_t *, molecule_t *);
complex_t **rotational_hamiltonian(system_t *, molecule_t *, int, int);
int determine_rotational_eigensymmetry(molecule_t *, int, int);
double rotational_basis(int, int, int, double, double);
double rotational_potential(system_t *, molecule_t *, double, double);
double hindered_potential(double);
double rotational_integrate(system_t *, molecule_t *, int, int, int, int, int);
#endif /* QM_ROTATION */

#ifdef DEBUG
void test_pairs(molecule_t *);
void test_molecule(molecule_t *);
void test_list(molecule_t *);
void test_cavity_grid(system_t *);
void test_lj(system_t *);
void test_q(system_t *);
#endif  /* DEBUG */
#endif  // FXN_PROTOTYPES_H

//fugacity functions
double h2_fugacity(double, double);
double h2_fugacity_back(double, double);
double h2_comp_back(double, double);
double h2_fugacity_shaw(double, double);
double h2_fugacity_zhou(double, double);
double ch4_fugacity(double, double);
double ch4_fugacity_back(double, double);
double ch4_comp_back(double, double);
double ch4_fugacity_PR(double, double);
double n2_fugacity(double, double);
double n2_fugacity_back(double, double);
double n2_comp_back(double, double);
double n2_fugacity_PR(double, double);
double n2_fugacity_zhou(double, double);
double co2_fugacity(double, double);

double get_rand(system_t *system);

//useful math calls
double dddotprod(double *, double *);
double didotprod(double *, int *);
int iidotprod(int *, int *);
double min(double a, double b);
