/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

/* this file is rather sloppy, but got the job done */

#include <mc.h>
#include <quaternion.h>

/* calculate the energy for the surface scan - same as function energy() but without observables */
double surface_energy(system_t *system, int energy_type) {
    double rd_energy, coulombic_energy, polar_energy, vdw_energy;

    /* zero the initial values */
    rd_energy = 0;
    coulombic_energy = 0;
    polar_energy = 0;
    vdw_energy = 0;

    /* get the pairwise terms necessary for the energy calculation */
    pairs(system);

    // energy_type = ENERGY_TOTAL
    // all other params are 0

    switch (energy_type) {
        case ENERGY_TOTAL:
            if (!(system->sg || system->rd_only)) coulombic_energy = coulombic_nopbc(system->molecules);
            if (system->polarization) polar_energy = polar(system);
            if (system->polarvdw) vdw_energy = vdw(system);
            if (system->sg)
                rd_energy = sg_nopbc(system->molecules);
            else if (system->cdvdw_exp_repulsion)
                rd_energy = exp_repulsion(system);
            else if (system->dreiding)
                rd_energy = dreiding_nopbc(system->molecules);
            else if (system->lj_buffered_14_7)
                rd_energy = lj_buffered_14_7_nopbc(system);
            else if (system->disp_expansion)
                rd_energy = disp_expansion_nopbc(system);
            else
                rd_energy = lj_nopbc(system);

            break;
        case ENERGY_ES:
            if (!(system->sg || system->rd_only)) coulombic_energy = coulombic_nopbc(system->molecules);
            break;
        case ENERGY_RD:
            if (system->sg) {
                rd_energy = sg_nopbc(system->molecules);
            } else if (system->dreiding) {
                rd_energy = dreiding_nopbc(system->molecules);
            } else if (system->lj_buffered_14_7) {
                rd_energy = lj_buffered_14_7_nopbc(system);
            } else if (system->disp_expansion) {
                rd_energy = disp_expansion_nopbc(system);
            } else {
                rd_energy = lj_nopbc(system);
            }
            break;
        case ENERGY_POLAR:
            if (system->polarization) polar_energy = polar(system);
            break;
        case ENERGY_VDW:
            if (system->polarvdw) vdw_energy = vdw(system);
            break;
    }

    /* return the total potential energy */
    return rd_energy + coulombic_energy + polar_energy + vdw_energy;
}

/* rotate a molecule about three Euler angles - same as normal function but for a single mol and without random angles */
void molecule_rotate_euler(molecule_t *molecule, double alpha, double beta, double gamma, int reverse_flag) {
    atom_t *atom_ptr;
    double rotation_matrix[3][3];
    double com[3];
    int i, ii, n;
    double *new_coord_array;

    //reverse the rotation?
    if (reverse_flag) {
        //switch around alpha and gamma
        double temp;
        temp = alpha;
        alpha = gamma;
        gamma = temp;
        alpha *= -1.;
        beta *= -1.;
        gamma *= -1.;
    }

    /* count the number of atoms in a molecule, and allocate new coords array */
    for (atom_ptr = molecule->atoms, n = 0; atom_ptr; atom_ptr = atom_ptr->next)
        ++n;
    new_coord_array = calloc(n * 3, sizeof(double));
    memnullcheck(new_coord_array, n * 3 * sizeof(double), __LINE__ - 1, __FILE__);

    /* save the com coordinate */
    com[0] = molecule->com[0];
    com[1] = molecule->com[1];
    com[2] = molecule->com[2];

    /* translate the molecule to the origin */
    for (atom_ptr = molecule->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
        atom_ptr->pos[0] -= com[0];
        atom_ptr->pos[1] -= com[1];
        atom_ptr->pos[2] -= com[2];
    }

    /* construct the 3D rotation matrix */
    rotation_matrix[0][0] = cos(gamma) * cos(beta) * cos(alpha) - sin(gamma) * sin(alpha);
    rotation_matrix[0][1] = sin(gamma) * cos(beta) * cos(alpha) + cos(gamma) * sin(alpha);
    rotation_matrix[0][2] = -sin(beta) * cos(alpha);
    rotation_matrix[1][0] = -cos(gamma) * cos(beta) * sin(alpha) - sin(gamma) * cos(alpha);
    rotation_matrix[1][1] = -sin(gamma) * cos(beta) * sin(alpha) + cos(gamma) * cos(alpha);
    rotation_matrix[1][2] = sin(beta) * sin(alpha);
    rotation_matrix[2][0] = cos(gamma) * sin(beta);
    rotation_matrix[2][1] = sin(gamma) * sin(beta);
    rotation_matrix[2][2] = cos(beta);

    /* matrix multiply */
    for (atom_ptr = molecule->atoms, i = 0; atom_ptr; atom_ptr = atom_ptr->next, i++) {
        ii = i * 3;
        new_coord_array[ii + 0] = rotation_matrix[0][0] * atom_ptr->pos[0] + rotation_matrix[0][1] * atom_ptr->pos[1] + rotation_matrix[0][2] * atom_ptr->pos[2];
        new_coord_array[ii + 1] = rotation_matrix[1][0] * atom_ptr->pos[0] + rotation_matrix[1][1] * atom_ptr->pos[1] + rotation_matrix[1][2] * atom_ptr->pos[2];
        new_coord_array[ii + 2] = rotation_matrix[2][0] * atom_ptr->pos[0] + rotation_matrix[2][1] * atom_ptr->pos[1] + rotation_matrix[2][2] * atom_ptr->pos[2];
    }

    /* set the new coordinates and then translate back from the origin */
    for (atom_ptr = molecule->atoms, i = 0; atom_ptr; atom_ptr = atom_ptr->next, i++) {
        ii = i * 3;
        atom_ptr->pos[0] = new_coord_array[ii + 0];
        atom_ptr->pos[1] = new_coord_array[ii + 1];
        atom_ptr->pos[2] = new_coord_array[ii + 2];

        atom_ptr->pos[0] += com[0];
        atom_ptr->pos[1] += com[1];
        atom_ptr->pos[2] += com[2];
    }

    /* free our temporary array */
    free(new_coord_array);
}

/* rotate a molecule about three Euler angles - same as normal function but for a single mol and without random angles */
/* changing it so it rotates around the three principal axises - Adam Hogan */
void molecule_rotate_quaternion(molecule_t *molecule, double alpha, double beta, double gamma, int reverse_flag) {
    atom_t *atom_ptr;
    double com[3];
    int i, ii, n;
    double *new_coord_array;

    /* count the number of atoms in a molecule, and allocate new coords array */
    for (atom_ptr = molecule->atoms, n = 0; atom_ptr; atom_ptr = atom_ptr->next)
        ++n;
    new_coord_array = calloc(n * 3, sizeof(double));
    memnullcheck(new_coord_array, n * 3 * sizeof(double), __LINE__ - 1, __FILE__);

    /* save the com coordinate */
    com[0] = molecule->com[0];
    com[1] = molecule->com[1];
    com[2] = molecule->com[2];

    /* translate the molecule to the origin */
    for (atom_ptr = molecule->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
        atom_ptr->pos[0] -= com[0];
        atom_ptr->pos[1] -= com[1];
        atom_ptr->pos[2] -= com[2];
    }

    /* construct the rotation quaternions */
    struct quaternion first, second, third, total, total_conjugate, position, answer;
    //total = quaternion with no rotation
    quaternion_construct_xyzw(&total, 0., 0., 0., 1.);
    //first = rotation around z axis
    quaternion_construct_axis_angle_degree(&first, 0., 0., 1., alpha);
    //second = rotation around y axis
    quaternion_construct_axis_angle_degree(&second, 0., 1., 0., beta);
    //third = rotation around z axis
    quaternion_construct_axis_angle_degree(&third, 1., 0., 0., gamma);
    //combine all rotations into total
    quaternion_multiplication(&first, &total, &total);
    quaternion_multiplication(&second, &total, &total);
    quaternion_multiplication(&third, &total, &total);
    //if we are undoing a rotation we need to conjugate total
    if (reverse_flag) quaternion_conjugate(&total, &total);
    quaternion_conjugate(&total, &total_conjugate);

    for (atom_ptr = molecule->atoms, i = 0; atom_ptr; atom_ptr = atom_ptr->next, i++) {
        ii = i * 3;
        //position = position vector
        quaternion_construct_xyzw(&position, atom_ptr->pos[0], atom_ptr->pos[1], atom_ptr->pos[2], 0.);

        //answer = total*(position*total_conjugate)
        quaternion_multiplication(&position, &total_conjugate, &answer);
        quaternion_multiplication(&total, &answer, &answer);

        //set the new coords
        new_coord_array[ii + 0] = answer.x;
        new_coord_array[ii + 1] = answer.y;
        new_coord_array[ii + 2] = answer.z;
    }

    /* set the new coordinates and then translate back from the origin */
    for (atom_ptr = molecule->atoms, i = 0; atom_ptr; atom_ptr = atom_ptr->next, i++) {
        ii = i * 3;
        atom_ptr->pos[0] = new_coord_array[ii + 0];
        atom_ptr->pos[1] = new_coord_array[ii + 1];
        atom_ptr->pos[2] = new_coord_array[ii + 2];

        atom_ptr->pos[0] += com[0];
        atom_ptr->pos[1] += com[1];
        atom_ptr->pos[2] += com[2];
    }

    /* free our temporary array */
    free(new_coord_array);
}

/* calculate the isotropic potential energy surface of a molecule */
int surface(system_t *system) {
    int ao, bo, go, am, bm, gm, surf_print; /* for surf_output */
    int avg_counter;
    double avg_factor;
    double r, pe_es, pe_rd, pe_polar, pe_vdw, pe_total;
    double pe_total_avg, pe_es_avg, pe_rd_avg, pe_polar_avg, pe_vdw_avg;
    double alpha_origin, beta_origin, gamma_origin;
    double alpha_move, beta_move, gamma_move;

    /* open surface trajectory file */
    if (system->surf_output) {
        if (open_surf_traj_file(system) < 0)
            error(
                "Surface: could not open surface trajectory file\n");
    }

    /* output the potential energy curve */
    if (system->surf_preserve) { /* preserve the orientation and only calculate based on displacement */

        //apply rotation if given in input file
        if (system->surf_preserve_rotation_on != NULL) {
            surface_dimer_geometry(system, 0.0,
                                   system->surf_preserve_rotation_on->alpha1,
                                   system->surf_preserve_rotation_on->beta1,
                                   system->surf_preserve_rotation_on->gamma1,
                                   system->surf_preserve_rotation_on->alpha2,
                                   system->surf_preserve_rotation_on->beta2,
                                   system->surf_preserve_rotation_on->gamma2, 0);
            printf(
                "SURFACE: Setting preserve angles (radians) for molecule 1: %lf %lf %lf\n",
                system->surf_preserve_rotation_on->alpha1,
                system->surf_preserve_rotation_on->beta1,
                system->surf_preserve_rotation_on->gamma1);
            printf(
                "SURFACE: Setting preserve angles (radians) for molecule 2: %lf %lf %lf\n",
                system->surf_preserve_rotation_on->alpha2,
                system->surf_preserve_rotation_on->beta2,
                system->surf_preserve_rotation_on->gamma2);
        }

        for (r = system->surf_min; r <= system->surf_max; r += system->surf_inc) {
            /* calculate the energy */
            surface_dimer_geometry(system, r, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0);

            if (system->surf_decomp) {
                pe_es = surface_energy(system, ENERGY_ES);
                pe_rd = surface_energy(system, ENERGY_RD);
                pe_polar = surface_energy(system, ENERGY_POLAR);
                pe_vdw = surface_energy(system, ENERGY_VDW);
                printf(
                    "%.5f %.5f %.5f %.5f %.5f\n", r, pe_es, pe_rd, pe_polar, pe_vdw);
                fflush(stdout);
            } else {
                pe_total = surface_energy(system, ENERGY_TOTAL);
                printf(
                    "%.5f %.5f\n", r, pe_total);
                fflush(stdout);
            }
        }

    } else { /* default is to do isotropic averaging */

        for (r = system->surf_min; r <= system->surf_max; r += system->surf_inc) {
            /* zero the averages for this r value */
            avg_counter = 0;
            pe_total_avg = 0;
            pe_es_avg = 0;
            pe_rd_avg = 0;
            pe_polar_avg = 0;
            pe_vdw_avg = 0;
            ao = 0;
            bo = 0;
            go = 0;
            am = 0;
            bm = 0;
            gm = 0;
            surf_print = 0;

            /* average over the angles */
            for (alpha_origin = 0; alpha_origin <= 2.0 * M_PI; alpha_origin += system->surf_ang) {
                for (beta_origin = 0; beta_origin <= M_PI; beta_origin += system->surf_ang) {
                    for (gamma_origin = 0; gamma_origin <= 2.0 * M_PI; gamma_origin += system->surf_ang) {
                        for (alpha_move = 0; alpha_move <= 2.0 * M_PI; alpha_move += system->surf_ang) {
                            for (beta_move = 0; beta_move <= M_PI; beta_move += system->surf_ang) {
                                for (gamma_move = 0; gamma_move <= 2.0 * M_PI; gamma_move += system->surf_ang) {
                                    ++gm;
                                    ++avg_counter;
                                    avg_factor = ((double)(avg_counter - 1)) / ((double)(avg_counter));

                                    surface_dimer_geometry(system, r, alpha_origin, beta_origin, gamma_origin, alpha_move, beta_move, gamma_move, 0);

                                    if (system->surf_print_level == 6)
                                        surf_print = 1;

                                    if (system->surf_output && surf_print == 1)
                                        write_surface_traj(system->file_pointers.fp_surf, system);

                                    surf_print = 0; /* turn off printing */

                                    if (system->surf_decomp) {
                                        pe_es = surface_energy(system, ENERGY_ES);
                                        pe_rd = surface_energy(system, ENERGY_RD);
                                        pe_polar = surface_energy(system, ENERGY_POLAR);
                                        pe_vdw = surface_energy(system, ENERGY_VDW);

                                        pe_es_avg = pe_es_avg * avg_factor + (pe_es / ((double)avg_counter));
                                        pe_rd_avg = pe_rd_avg * avg_factor + (pe_rd / ((double)avg_counter));
                                        pe_polar_avg = pe_polar_avg * avg_factor + (pe_polar / ((double)avg_counter));
                                        pe_vdw_avg = pe_vdw_avg * avg_factor + (pe_vdw / ((double)avg_counter));

                                    } else {
                                        pe_total = surface_energy(system, ENERGY_TOTAL);
                                        pe_total_avg = pe_total_avg * avg_factor + (pe_total / ((double)avg_counter));
                                    }
                                    //unrotate (not the most efficent thing to be doing, but oh well)
                                    surface_dimer_geometry(system, 0.0, alpha_origin, beta_origin, gamma_origin, alpha_move, beta_move, gamma_move, 1);

                                } /* end gamma_move */
                                ++bm;
                                if (system->surf_print_level == 5)
                                    surf_print = 1;
                            } /* end beta_move */
                            ++am;
                            if (system->surf_print_level == 4)
                                surf_print = 1;
                        } /* end alpha_move */
                        ++go;
                        if (system->surf_print_level == 3)
                            surf_print = 1;
                    } /* end gamma_origin */
                    ++bo;
                    if (system->surf_print_level == 2)
                        surf_print = 1;
                } /* end beta_origin */
                ++ao;
                if (system->surf_print_level == 1)
                    surf_print = 1;
            } /* end alpha_origin */

            if (system->surf_decomp) {
                printf(
                    "%.5f %.5f %.5f %.5f\n", r, pe_es_avg, pe_rd_avg, pe_polar_avg);
                fflush(stdout);

            } else {
                printf(
                    "%.5f %.5f\n", r, pe_total_avg);
                fflush(stdout);
            }
        } /* end r */
    }

    return (0);
}

/* set the distance and angles for a particular dimer orientation */
int surface_dimer_geometry(system_t *system, double r, double alpha_origin, double beta_origin, double gamma_origin, double alpha_move, double beta_move, double gamma_move, int reverse_flag) {
    int i;
    molecule_t *molecule_origin, *molecule_move;
    atom_t *atom_ptr;

    /* make sure that there are only two molecules */
    molecule_origin = system->molecules;
    molecule_move = molecule_origin->next;
    if (!molecule_move) {
        error(
            "SURFACE: the input PQR has only a single molecule\n");
        return (-1);
    }
    if (molecule_move->next) {
        error(
            "SURFACE: the input PQR must contain exactly two molecules\n");
        return (-1);
    }

    /* relocate both molecules to the origin */
    for (atom_ptr = molecule_origin->atoms; atom_ptr; atom_ptr = atom_ptr->next)
        for (i = 0; i < 3; i++)
            atom_ptr->pos[i] -= molecule_origin->com[i];
    for (i = 0; i < 3; i++)
        molecule_origin->com[i] = 0;

    for (atom_ptr = molecule_move->atoms; atom_ptr; atom_ptr = atom_ptr->next)
        for (i = 0; i < 3; i++)
            atom_ptr->pos[i] -= molecule_move->com[i];
    for (i = 0; i < 3; i++)
        molecule_move->com[i] = 0;

    /* relocate the moveable molecule r distance away along the x axis */
    molecule_move->com[0] = r;
    molecule_move->com[1] = 0;
    molecule_move->com[2] = 0;
    for (atom_ptr = molecule_move->atoms; atom_ptr; atom_ptr = atom_ptr->next)
        for (i = 0; i < 3; i++)
            atom_ptr->pos[i] += molecule_move->com[i];

    /* apply the rotation to the second molecule */
    if (system->surf_global_axis_on) {
        molecule_rotate_quaternion(molecule_origin, alpha_origin, beta_origin, gamma_origin, reverse_flag);
        molecule_rotate_quaternion(molecule_move, alpha_move, beta_move, gamma_move, reverse_flag);
    } else {
        molecule_rotate_euler(molecule_origin, alpha_origin, beta_origin, gamma_origin, reverse_flag);
        molecule_rotate_euler(molecule_move, alpha_move, beta_move, gamma_move, reverse_flag);
    }

    return (0);
}

// Traverse the atoms in a molecule searching for a specific id
// When the id is found, return a pointer to that atom, return NULL if
// unable to locate.
atom_t *find_atom_by_id(system_t *system, int target_id) {
    molecule_t *molecule_ptr;  // Ptr to traverse the system molecule linked list
    atom_t *atom_ptr;          // Ptr to traverse the atom list in each molecule

    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
            // Look for the bond_id, that is, the id specified in the pqr input
            // file, NOT the id stored in the id field, which is just its position
            // in the input file.
            if (atom_ptr->bond_id == target_id)
                return atom_ptr;

    return NULL;
}

// Given a vector with components (x,y,z) and a distance, dr, find_scale_factor()
// will return a scalar such that if the components of the original vector (x,y,z)
// are each multiplied by the value returned, the length of the resultant vector
// will be longer (or shorter, if dr is < 0) than the original by dr. If the
// magnitude of the original vector is 0, the function returns 0.
double find_scale_factor(double x, double y, double z, double dr) {
    // Equation was determined by setting
    // sqrt( (sx)^2 + (sy)^2 + (sz)^2 )  =  sqrt( x^2 + y^2 + z^2 ) + dr
    // and solving for s.
    double R = x * x + y * y + z * z;
    return (R ? (sqrt(1 + 2 * dr / sqrt(R) + dr * dr / R)) : 0);
}

int surface_dimer_parameters(system_t *system, param_g *params) {
    int i = 0;  // generic counter

    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    param_t *param_ptr;

    // Default bond partner for atoms without explicit site neighbor
    atom_t *origin = calloc(sizeof(atom_t), 1);  // calloc sets pos[] fields all to 0
    memnullcheck(origin, sizeof(atom_t), __LINE__ - 1, __FILE__);

    system->polar_damp = params->alpha;

    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            for (param_ptr = params->type_params; param_ptr; param_ptr = param_ptr->next) {
                if (!strcasecmp(param_ptr->atomtype, atom_ptr->atomtype)) {
                    atom_ptr->charge = param_ptr->charge;
                    atom_ptr->epsilon = param_ptr->epsilon;
                    atom_ptr->sigma = param_ptr->sigma;
                    atom_ptr->omega = param_ptr->omega;
                    atom_ptr->polarizability = param_ptr->pol;
                    atom_ptr->c6 = param_ptr->c6;
                    atom_ptr->c8 = param_ptr->c8;
                    atom_ptr->c10 = param_ptr->c10;

                    double perturbation_vector[3];

                    // Find the the atom that defines the bond axis and store for easy reference
                    atom_t *snPtr = 0;
                    if (atom_ptr->site_neighbor_id == 0)
                        snPtr = origin;
                    else {
                        snPtr = find_atom_by_id(system, atom_ptr->site_neighbor_id);
                        if (!snPtr) {
                            snPtr = origin;
                            printf(
                                "\n*** WARNING ***\n");
                            printf(
                                "Bonding partner for atom %d not found ", atom_ptr->bond_id);
                            printf(
                                "--> %d, possibly invalid. Using origin.\n\n", atom_ptr->site_neighbor_id);
                        }
                    }

                    // Generate the vector along which the dr perturbations will occur
                    for (i = 0; i < 3; i++)
                        perturbation_vector[i] = atom_ptr->pos[i] - snPtr->pos[i];

                    // Determine how much the bond vector should be scaled to increase the length by dr.
                    double scale = find_scale_factor(perturbation_vector[0], perturbation_vector[1], perturbation_vector[2], param_ptr->dr);

                    // Scale the perturbation vector
                    for (i = 0; i < 3; i++)
                        perturbation_vector[i] *= scale;

                    // Add the vector back to the bond partner coordinates to get the new position.
                    for (i = 0; i < 3; i++)
                        atom_ptr->pos[i] = snPtr->pos[i] + perturbation_vector[i];
                }
            }  // for param_ptr
        }      // for atom_ptr
    }          // for molecule_ptr

    free(origin);

    return (0);
}

/* fit potential energy parameters via simulated annealing */
int surface_fit(system_t *system) {
    //override default values if specified in input file
    double temperature = ((system->fit_start_temp) ? system->fit_start_temp : TEMPERATURE);
    double max_energy = ((system->fit_max_energy) ? system->fit_max_energy : MAX_ENERGY);
    double schedule = ((system->fit_schedule) ? system->fit_schedule : SCHEDULE);

    //used only if qshift is on
    double quadrupole = 0;

    int i, j,     // generic counters
        nCurves,  // number of curves to fit against
        nPoints,  // number of points in each curve
        nSteps;
    double *r_input = 0;     // Master list of r-values
    curveData_t *curve = 0;  // Array to hold curve point data
    param_t *param_ptr;
    param_g *params;
    double r_min, r_max, r_inc;
    double current_error, last_error, global_minimum;
    qshiftData_t *qshiftData = NULL;  //used only with qshift

    // ee_local variables
    int a = 0;
    int globalFound = 0;
    double oldGlobal = 0;
    double initEps[2];
    double initSig[2];
    double globalEps[2];
    double globalSig[2];

    // Record number of curves for convenient reference
    nCurves = system->fit_input_list.data.count;

    // Read the curve data from the input files into a local array
    curve = readFitInputFiles(system, nCurves);
    if (!curve) {
        error(
            "SURFACE: An error occurred while reading the surface-fit input files.\n");
        return (-1);
    }

    // Record number of points in each curve, for convenient reference
    nPoints = curve[0].nPoints;

    // Make a master list of R-values
    r_input = curve[0].r;
    curve[0].r = 0;

    // Free the memory used by all the remaining lists (i.e. all lists
    // other than the one whose address was just transferred).
    for (i = 1; i < nCurves; i++)
        free(curve[i].r);

    // Allocate memory for the output and global arrays
    if (-1 == alloc_curves(nCurves, nPoints, curve)) {
        free(r_input);
        return -1;
    }

    // Determine Independent Domain
    //////////////////////////////////
    r_min = r_input[0];
    r_max = r_input[nPoints - 1];
    r_inc = r_input[1] - r_input[0];

    // Output the geometry for visual verification
    output_pqrs(system, nCurves, curve);

    // Get the Initial Error
    //////////////////////////////////

    //obtain curves
    get_curves(system, nCurves, curve, r_min, r_max, r_inc);

    //calc current error
    current_error = error_calc(system, nCurves, nPoints, curve, max_energy);

    //record parameters. we will determine which ones need to be adjusted later
    params = record_params(system);

    // print some header info
    for (param_ptr = params->type_params; param_ptr; param_ptr = param_ptr->next) {
        printf(
            "SURFACE: Atom type: %d @ %s\n", param_ptr->ntypes, param_ptr->atomtype);

        if (param_ptr->epsilon != 0) {
            initEps[a] = param_ptr->epsilon;
            initSig[a] = param_ptr->sigma;
            a++;
        }
    }
    if (!system->fit_boltzmann_weight)
        printf(
            "*** any input energy values greater than %f K will not contribute to the fit ***\n", max_energy);

    //write params and current_error to stdout
    printf(
        "*** Initial Fit: \n");
    output_params(temperature, current_error, params);
    printf(
        "*****************\n");

    //set the minimum error we've found
    global_minimum = current_error;
    last_error = current_error;
    for (j = 0; j < nCurves; j++)
        for (i = 0; i < nPoints; i++)
            curve[j].global[i] = curve[j].output[i];

    //initialize stuff for qshift
    if (system->surf_qshift_on) {
        qshiftData = malloc(sizeof(qshiftData_t));
        quadrupole = calcquadrupole(system);
    }

    // ANNEALING
    // Loop over temperature. When the temperature reaches MIN_TEMPERATURE
    // we quit. Assuming we find an error smaller than the initial error
    // we'll spit out a fit_geometry.pqr file, and you'll want to re-anneal.
    for (nSteps = 0; temperature > MIN_TEMPERATURE; ++nSteps) {
        // randomly perturb the parameters
        surf_perturb(system, quadrupole, qshiftData, params);

        // apply the trial parameters
        surface_dimer_parameters(system, params);
        get_curves(system, nCurves, curve, r_min, r_max, r_inc);

        // calc error
        current_error = error_calc(system, nCurves, nPoints, curve, max_energy);

        int condition = 0;
        if (system->surf_descent)
            condition = current_error < last_error;
        else
            condition = get_rand(system) < exp(-(current_error - last_error) / temperature);

        //  DO MC at this 'temperature'
        if (condition) {
            //ACCEPT MOVE /////
            apply_new_parameters(params);
            last_error = current_error;

            if (nSteps >= EQUIL) {
                nSteps = 0;

                //write params and current_error to stdout
                output_params(temperature, current_error, params);

                // output the new global minimum
                if (current_error < global_minimum) {
                    system->fit_best_square_error = global_minimum = current_error;
                    new_global_min(system, nCurves, nPoints, curve);

                    // Store the LJ parameters corresponding to global min
                    a = 0;
                    globalFound = 1;
                    oldGlobal = current_error;  // keep track for ee_local
                    for (param_ptr = params->type_params; param_ptr; param_ptr = param_ptr->next) {
                        if (param_ptr->epsilon != 0) {
                            globalEps[a] = param_ptr->epsilon;
                            globalSig[a] = param_ptr->sigma;
                            a++;
                        }
                    }
                }
            }
        } else {
            // MOVE REJECTED ////
            revert_parameters(system, params);
        }  //END DO MONTE CARLO

        // decrement the temperature
        temperature = temperature * schedule;  // schedule
    }                                          //END ANNEALING

    // Output the Fit Curves
    output_fit(nCurves, nPoints, curve, max_energy, r_input);

    // Exhaustive Enumeration -- Chris Cioce
    if (system->ee_local) {
        int b, c, d, nE1, nE2, nS1, nS2;
        double currEps[2];
        double currSig[2];
        double newGlobal = 0;
        double E1start, E1stop, E1incr, E2start, E2stop, E2incr, E1c, E2c;  // E1, E2 (unused)
        double S1start, S1stop, S1incr, S2start, S2stop, S2incr, S1c, S2c;  // S1, S2 (unused)

        if (globalFound == 1) {
            printf(
                "\nEE_LOCAL ~> SA found an improved surface. Using current global minimum parameters for EE.\n");
            currEps[0] = globalEps[0];
            currEps[1] = globalEps[1];
            currSig[0] = globalSig[0];
            currSig[1] = globalSig[1];
        } else {
            printf(
                "\nEE_LOCAL ~> SA did not find an improved surface. Using initial parameters for EE.\n");
            currEps[0] = initEps[0];
            currEps[1] = initEps[1];
            currSig[0] = initSig[0];
            currSig[1] = initSig[1];
        }

        for (a = 0; a < 2; a++) {
            printf(
                "EE_LOCAL ~> InitEps[%d]: %f\tInitSig[%d]: %f\tGlobalEps[%d]: %f\tGlobalSig[%d]: %f\n", a, initEps[a], a, initSig[a], a, globalEps[a], a, globalSig[a]);
        }

        E1start = currEps[0] - (currEps[0] * system->range_eps);
        E1stop = currEps[0] + (currEps[0] * system->range_eps);
        E1incr = currEps[0] * system->step_eps;
        E2start = currEps[1] - (currEps[1] * system->range_eps);
        E2stop = currEps[1] + (currEps[1] * system->range_eps);
        E2incr = currEps[1] * system->step_eps;
        S1start = currSig[0] - (currSig[0] * system->range_sig);
        S1stop = currSig[0] + (currSig[0] * system->range_sig);
        S1incr = currSig[0] * system->step_sig;
        S2start = currSig[1] - (currSig[1] * system->range_sig);
        S2stop = currSig[1] + (currSig[1] * system->range_sig);
        S2incr = currSig[1] * system->step_sig;
        nE1 = ((E1stop - E1start) / E1incr) + 1;
        nE2 = ((E2stop - E2start) / E2incr) + 1;
        nS1 = ((S1stop - S1start) / S1incr) + 1;
        nS2 = ((S2stop - S2start) / S2incr) + 1;
        printf(
            "EE_LOCAL ~> E1start = %f\tE1stop = %f\tE1incr = %f\tnE1 = %d\n", E1start, E1stop, E1incr, nE1);
        printf(
            "EE_LOCAL ~> E2start = %f\tE2stop = %f\tE2incr = %f\tnE2 = %d\n", E2start, E2stop, E2incr, nE2);
        printf(
            "EE_LOCAL ~> S1start = %f\tS1stop = %f\tS1incr = %f\tnS1 = %d\n", S1start, S1stop, S1incr, nS1);
        printf(
            "EE_LOCAL ~> S2start = %f\tS2stop = %f\tS2incr = %f\tnS2 = %d\n", S2start, S2stop, S2incr, nS2);
        printf(
            "EE_LOCAL ~> Will now perform %d total iterations in parameter subspace...\n", (nE1 * nE2 * nS1 * nS2));

        E1c = E1start;
        E2c = E2start;
        S1c = S1start;
        S2c = S2start;

        // Set initial parameters
        for (param_ptr = params->type_params; param_ptr; param_ptr = param_ptr->next) {
            if (param_ptr->epsilon == currEps[0]) {
                param_ptr->epsilon = E1c;
            }
            if (param_ptr->epsilon == currEps[1]) {
                param_ptr->epsilon = E2c;
            }
            if (param_ptr->sigma == currSig[0]) {
                param_ptr->sigma = S1c;
            }
            if (param_ptr->sigma == currSig[1]) {
                param_ptr->sigma = S2c;
            }
        }

        // apply initial parameters
        surface_dimer_parameters(system, params);
        get_curves(system, nCurves, curve, r_min, r_max, r_inc);

        //calc current error
        current_error = error_calc(system, nCurves, nPoints, curve, max_energy);

        //record parameters. we will determine which ones need to be adjusted later
        //params = record_params ( system );

        //write params and current_error to stdout <-- just for visual verification
        printf(
            "*** Initial Fit: \n");
        output_params(0.0, current_error, params);
        printf(
            "*****************\n");

        // Main loops
        globalFound = 0;
        int master = 1;
        for (a = 0; a <= nE1; a++) {
            for (b = 0; b <= nE2; b++) {
                for (c = 0; c <= nS1; c++) {
                    for (d = 0; d <= nS2; d++) {
                        //printf("MASTER: %d\ta: %d\t(%f)\tb: %d\t(%f)\tc: %d\t(%f)\td: %d\t(%f)\n", master, a, E1c, b, E2c, c, S1c, d, S2c);	// Debug output

                        // apply parameters
                        surface_dimer_parameters(system, params);
                        get_curves(system, nCurves, curve, r_min, r_max, r_inc);

                        // calc current error
                        current_error = error_calc(system, nCurves, nPoints, curve, max_energy);

                        //output_params ( 0.0, current_error, params );	// Debug output

                        // check for new global minimum
                        if (current_error < global_minimum) {
                            system->fit_best_square_error = global_minimum = current_error;
                            new_global_min(system, nCurves, nPoints, curve);

                            // Store the LJ parameters corresponding to global min
                            a = 0;
                            globalFound = 1;
                            newGlobal = current_error;
                            for (param_ptr = params->type_params; param_ptr; param_ptr = param_ptr->next) {
                                if (param_ptr->epsilon != 0) {
                                    globalEps[a] = param_ptr->epsilon;
                                    globalSig[a] = param_ptr->sigma;
                                    a++;
                                }
                            }
                        }

                        // Set updated S2 parameter
                        for (param_ptr = params->type_params; param_ptr; param_ptr = param_ptr->next) {
                            if (param_ptr->sigma == S2c) {
                                param_ptr->sigma = (S2c + S2incr);
                            }
                        }
                        S2c += S2incr;

                        master++;
                    }  // END S2 loop

                    // Set updated S1 & S2 parameters
                    for (param_ptr = params->type_params; param_ptr; param_ptr = param_ptr->next) {
                        if (param_ptr->sigma == S1c) {
                            param_ptr->sigma = (S1c + S1incr);
                        }
                        if (param_ptr->sigma == S2c) {
                            param_ptr->sigma = S2start;
                        }
                    }

                    S1c += S1incr;
                    S2c = S2start;

                }  // END S1 loop

                // Set updated E2 & S1 parameters
                for (param_ptr = params->type_params; param_ptr; param_ptr = param_ptr->next) {
                    if (param_ptr->epsilon == E2c) {
                        param_ptr->epsilon = (E2c + E2incr);
                    }
                    if (param_ptr->sigma == S1c) {
                        param_ptr->sigma = S1start;
                    }
                }

                E2c += E2incr;
                S1c = S1start;

            }  // END E2 loop

            // Set updated E1 & E2 parameters
            for (param_ptr = params->type_params; param_ptr; param_ptr = param_ptr->next) {
                if (param_ptr->epsilon == E1c) {
                    param_ptr->epsilon = (E1c + E1incr);
                }
                if (param_ptr->epsilon == E2c) {
                    param_ptr->epsilon = E2start;
                }
            }

            E1c += E1incr;
            E2c = E2start;

        }  // END E1 loop

        // If a new global min is found, write to stdout
        // TODO: Adjust printf statements to prettify final output :D
        if (globalFound == 1) {
            printf(
                "\nEE_LOCAL ~> EE found an improved surface. Parameters written to fit_geometry.pqr\n");
            //printf("EE_LOCAL ~> E1: %f\tS1: %f\tE2: %f\tS2: %f\n", globalEps[0], globalSig[0], globalEps[1], globalSig[1]);
            printf(
                "EE_LOCAL ~> SA Global Min: %f\n", oldGlobal);
            printf(
                "EE_LOCAL ~> EE Global Min: %f\n", newGlobal);
            printf(
                "            ----------------------------\n");
            printf(
                "            ----------------------------\n");
            printf(
                "               ====> %.2f %% improvement\n\n", ((oldGlobal - newGlobal) / oldGlobal) * 100);
        } else {
            printf(
                "EE_LOCAL ~> EE was unable to find an improved surface versus SA. Oh well, at least you tried...and for that, we thank you!\n\n");
        }

    }  // END ee_local

    // Return all memory back to the system
    free_all_mem(nCurves, curve, params, qshiftData, r_input);
    return (0);
}
