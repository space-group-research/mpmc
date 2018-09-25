#include <mc.h>

int read_pqr_box(FILE *fp, system_t *system) {
    char buffer[MAXLINE], token[7][MAXLINE];
    char msg[MAXLINE];
    int basis_set[3];

    output(
        "INPUT: (read_pqr_box) checking input pqr for basis info\n");

    //flags to make sure we set all basis vectors
    basis_set[0] = basis_set[1] = basis_set[2] = 0;
    while (fgets(buffer, MAXLINE, fp) != NULL) {
        sscanf(buffer,
               "%s %s %s %s %s %s %s",
               token[0], token[1], token[2], token[3], token[4], token[5], token[6]);

        if ((!strncmp(token[0],
                      "END", 3))) break;  //if end of molecule, then stop searching

        if ((!strcmp(token[0],
                     "REMARK")) &&
            (!strcmp(token[1],
                     "BOX")) &&
            (!strcmp(token[3],
                     "="))) {
            if (!strcmp(token[2],
                        "BASIS[0]")) {  //set basis[0]
                {
                    if (safe_atof(token[4], &(system->pbc->basis[0][0]))) continue;
                }  //make sure each conversion is successful
                {
                    if (safe_atof(token[5], &(system->pbc->basis[0][1]))) continue;
                }
                {
                    if (safe_atof(token[6], &(system->pbc->basis[0][2]))) continue;
                }
                //if we get this far, then we've successfully read in the basis vector
                basis_set[0] = 1;
            }
            if (!strcmp(token[2],
                        "BASIS[1]")) {  //set basis[0]
                {
                    if (safe_atof(token[4], &(system->pbc->basis[1][0]))) continue;
                }  //make sure each conversion is successful
                {
                    if (safe_atof(token[5], &(system->pbc->basis[1][1]))) continue;
                }
                {
                    if (safe_atof(token[6], &(system->pbc->basis[1][2]))) continue;
                }
                //if we get this far, then we've successfully read in the basis vector
                basis_set[1] = 1;
            }
            if (!strcmp(token[2],
                        "BASIS[2]")) {  //set basis[0]
                {
                    if (safe_atof(token[4], &(system->pbc->basis[2][0]))) continue;
                }  //make sure each conversion is successful
                {
                    if (safe_atof(token[5], &(system->pbc->basis[2][1]))) continue;
                }
                {
                    if (safe_atof(token[6], &(system->pbc->basis[2][2]))) continue;
                }
                //if we get this far, then we've successfully read in the basis vector
                basis_set[2] = 1;
            } else
                continue;
        } else
            continue;
    }

    if (basis_set[0] == 1) {
        sprintf(msg,
                "INPUT: basis[0] successfully read from pqr {%.5lf %.5lf %.5lf}\n",
                system->pbc->basis[0][0], system->pbc->basis[0][1], system->pbc->basis[0][2]);
        output(msg);
    } else {
        sprintf(msg,
                "INPUT: unable to read basis[0] from pqr file.\n");
        error(msg);
    }
    if (basis_set[1] == 1) {
        sprintf(msg,
                "INPUT: basis[1] successfully read from pqr {%.5lf %.5lf %.5lf}\n",
                system->pbc->basis[1][0], system->pbc->basis[1][1], system->pbc->basis[1][2]);
        output(msg);
    } else {
        sprintf(msg,
                "INPUT: unable to read basis[1] from pqr file.\n");
        error(msg);
    }
    if (basis_set[2] == 1) {
        sprintf(msg,
                "INPUT: basis[2] successfully read from pqr {%.5lf %.5lf %.5lf}\n",
                system->pbc->basis[2][0], system->pbc->basis[2][1], system->pbc->basis[2][2]);
        output(msg);
    } else {
        sprintf(msg,
                "INPUT: unable to read basis[2] from pqr file.\n");
        error(msg);
    }

    return 0;
}

#ifdef QM_ROTATION
void allocqmrotation(system_t *system, molecule_t *molecule_ptr) {
    int i;

    molecule_ptr->quantum_rotational_energies = calloc(system->quantum_rotation_level_max, sizeof(double));
    memnullcheck(molecule_ptr->quantum_rotational_energies,
                 system->quantum_rotation_level_max * sizeof(double), __LINE__ - 1, __FILE__);
    molecule_ptr->quantum_rotational_eigenvectors = calloc(system->quantum_rotation_level_max, sizeof(complex_t *));
    memnullcheck(molecule_ptr->quantum_rotational_eigenvectors,
                 system->quantum_rotation_level_max * sizeof(complex_t *), __LINE__ - 1, __FILE__);
    for (i = 0; i < system->quantum_rotation_level_max; i++) {
        molecule_ptr->quantum_rotational_eigenvectors[i] =
            calloc((system->quantum_rotation_l_max + 1) * (system->quantum_rotation_l_max + 1), sizeof(complex_t));
        memnullcheck(molecule_ptr->quantum_rotational_eigenvectors[i],
                     (system->quantum_rotation_l_max + 1) * (system->quantum_rotation_l_max + 1) * sizeof(complex_t), __LINE__ - 1, __FILE__);
    }
    molecule_ptr->quantum_rotational_eigensymmetry = calloc(system->quantum_rotation_level_max, sizeof(int));
    memnullcheck(molecule_ptr->quantum_rotational_eigensymmetry,
                 system->quantum_rotation_level_max * sizeof(int), __LINE__ - 1, __FILE__);

    return;
}
#endif

/* unused
void allocqmvibration(system_t * system, molecule_t * molecule_ptr) {

	molecule_ptr->quantum_vibrational_energies = calloc(system->quantum_vibration_level_max, sizeof(double));
	memnullcheck(molecule_ptr->quantum_vibrational_energies,
		system->quantum_vibration_level_max*sizeof(double),__LINE__-1, __FILE__);
	molecule_ptr->quantum_vibrational_eigenvectors = calloc(system->quantum_vibration_level_max, sizeof(complex_t *));
	memnullcheck(molecule_ptr->quantum_vibrational_eigenvectors,	
		system->quantum_vibration_level_max*sizeof(complex_t *),__LINE__-1, __FILE__);
	for(i = 0; i < system->quantum_vibration_level_max; i++) {
		molecule_ptr->quantum_vibrational_eigenvectors[i] = 
			calloc((system->quantum_vibration_l_max+1)*(system->quantum_vibration_l_max+1), sizeof(complex_t));
		memnullcheck(molecule_ptr->quantum_vibrational_eigenvectors[i],
			(system->quantum_vibration_l_max+1)*(system->quantum_vibration_l_max+1)*sizeof(complex_t),__LINE__-1, __FILE__);
	}
	molecule_ptr->quantum_vibrational_eigensymmetry = calloc(system->quantum_vibration_level_max, sizeof(int));
	memnullcheck(molecule_ptr->quantum_vibrational_eigensymmetry,
		system->quantum_vibration_level_max*sizeof(int),__LINE__-1, __FILE__);

	return;
}
*/

molecule_t *read_molecules(FILE *fp, system_t *system) {
    molecule_t *molecules, *molecule_ptr;
    atom_t *atom_ptr, *prev_atom_ptr;
    char linebuf[MAXLINE];
    char token_atom[MAXLINE], token_atomid[MAXLINE], token_atomtype[MAXLINE], token_moleculetype[MAXLINE],
        token_frozen[MAXLINE], token_moleculeid[MAXLINE], token_x[MAXLINE], token_y[MAXLINE], token_z[MAXLINE],
        token_mass[MAXLINE], token_charge[MAXLINE], token_alpha[MAXLINE], token_epsilon[MAXLINE],
        token_sigma[MAXLINE], token_omega[MAXLINE], token_gwp_alpha[MAXLINE], token_c6[MAXLINE], token_c8[MAXLINE], token_c10[MAXLINE], token_c9[MAXLINE];
    int current_frozen, current_adiabatic, current_spectre, current_target, current_moleculeid,
        current_atomid, current_site_neighbor, moveable, spectres, targets, atom_counter;
    double current_x, current_y, current_z,
        current_mass, current_charge, current_alpha, current_epsilon,
        current_sigma, current_omega, current_gwp_alpha, current_c6,
        current_c8, current_c10, current_c9;  //, current_molecule_mass; (unused variable)

    fpos_t file_pos;
    fgetpos(fp, &file_pos);  //get file pointer position, we will restore this when done

    /* allocate the start of the list */
    molecules = calloc(1, sizeof(molecule_t));
    memnullcheck(molecules, sizeof(molecule_t), __LINE__ - 1, __FILE__);
    molecule_ptr = molecules;
    molecule_ptr->id = 1;
    molecule_ptr->atoms = calloc(1, sizeof(atom_t));
    memnullcheck(molecule_ptr->atoms, sizeof(atom_t), __LINE__ - 1, __FILE__);
    atom_ptr = molecule_ptr->atoms;
    prev_atom_ptr = atom_ptr;

    /* clear the linebuffer and read the tokens in */
    atom_counter = 0;
    memset(linebuf, 0, MAXLINE);

    while (fgets(linebuf, MAXLINE, fp)) {
        /* clear the tokens */
        memset(token_atom, 0, MAXLINE);
        memset(token_atomid, 0, MAXLINE);
        memset(token_atomtype, 0, MAXLINE);
        memset(token_moleculetype, 0, MAXLINE);
        memset(token_frozen, 0, MAXLINE);
        memset(token_moleculeid, 0, MAXLINE);
        memset(token_x, 0, MAXLINE);
        memset(token_y, 0, MAXLINE);
        memset(token_z, 0, MAXLINE);
        memset(token_mass, 0, MAXLINE);
        memset(token_charge, 0, MAXLINE);
        memset(token_alpha, 0, MAXLINE);
        memset(token_epsilon, 0, MAXLINE);
        memset(token_sigma, 0, MAXLINE);
        memset(token_omega, 0, MAXLINE);
        memset(token_gwp_alpha, 0, MAXLINE);
        memset(token_c6, 0, MAXLINE);
        memset(token_c8, 0, MAXLINE);
        memset(token_c10, 0, MAXLINE);
        memset(token_c9, 0, MAXLINE);

        /* parse the line */
        sscanf(linebuf,
               "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
               token_atom, token_atomid, token_atomtype, token_moleculetype, token_frozen,
               token_moleculeid, token_x, token_y, token_z, token_mass, token_charge,
               token_alpha, token_epsilon, token_sigma, token_omega, token_gwp_alpha,
               token_c6, token_c8, token_c10, token_c9);

        if (!strncasecmp(token_atom,
                         "END", 3)) break;  //we've reached the end of the current molecule, quit looping

        if (!strcasecmp(token_atom,
                        "ATOM") &&
            strcasecmp(token_moleculetype,
                       "BOX")) {
            current_frozen = 0;
            current_adiabatic = 0;
            current_spectre = 0;
            current_target = 0;
            if (!strcasecmp(token_frozen,
                            "F"))
                current_frozen = 1;
            if (!strcasecmp(token_frozen,
                            "A"))
                current_adiabatic = 1;
            if (!strcasecmp(token_frozen,
                            "S"))
                current_spectre = 1;
            if (!strcasecmp(token_frozen,
                            "T"))
                current_target = 1;

            current_moleculeid = atoi(token_moleculeid);
            current_atomid = atoi(token_atomid);
            current_x = atof(token_x);
            current_y = atof(token_y);
            current_z = atof(token_z);
            current_mass = atof(token_mass); /* mass in amu */
            current_charge = atof(token_charge);
            current_charge *= E2REDUCED; /* convert charge into reduced units  */
            current_alpha = atof(token_alpha);
            current_epsilon = atof(token_epsilon);
            current_sigma = atof(token_sigma);
            current_omega = atof(token_omega);
            current_gwp_alpha = atof(token_gwp_alpha);
            current_c6 = atof(token_c6);
            current_c8 = atof(token_c8);
            current_c10 = atof(token_c10);
            current_c9 = atof(token_c9);

            if (system->cdvdw_sig_repulsion) {
                if (current_epsilon != 1.0) {
                    error(
                        "warning: setting epsilon to 1.0 (due to sig_repulsion)\n");
                    current_epsilon = 1.0;
                }
            } else if (system->polarvdw && !system->cdvdw_exp_repulsion) {
                if (current_sigma != 1.0) {
                    error(
                        "warning: setting sigma to 1.0 (due to polarvdw)\n");
                    current_sigma = 1.0;
                }
            }
            // Functionality of site_neighbor disabled in favor of omega/gwp_alpha parameters
            // Current behavior is to default to atom 0, typically the center of mass for
            // the molecule.
            current_site_neighbor = 0;  //atoi( token_site_neighbor );

            if (current_frozen)
                current_charge *= system->scale_charge;

            if (molecule_ptr->id != current_moleculeid) {
                molecule_ptr->next = calloc(1, sizeof(molecule_t));
                memnullcheck(molecule_ptr, sizeof(molecule_t), __LINE__ - 1, __FILE__);
                molecule_ptr = molecule_ptr->next;
                molecule_ptr->atoms = calloc(1, sizeof(atom_t));
                memnullcheck(molecule_ptr->atoms, sizeof(atom_t), __LINE__ - 1, __FILE__);
                prev_atom_ptr->next = NULL;
                free(atom_ptr);
                atom_ptr = molecule_ptr->atoms;
            }

            strcpy(molecule_ptr->moleculetype, token_moleculetype);

            molecule_ptr->id = current_moleculeid;
            molecule_ptr->frozen = current_frozen;
            molecule_ptr->adiabatic = current_adiabatic;
            molecule_ptr->spectre = current_spectre;
            molecule_ptr->target = current_target;
            molecule_ptr->mass += current_mass;

#ifdef QM_ROTATION
            /* if quantum rot calc. enabled, allocate the necessary structures */
            if (system->quantum_rotation && !molecule_ptr->frozen && !molecule_ptr->quantum_rotational_energies)
                allocqmrotation(system, molecule_ptr);
#endif /* QM_ROTATION */

#ifdef XXX
            /* if quantum vib calc. enabled, allocate the necessary structures */
            if (system->quantum_vibration && !molecule_ptr->frozen)
                allocqmvibration(system, molecule_ptr);
#endif /* XXX */

            ++atom_counter;
            atom_ptr->id = atom_counter;
            atom_ptr->bond_id = current_atomid;
            memset(atom_ptr->atomtype, 0, MAXLINE);
            strcpy(atom_ptr->atomtype, token_atomtype);
            atom_ptr->frozen = current_frozen;
            atom_ptr->adiabatic = current_adiabatic;
            atom_ptr->spectre = current_spectre;
            atom_ptr->target = current_target;
            atom_ptr->pos[0] = current_x;
            atom_ptr->pos[1] = current_y;
            atom_ptr->pos[2] = current_z;
            atom_ptr->mass = current_mass;
            atom_ptr->charge = current_charge;
            atom_ptr->polarizability = current_alpha;
            atom_ptr->epsilon = current_epsilon;
            atom_ptr->sigma = current_sigma;
            atom_ptr->omega = current_omega;
            atom_ptr->gwp_alpha = current_gwp_alpha;
            atom_ptr->c6 = current_c6;
            atom_ptr->c8 = current_c8;
            atom_ptr->c10 = current_c10;
            atom_ptr->c9 = current_c9;
            if (current_gwp_alpha != 0.)
                atom_ptr->gwp_spin = 1;
            else
                atom_ptr->gwp_spin = 0;

            atom_ptr->site_neighbor_id = current_site_neighbor;
            atom_ptr->next = calloc(1, sizeof(atom_t));
            memnullcheck(atom_ptr->next, sizeof(atom_t), __LINE__ - 1, __FILE__);
            prev_atom_ptr = atom_ptr;
            atom_ptr = atom_ptr->next;
        }

        memset(linebuf, 0, MAXLINE);
    }

    /* terminate the atom list */
    prev_atom_ptr->next = NULL;
    free(atom_ptr);

    /* scan the list, make sure that there is at least one moveable molecule */
    moveable = 0;
    spectres = 0;
    targets = 0;
    for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        if (!molecule_ptr->frozen) ++moveable;
        if (molecule_ptr->target) ++targets;
        if (molecule_ptr->spectre) ++spectres;
    }

    if (system->spectre) {
        if (!spectres || !targets) {
            error(
                "INPUT: either no targets or spectres found\n");
            return (NULL);
        }
    } else {
        if (!moveable) {
            error(
                "INPUT: no moveable molecules found, there must be at least one in your PQR file\n");
            return (NULL);
        }
    }
    if (!atom_counter) {
        free(molecules);
        free(molecule_ptr);
        return (NULL);
    }

    //if we're going to read box info, we need to rewind our file
    if (system->read_pqr_box_on)
        fsetpos(fp, &file_pos);

    return (molecules);
}

//  read_insertion_molecules( system_t * ) was a cut/paste of read_molecules
//  modified to read the candidate insertion molecules from a separate file.
///////////////////////////////////////////////////////////////////////////////

// helper functions which will test for presence of a sorbate in
// the sorbate list and which will add the sorbate if necessary.
int sorbateIs_Not_InList(system_t *, char *);
void addSorbateToList(system_t *, char *);

molecule_t *read_insertion_molecules(system_t *system) {
    int j;

    molecule_t *molecules,
        *molecule_ptr;

    atom_t *atom_ptr,
        *prev_atom_ptr;

    char linebuf[MAXLINE], *n;

    FILE *fp;
    char token_atom[MAXLINE], token_atomid[MAXLINE], token_atomtype[MAXLINE],
        token_moleculetype[MAXLINE],
        token_frozen[MAXLINE],
        token_moleculeid[MAXLINE],
        token_x[MAXLINE], token_y[MAXLINE], token_z[MAXLINE],
        token_mass[MAXLINE],
        token_charge[MAXLINE],
        token_alpha[MAXLINE], token_epsilon[MAXLINE], token_sigma[MAXLINE],
        token_omega[MAXLINE], token_gwp_alpha[MAXLINE],
        token_c6[MAXLINE], token_c8[MAXLINE], token_c10[MAXLINE], token_c9[MAXLINE];

    int current_frozen,
        current_adiabatic,
        current_spectre,
        current_target,
        current_moleculeid,
        current_atomid,
        current_site_neighbor;
    double current_x, current_y, current_z,
        current_mass, current_charge,
        current_alpha, current_epsilon, current_sigma, current_omega, current_gwp_alpha,
        current_c6, current_c8, current_c10, current_c9;

    int atom_counter;

    // allocate the start of the list
    molecules = calloc(1, sizeof(molecule_t));
    memnullcheck(molecules, sizeof(molecule_t), __LINE__ - 1, __FILE__);
    molecule_ptr = molecules;
    molecule_ptr->id = 1;
    molecule_ptr->atoms = calloc(1, sizeof(atom_t));
    memnullcheck(molecule_ptr->atoms, sizeof(atom_t), __LINE__ - 1, __FILE__);
    atom_ptr = molecule_ptr->atoms;
    prev_atom_ptr = atom_ptr;

    // open the molecule input file
    fp = fopen(system->insert_input,
               "r");
    filecheck(fp, system->insert_input, READ);

    // clear the linebuffer and read the tokens in
    atom_counter = 0;
    memset(linebuf, 0, MAXLINE);
    n = fgets(linebuf, MAXLINE, fp);
    while (n) {
        // clear the tokens
        memset(token_atom, 0, MAXLINE);
        memset(token_atomid, 0, MAXLINE);
        memset(token_atomtype, 0, MAXLINE);
        memset(token_moleculetype, 0, MAXLINE);
        memset(token_frozen, 0, MAXLINE);
        memset(token_moleculeid, 0, MAXLINE);
        memset(token_x, 0, MAXLINE);
        memset(token_y, 0, MAXLINE);
        memset(token_z, 0, MAXLINE);
        memset(token_mass, 0, MAXLINE);
        memset(token_charge, 0, MAXLINE);
        memset(token_alpha, 0, MAXLINE);
        memset(token_epsilon, 0, MAXLINE);
        memset(token_sigma, 0, MAXLINE);
        memset(token_omega, 0, MAXLINE);
        memset(token_gwp_alpha, 0, MAXLINE);
        memset(token_c6, 0, MAXLINE);
        memset(token_c8, 0, MAXLINE);
        memset(token_c10, 0, MAXLINE);
        memset(token_c9, 0, MAXLINE);

        // parse the input line
        sscanf(linebuf,
               "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
               token_atom, token_atomid, token_atomtype,
               token_moleculetype,
               token_frozen,
               token_moleculeid,
               token_x, token_y, token_z,
               token_mass, token_charge,
               token_alpha, token_epsilon, token_sigma, token_omega, token_gwp_alpha,
               token_c6, token_c8, token_c10, token_c9);

        if (!strcasecmp(token_atom,
                        "ATOM") &&
            strcasecmp(token_moleculetype,
                       "BOX")) {
            current_frozen = 0;
            current_adiabatic = 0;
            current_spectre = 0;
            current_target = 0;
            if (!strcasecmp(token_frozen,
                            "F"))
                current_frozen = 1;
            if (!strcasecmp(token_frozen,
                            "A"))
                current_adiabatic = 1;
            if (!strcasecmp(token_frozen,
                            "S"))
                current_spectre = 1;
            if (!strcasecmp(token_frozen,
                            "T"))
                current_target = 1;

            current_moleculeid = atoi(token_moleculeid);
            current_atomid = atoi(token_atomid);
            current_x = atof(token_x);
            current_y = atof(token_y);
            current_z = atof(token_z);
            current_mass = atof(token_mass);  // mass in amu
            current_charge = atof(token_charge);
            current_charge *= E2REDUCED;  // convert charge into reduced units
            current_alpha = atof(token_alpha);
            current_epsilon = atof(token_epsilon);
            current_sigma = atof(token_sigma);
            current_omega = atof(token_omega);
            current_gwp_alpha = atof(token_gwp_alpha);
            current_c6 = atof(token_c6);
            current_c8 = atof(token_c8);
            current_c10 = atof(token_c10);
            current_c9 = atof(token_c9);

            if (system->cdvdw_sig_repulsion) {
                if (current_epsilon != 1.0) {
                    error(
                        "warning: setting epsilon to 1.0 (due to sig_repulsion)\n");
                    current_epsilon = 1.0;
                }
            } else if (system->polarvdw) {
                if (current_sigma != 1.0) {
                    error(
                        "warning: setting sigma to 1.0 (due to polarvdw)\n");
                    current_sigma = 1.0;
                }
            }

            // Functionality of site_neighbor disabled in favor of omega/gwp_alpha parameters
            // Current behavior is to default to atom 0, typically the center of mass for
            // the molecule.
            current_site_neighbor = 0;

            if (current_frozen)
                current_charge *= system->scale_charge;

            if (molecule_ptr->id != current_moleculeid) {
                molecule_ptr->next = calloc(1, sizeof(molecule_t));
                memnullcheck(molecule_ptr->next, sizeof(molecule_t), __LINE__ - 1, __FILE__);
                molecule_ptr = molecule_ptr->next;
                molecule_ptr->atoms = calloc(1, sizeof(atom_t));
                memnullcheck(molecule_ptr->atoms, sizeof(atom_t), __LINE__ - 1, __FILE__);
                prev_atom_ptr->next = NULL;
                free(atom_ptr);
                atom_ptr = molecule_ptr->atoms;
            }
            strcpy(molecule_ptr->moleculetype, token_moleculetype);

            molecule_ptr->id = current_moleculeid;
            molecule_ptr->frozen = current_frozen;
            molecule_ptr->adiabatic = current_adiabatic;
            molecule_ptr->spectre = current_spectre;
            molecule_ptr->target = current_target;
            molecule_ptr->mass += current_mass;

#ifdef QM_ROTATION
            /* if quantum rot calc. enabled, allocate the necessary structures */
            if (system->quantum_rotation && !molecule_ptr->frozen && !molecule_ptr->quantum_rotational_energies)
                allocqmrotation(system, molecule_ptr);
#endif /* QM_ROTATION */

#ifdef XXX
            /* if quantum vib calc. enabled, allocate the necessary structures */
            if (system->quantum_vibration && !molecule_ptr->frozen)
                allocqmvibration(system, molecule_ptr);
#endif /* XXX */

            ++atom_counter;
            atom_ptr->id = atom_counter;
            atom_ptr->bond_id = current_atomid;
            memset(atom_ptr->atomtype, 0, MAXLINE);
            strcpy(atom_ptr->atomtype, token_atomtype);
            atom_ptr->frozen = current_frozen;
            atom_ptr->adiabatic = current_adiabatic;
            atom_ptr->spectre = current_spectre;
            atom_ptr->target = current_target;
            atom_ptr->pos[0] = current_x;
            atom_ptr->pos[1] = current_y;
            atom_ptr->pos[2] = current_z;
            atom_ptr->mass = current_mass;
            atom_ptr->charge = current_charge;
            atom_ptr->polarizability = current_alpha;
            atom_ptr->epsilon = current_epsilon;
            atom_ptr->sigma = current_sigma;
            atom_ptr->omega = current_omega;
            atom_ptr->gwp_alpha = current_gwp_alpha;
            atom_ptr->c6 = current_c6;
            atom_ptr->c8 = current_c8;
            atom_ptr->c10 = current_c10;
            atom_ptr->c9 = current_c9;
            if (current_gwp_alpha != 0.)
                atom_ptr->gwp_spin = 1;
            else
                atom_ptr->gwp_spin = 0;

            atom_ptr->site_neighbor_id = current_site_neighbor;
            atom_ptr->next = calloc(1, sizeof(atom_t));
            memnullcheck(atom_ptr->next, sizeof(atom_t), __LINE__ - 1, __FILE__);
            prev_atom_ptr = atom_ptr;
            atom_ptr = atom_ptr->next;
        }

        memset(linebuf, 0, MAXLINE);
        n = fgets(linebuf, MAXLINE, fp);
    }

    // terminate the atom list
    prev_atom_ptr->next = NULL;
    free(atom_ptr);

    // Count the molecules and create an array of pointers, where each pointer
    // points directly to a molecule in the linked list. While counting the
    // molecules, check to see whether or not each sorbate is included in the
    // sorbate statistics list, if it is not, we will create an entry for it
    // so that each unique class of sorbate will have its own node in the stats
    // list.
    /////////////////////////////////////////////////////////////////////////////

    // count
    int molecule_counter = 0;
    for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        molecule_counter++;
        // Check to see if this molecule is a new sorbate.
        if (sorbateIs_Not_InList(system, molecule_ptr->moleculetype)) {
            system->sorbateCount++;
            addSorbateToList(system, molecule_ptr->moleculetype);
            for (j = 0; j < system->sorbateCount; j++) {
                if (!strcasecmp(system->sorbateInfo[j].id, molecule_ptr->moleculetype)) {
                    system->sorbateInfo[j].mass = molecule_ptr->mass;
                    break;
                }
            }
        }
    }

    // allocate space for array that will give direct access to insertion-candidate
    // molecules in the linked list.
    system->insertion_molecules_array = malloc(molecule_counter * sizeof(molecule_t *));
    memnullcheck(system->insertion_molecules_array, molecule_counter * sizeof(molecule_t *), __LINE__ - 1, __FILE__);

    // point array pointers to their corresponding molecules
    molecule_counter = 0;
    for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        system->insertion_molecules_array[molecule_counter] = molecule_ptr;
        molecule_counter++;
    }
    system->num_insertion_molecules = molecule_counter;

    //  ERROR CHECKING OF SOME SORT MAY NEED TO GO HERE
    //   WHAT'S ALLOWED/DISALLOWED FOR INSERTED MOLECULES?
    /////////////////////////////////////////////////////////

    fclose(fp);

    return (molecules);
}

// sorbateIs_Not_InList() examines a sorbate stat list and a sorbate ID.
// Returns true if this sorbate is NOT in the list, false otherwise
int sorbateIs_Not_InList(system_t *system, char *type) {
    int i;

    // Look at each item in the list until we find a match
    // with the sorbate id about which the inquiry was made.
    for (i = 0; i < system->sorbateCount; i++) {
        if (!strcasecmp(type, system->sorbateInfo[i].id))
            // sorbate is already accounted for
            return 0;
    }

    // If end of list is reached w/o encountering the sorbate, then
    // the sorbate is new... return 1, i.e. the statement "sorbate
    // is NOT in list" is TRUE.
    return 1;
}

void addSorbateToList(system_t *system, char *sorbate_type) {
    // grow array
    system->sorbateInfo = realloc(system->sorbateInfo, system->sorbateCount * sizeof(sorbateInfo_t));
    memnullcheck(system->sorbateInfo, system->sorbateCount * sizeof(sorbateInfo_t), __LINE__ - 1, __FILE__);

    // set sorbate type for the new element
    strcpy(system->sorbateInfo[system->sorbateCount - 1].id, sorbate_type);

    // send a status update to stdout
    char buffer[MAXLINE];
    sprintf(buffer,
            "INPUT: System will track individual stats for sorbate %s.\n",
            system->sorbateInfo[system->sorbateCount - 1].id);
    output(buffer);

    return;
}

#ifdef DEBUG
void test_list(molecule_t *molecules) {
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr;

    for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        printf(
            "DEBUG_LIST: moleculeid = %d\n", molecule_ptr->id);
        printf(
            "DEBUG_LIST: moleculetype = %s\n", molecule_ptr->moleculetype);
        printf(
            "DEBUG_LIST: molecule_frozen = %d\n", molecule_ptr->frozen);
        printf(
            "DEBUG_LIST: molecule_mass = %f\n", molecule_ptr->mass);
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            printf(
                "DEBUG_LIST: atomtype = %s x = %f y = %f z = %f\n", atom_ptr->atomtype, atom_ptr->pos[0], atom_ptr->pos[1], atom_ptr->pos[2]);
            printf(
                "DEBUG_LIST: atom frozen = %d mass = %f, charge = %f, alpha = %f, eps = %f, sig = %f\n", atom_ptr->frozen, atom_ptr->mass, atom_ptr->charge, atom_ptr->polarizability, atom_ptr->epsilon, atom_ptr->sigma);
            for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next)
                if (!(pair_ptr->rd_excluded || pair_ptr->es_excluded || pair_ptr->frozen)) printf(
                    "DEBUG_LIST: pair = 0x%lx eps = %f sig = %f\n",
                    (long unsigned int)pair_ptr, pair_ptr->epsilon, pair_ptr->sigma);
        }
    }

    fflush(stdout);
}

void test_molecule(molecule_t *molecule) {
    atom_t *atom_ptr;
    pair_t *pair_ptr;

    printf(
        "moleculeid = %d\n", molecule->id);
    printf(
        "moleculetype = %s\n", molecule->moleculetype);
    printf(
        "molecule_frozen = %d\n", molecule->frozen);
    printf(
        "molecule_mass = %f\n", molecule->mass);
    for (atom_ptr = molecule->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
        printf(
            "atomtype = %s x = %f y = %f z = %f\n", atom_ptr->atomtype, atom_ptr->pos[0], atom_ptr->pos[1], atom_ptr->pos[2]);
        printf(
            "atom frozen = %d mass = %f, charge = %f, alpha = %f, eps = %f, sig = %f\n", atom_ptr->frozen, atom_ptr->mass, atom_ptr->charge, atom_ptr->polarizability, atom_ptr->epsilon, atom_ptr->sigma);
        for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
            printf(
                "pair at 0x%lx\n", (long unsigned int)pair_ptr);
            fflush(stdout);
        }
    }

    printf(
        "...finished\n");
    fflush(stdout);
}
#endif /* DEBUG */
