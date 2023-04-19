#include <stdbool.h>
#include <surface_multi_fit.h>

// Copyright Adam Hogan 2016-2021 - GNU GPL v3

void load_initial_multi_params(system_t *system, multiParamData_t *params) {
    // The PQR input file should have one atom per atom type with the initial parameters
    int nAtoms = 0, i = 0, j;

    // Count the number of atoms
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
            nAtoms++;

    // Allocate param arrays
    params->atomtype = calloc(nAtoms, sizeof(char *));
    memnullcheck(params->atomtype, nAtoms * sizeof(char *), __LINE__ - 1, __FILE__);
    params->c6 = calloc(nAtoms, sizeof(double));
    memnullcheck(params->c6, nAtoms * sizeof(double), __LINE__ - 1, __FILE__);
    params->last_c6 = calloc(nAtoms, sizeof(double));
    memnullcheck(params->last_c6, nAtoms * sizeof(double), __LINE__ - 1, __FILE__);
    params->c8 = calloc(nAtoms, sizeof(double));
    memnullcheck(params->c8, nAtoms * sizeof(double), __LINE__ - 1, __FILE__);
    params->last_c8 = calloc(nAtoms, sizeof(double));
    memnullcheck(params->last_c8, nAtoms * sizeof(double), __LINE__ - 1, __FILE__);
    params->c10 = calloc(nAtoms, sizeof(double));
    memnullcheck(params->c10, nAtoms * sizeof(double), __LINE__ - 1, __FILE__);
    params->last_c10 = calloc(nAtoms, sizeof(double));
    memnullcheck(params->last_c10, nAtoms * sizeof(double), __LINE__ - 1, __FILE__);
    params->epsilon = calloc(nAtoms, sizeof(double));
    memnullcheck(params->epsilon, nAtoms * sizeof(double), __LINE__ - 1, __FILE__);
    params->last_epsilon = calloc(nAtoms, sizeof(double));
    memnullcheck(params->last_epsilon, nAtoms * sizeof(double), __LINE__ - 1, __FILE__);
    params->sigma = calloc(nAtoms, sizeof(double));
    memnullcheck(params->sigma, nAtoms * sizeof(double), __LINE__ - 1, __FILE__);
    params->last_sigma = calloc(nAtoms, sizeof(double));
    memnullcheck(params->last_sigma, nAtoms * sizeof(double), __LINE__ - 1, __FILE__);
    params->omega = calloc(nAtoms, sizeof(double));
    memnullcheck(params->omega, nAtoms * sizeof(double), __LINE__ - 1, __FILE__);
    params->last_omega = calloc(nAtoms, sizeof(double));
    memnullcheck(params->last_omega, nAtoms * sizeof(double), __LINE__ - 1, __FILE__);
    params->polarizability = calloc(nAtoms, sizeof(double));
    memnullcheck(params->polarizability, nAtoms * sizeof(double), __LINE__ - 1, __FILE__);
    params->last_polarizability = calloc(nAtoms, sizeof(double));
    memnullcheck(params->last_polarizability, nAtoms * sizeof(double), __LINE__ - 1, __FILE__);

    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            if (i >= nAtoms) {
                error("SURFACE MULTI FIT: Error loading initial parameters %s:%s\n");
                die(1);
            }
            params->atomtype[i] = calloc(MAXLINE, sizeof(char));
            memnullcheck(params->atomtype[i], MAXLINE * sizeof(char), __LINE__ - 1, __FILE__);
            strcpy(params->atomtype[i], atom_ptr->atomtype);
            params->c6[i] = atom_ptr->c6;
            params->last_c6[i] = params->c6[i];
            params->c8[i] = atom_ptr->c8;
            params->last_c8[i] = params->c8[i];
            params->c10[i] = atom_ptr->c10;
            params->last_c10[i] = params->c10[i];
            params->epsilon[i] = atom_ptr->epsilon;
            params->last_epsilon[i] = params->epsilon[i];
            params->sigma[i] = atom_ptr->sigma;
            params->last_sigma[i] = params->sigma[i];
            params->omega[i] = atom_ptr->omega;
            params->last_omega[i] = params->omega[i];
            params->polarizability[i] = atom_ptr->polarizability;
            params->last_polarizability[i] = params->polarizability[i];
            i++;
            params->nParams = i;
        }
    }

    for (i = 0; i < params->nParams; i++) {
        for (j = 0; j < params->nParams; j++) {
            if (!(i == j)) {
                if (strcasecmp(params->atomtype[i], params->atomtype[j]) == 0) {
                    printf("SURFACE MULTI FIT: Atomic parameter (%s) (in input pqr) duplicated.\n", params->atomtype[i]);
                    die(-1);
                }
            }
        }
    }

    return;
}

void apply_multi_params(multiParamData_t *params, multiConfigData_t *configs) {
    int i, j;
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    for (i = 0; i < configs->nConfigs; i++) {
        for (molecule_ptr = configs->molecules[i]; molecule_ptr; molecule_ptr = molecule_ptr->next) {
            for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
                for (j = 0; j < params->nParams; j++) {
                    if (strcasecmp(atom_ptr->atomtype, params->atomtype[j]) == 0) {
                        atom_ptr->c6 = params->c6[j];
                        atom_ptr->c8 = params->c8[j];
                        atom_ptr->c10 = params->c10[j];
                        atom_ptr->epsilon = params->epsilon[j];
                        atom_ptr->sigma = params->sigma[j];
                        atom_ptr->omega = params->omega[j];
                        atom_ptr->polarizability = params->polarizability[j];
                        break;
                    }
                }
            }
        }
    }
    return;
}

void accept_multi_params(multiParamData_t *params, multiConfigData_t *configs) {
    int i;
    for (i = 0; i < params->nParams; i++) {
        params->last_c6[i] = params->c6[i];
        params->last_c8[i] = params->c8[i];
        params->last_c10[i] = params->c10[i];
        params->last_epsilon[i] = params->epsilon[i];
        params->last_sigma[i] = params->sigma[i];
        params->last_omega[i] = params->omega[i];
        params->last_polarizability[i] = params->polarizability[i];
    }
    for (i = 0; i < configs->nConfigs; i++) {
        configs->lastFitEnergy[i] = configs->fitEnergy[i];
    }
    return;
}

void undo_multi_params(multiParamData_t *params, multiConfigData_t *configs) {
    int i;
    for (i = 0; i < params->nParams; i++) {
        params->c6[i] = params->last_c6[i];
        params->c8[i] = params->last_c8[i];
        params->c10[i] = params->last_c10[i];
        params->epsilon[i] = params->last_epsilon[i];
        params->sigma[i] = params->last_sigma[i];
        params->omega[i] = params->last_omega[i];
        params->polarizability[i] = params->last_polarizability[i];
    }
    for (i = 0; i < configs->nConfigs; i++) {
        configs->fitEnergy[i] = configs->lastFitEnergy[i];
    }
    return;
}

void read_multi_configs(system_t *system, multiConfigData_t *configs, multiParamData_t *params) {
    /*
    config input format, any number of following with newlines separating:
    line 0: comment technically, preferably labeled "Configuration #"
    line 1: single double with no whitespace, the ab initio energy
    line 2-N+2: any number of atoms, which consist of a string, int, double, double, double, double
        where the first string is the atomtype, first int is the molecule #, followed by the x, y and z positions and finally the charge
    */
    char *filename = system->multi_fit_input;
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        char error_message[MAXLINE];
        sprintf(error_message, "SURFACE MULTI FIT: Couldn't open multi-fit configuration file %s.\n", filename);
        error(error_message);
        die(-1);
    }
    char line[MAXLINE], *token, *return_value;
    const char s[] = " \t";

    return_value = fgets(line, sizeof(line), fp);
    if (!return_value) {
        error("SURFACE MULTI FIT: The multi-fit configuration file does not contain configurations.\n");
        die(-1);
    }
    printf("%s\n", s);
    token = strtok(line, s);
    if (strcmp("Configuration", token)) {
        error("SURFACE MULTI FIT: The multi-fit configuration file does not contain configurations.\n");
        die(-1);
    }

    // Count the number of configs in the file
    configs->nConfigs = 1;
    while (fgets(line, sizeof(line), fp)) {
        token = strtok(line, s);
        if (strcmp("Configuration", token) == 0) {
            configs->nConfigs += 1;
        }
    }

    configs->abInitioEnergy = calloc(configs->nConfigs, sizeof(double));
    memnullcheck(configs->abInitioEnergy, configs->nConfigs * sizeof(double), __LINE__ - 1, __FILE__);
    configs->fitEnergy = calloc(configs->nConfigs, sizeof(double));
    memnullcheck(configs->fitEnergy, configs->nConfigs * sizeof(double), __LINE__ - 1, __FILE__);
    configs->lastFitEnergy = calloc(configs->nConfigs, sizeof(double));
    memnullcheck(configs->lastFitEnergy, configs->nConfigs * sizeof(double), __LINE__ - 1, __FILE__);
    configs->molecules = calloc(configs->nConfigs, sizeof(molecule_t *));
    memnullcheck(configs->molecules, configs->nConfigs * sizeof(molecule_t *), __LINE__ - 1, __FILE__);
    configs->periodic = calloc(configs->nConfigs, sizeof(int));
    memnullcheck(configs->periodic, configs->nConfigs * sizeof(int), __LINE__ - 1, __FILE__);
    configs->pbc = calloc(configs->nConfigs, sizeof(pbc_t *));
    memnullcheck(configs->pbc, configs->nConfigs * sizeof(pbc_t *), __LINE__ - 1, __FILE__);
    int i, j, k;
    for (i=0;i<configs->nConfigs;i++)
    {
        configs->pbc[i] = calloc(1, sizeof(pbc_t));
        memnullcheck(configs->pbc[i], 1 * sizeof(pbc_t), __LINE__ - 1, __FILE__);
    }

    i = -1;
    rewind(fp);
    while (fgets(line, sizeof(line), fp)) {
        token = strtok(line, s);
        if (strcmp("Configuration", token) == 0) {
            // Setup new molecule linked list here
            i++;
            configs->molecules[i] = calloc(1, sizeof(molecule_t));
            memnullcheck(configs->molecules[i], sizeof(molecule_t), __LINE__ - 1, __FILE__);
            configs->molecules[i]->atoms = NULL;
            configs->molecules[i]->next = calloc(1, sizeof(molecule_t));
            memnullcheck(configs->molecules[i]->next, sizeof(molecule_t), __LINE__ - 1, __FILE__);
            configs->molecules[i]->next->atoms = NULL;
            fgets(line, sizeof(line), fp);
            token = strtok(line, s);
            if (safe_atof(token,&(configs->abInitioEnergy[i])))
            {
                error("SURFACE MULTI FIT: Unable to read ab initio energy in multi-fit configuration file.\n");
                die(-1);
            }
            configs->periodic[i] = false;
            system->pbc = configs->pbc[i];
            system->pbc->basis[0][0] = 100000.;
            system->pbc->basis[0][1] = 0.;
            system->pbc->basis[0][2] = 0.;
            system->pbc->basis[1][0] = 0.;
            system->pbc->basis[1][1] = 100000.;
            system->pbc->basis[1][2] = 0.;
            system->pbc->basis[2][0] = 0.;
            system->pbc->basis[2][1] = 0.;
            system->pbc->basis[2][2] = 100000.;
            pbc(system);
            token = strtok(NULL, s);
            if (token != NULL)
            {
                if (strcmp("periodic", token) == 0) {
                    configs->periodic[i] = true;
                    double a, b, c, alpha, beta, gamma;
                    if (safe_atof(strtok(NULL, s), &a))
                    {
                        error("SURFACE MULTI FIT: Improper periodic specification in multi-fit configuration file.\n");
                        die(-1);
                    }
                    if (safe_atof(strtok(NULL, s), &b))
                    {
                        error("SURFACE MULTI FIT: Improper periodic specification in multi-fit configuration file.\n");
                        die(-1);
                    }
                    if (safe_atof(strtok(NULL, s), &c))
                    {
                        error("SURFACE MULTI FIT: Improper periodic specification in multi-fit configuration file.\n");
                        die(-1);
                    }
                    if (safe_atof(strtok(NULL, s), &alpha))
                    {
                        error("SURFACE MULTI FIT: Improper periodic specification in multi-fit configuration file.\n");
                        die(-1);
                    }
                    if (safe_atof(strtok(NULL, s), &beta))
                    {
                        error("SURFACE MULTI FIT: Improper periodic specification in multi-fit configuration file.\n");
                        die(-1);
                    }
                    if (safe_atof(strtok(NULL, s), &gamma))
                    {
                        error("SURFACE MULTI FIT: Improper periodic specification in multi-fit configuration file.\n");
                        die(-1);
                    }
                    car2basis(system, a, b, c, alpha, beta, gamma);
                    pbc(system);
                }
            }
        } else {
            atom_t *new_atom = calloc(1, sizeof(atom_t));
            memnullcheck(new_atom, sizeof(atom_t), __LINE__ - 1, __FILE__);
            new_atom->next = NULL;
            strcpy(new_atom->atomtype, token);

            // Check that each atom has a corresponding atom in the parameter file
            bool name_exists = false;
            for (j = 0; j < params->nParams; j++)
            {
                if (strcasecmp(params->atomtype[j], new_atom->atomtype) == 0)
                {
                    name_exists = true;
                }
            }
            if (!name_exists)
            {
                printf("SURFACE MULTI FIT: Atom in multi-fit configuration file (%s) is not specified in the atomic parameter file (input pqr).\n", new_atom->atomtype);
                die(-1);
            }

            token = strtok(NULL, s);
            if (token == NULL)
            {
                error("SURFACE MULTI FIT: The multi-fit configuration file is not in the correct format.\n");
                die(-1);
            }
            int molecule_number = atoi(token);
            token = strtok(NULL, s);
            if (token == NULL)
            {
                error("SURFACE MULTI FIT: The multi-fit configuration file is not in the correct format.\n");
                die(-1);
            }
            safe_atof(token,&(new_atom->pos[0]));
            new_atom->wrapped_pos[0] = new_atom->pos[0];
            token = strtok(NULL, s);
            if (token == NULL)
            {
                error("SURFACE MULTI FIT: The multi-fit configuration file is not in the correct format.\n");
                die(-1);
            }
            safe_atof(token,&(new_atom->pos[1]));
            new_atom->wrapped_pos[1] = new_atom->pos[1];
            token = strtok(NULL, s);
            if (token == NULL)
            {
                error("SURFACE MULTI FIT: The multi-fit configuration file is not in the correct format.\n");
                die(-1);
            }
            safe_atof(token,&(new_atom->pos[2]));
            new_atom->wrapped_pos[2] = new_atom->pos[2];
            token = strtok(NULL, s);
            if (token == NULL)
            {
                printf("SURFACE MULTI FIT: Charge not specified in multi-fit configuration file, defaulting to 0 charge.\n");
                //die(-1);
                new_atom->charge = 0;  // actually just set the charge to zero instead of dieing
            }
            else
            {
                safe_atof(token,&(new_atom->charge));
                new_atom->charge *= E2REDUCED;
            }
            molecule_t *insert_in_this_molecule = configs->molecules[i];

            for (j=0; j<molecule_number-1; j++)
            {
                if (insert_in_this_molecule->next == NULL) {
                    insert_in_this_molecule->next = calloc(1, sizeof(molecule_t));
                    memnullcheck(insert_in_this_molecule->next, sizeof(molecule_t), __LINE__ - 1, __FILE__);
                    insert_in_this_molecule->next->atoms = NULL;
                }
                insert_in_this_molecule = insert_in_this_molecule->next;
            }

            if (insert_in_this_molecule->atoms == NULL)
                insert_in_this_molecule->atoms = new_atom;
            else {
                atom_t *last_atom = insert_in_this_molecule->atoms;
                while (last_atom->next != NULL) {
                    last_atom = last_atom->next;
                }
                last_atom->next = new_atom;
            }
        }
    }

    // Check that each atomic parameter has atleast one corresponding atom in the multi-fit config file
    for (k = 0; k < params->nParams; k++) {
        bool name_exists = false;
        for (j = 0; j < configs->nConfigs; j++) {
            molecule_t *molecule_ptr;
            atom_t *atom_ptr;
            for (molecule_ptr = configs->molecules[j]; molecule_ptr; molecule_ptr = molecule_ptr->next) {
                for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
                    if (strcasecmp(params->atomtype[k], atom_ptr->atomtype) == 0) {
                        name_exists = true;
                    }
                }
            }
        }
        if (!name_exists) {
            printf("SURFACE MULTI FIT: Atomic parameter (%s) (in input pqr) not found in multi-fit configuration file.\n", params->atomtype[k]);
            die(-1);
        }
    }

    // Setup pairs
    apply_multi_params(params, configs);
    // Keep track of the highest number of atoms in a config for the polarization matrices
    int highest_n = 0;
    for (i = 0; i < configs->nConfigs; i++) {
        configs->fitEnergy[i] = MAXVALUE;
        system->molecules = configs->molecules[i];
        molecule_t *molecule_ptr;
        atom_t *atom_ptr;
        pair_t *pair_ptr;
        int n = 0;

        for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
            for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
                n++;
                for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next)
                    pair_ptr->recalculate_energy = 1;
            }

        system->natoms = n;
        if (n > highest_n)
            highest_n = n;

        setup_pairs(system);
    }

    highest_n *= 3;

    // Allocate the A and B polarization matrices
    system->A_matrix = calloc(highest_n, sizeof(double *));
    memnullcheck(system->A_matrix, highest_n * sizeof(double *), __LINE__ - 1, __FILE__);

    for (i = 0; i < highest_n; i++) {
        system->A_matrix[i] = malloc(highest_n * sizeof(double));
        memnullcheck(system->A_matrix[i], highest_n * sizeof(double), __LINE__ - 1, __FILE__);
    }

    if (!system->polar_iterative) {
        system->B_matrix = calloc(highest_n, sizeof(double *));
        memnullcheck(system->B_matrix, highest_n * sizeof(double *), __LINE__ - 1, __FILE__);
        for (i = 0; i < highest_n; i++) {
            system->B_matrix[i] = malloc(highest_n * sizeof(double));
            memnullcheck(system->B_matrix[i], highest_n * sizeof(double), __LINE__ - 1, __FILE__);
        }
    }

    return;
}

void perturb_multi_params(system_t *system, multiParamData_t *params) {
    int i, j;

    //override default values if specified in input file
    double scale_epsilon = ((system->surf_scale_epsilon_on) ? system->surf_scale_epsilon : SCALE_EPSILON);
    double scale_sigma = ((system->surf_scale_sigma_on) ? system->surf_scale_sigma : SCALE_SIGMA);
    double scale_omega = ((system->surf_scale_omega_on) ? system->surf_scale_omega : SCALE_OMEGA);
    double scale_c6 = ((system->surf_scale_c6_on) ? system->surf_scale_c6 : 0.0);
    double scale_c8 = ((system->surf_scale_c8_on) ? system->surf_scale_c8 : 0.0);
    double scale_c10 = ((system->surf_scale_c10_on) ? system->surf_scale_c10 : 0.0);
    double scale_pol = ((system->surf_scale_pol_on) ? system->surf_scale_pol : 0.0);

    int do_not_fit_length;

    if (system->surf_do_not_fit_list == NULL) {
        do_not_fit_length = 0;
    } else {
        for (i = 0; strlen(system->surf_do_not_fit_list[i]) > 0; i++)
            ;
        do_not_fit_length = i;
    }

    for (i = 0; i < params->nParams; i++) {
        int fit = 1;
        for (j = 0; j < do_not_fit_length; j++) {
            if (strcasecmp(params->atomtype[i], system->surf_do_not_fit_list[j]) == 0)
                fit = 0;
        }

        if (fit) {
            if (params->epsilon[i] > 0.0) {
                if (system->polarvdw)
                    params->epsilon[i] += params->epsilon[i] * (scale_epsilon * (0.5 - get_rand(system)));
                else
                    params->epsilon[i] += scale_epsilon * (0.5 - get_rand(system));
                if (params->epsilon[i] < 0.0) params->epsilon[i] = params->last_epsilon[i];
            }

            //if polarvdw is on, also adjust omega
            if (system->polarvdw || system->disp_expansion_mbvdw) {
                if (params->omega[i] > 0.0)
                    params->omega[i] += scale_omega * (0.5 - get_rand(system));
                if (params->omega[i] < 0.0) params->omega[i] = params->last_omega[i];
                if (system->cdvdw_exp_repulsion) {
                    if (params->sigma[i] > 0.0)  //scale_sigma given as percentage when exp_rep is on
                        params->sigma[i] += params->sigma[i] * (scale_sigma * (0.5 - get_rand(system)));
                    if (params->sigma[i] < 0.0) params->sigma[i] = params->last_sigma[i];
                }
            }
            //otherwise adjust sigma (adjust sigma in the case of disp_expansion_mbvdw as well)
            else if (!system->polarvdw) {
                if (params->sigma[i] > 0.0)
                    params->sigma[i] += scale_sigma * (0.5 - get_rand(system));
                if (params->sigma[i] < 0.0) params->sigma[i] = params->last_sigma[i];
            }

            if (system->surf_scale_c6_on) {
                if (params->c6[i] > 0.0)
                    params->c6[i] += scale_c6 * (0.5 - get_rand(system));
                if (params->c6[i] < 0.0) params->c6[i] = params->last_c6[i];
            }

            if (system->surf_scale_c8_on) {
                if (params->c8[i] > 0.0)
                    params->c8[i] += scale_c8 * (0.5 - get_rand(system));
                if (params->c8[i] < 0.0) params->c8[i] = params->last_c8[i];
            }

            if (system->surf_scale_c10_on) {
                if (params->c10[i] > 0.0)
                    params->c10[i] += scale_c10 * (0.5 - get_rand(system));
                if (params->c10[i] < 0.0) params->c10[i] = params->last_c10[i];
            }

            if (system->surf_scale_pol_on) {
                if (params->polarizability[i] > 0.0)
                    params->polarizability[i] += scale_pol * (0.5 - get_rand(system));
                if (params->polarizability[i] < 0.0) params->polarizability[i] = params->last_polarizability[i];
            }
        }
    }

    return;
}

double calc_multi_configurational_energy(system_t *system) {
    double potential_energy, rd_energy, coulombic_energy, polar_energy, vdw_energy, three_body_energy;

    /* zero the initial values */
    potential_energy = 0;
    rd_energy = 0;
    coulombic_energy = 0;
    polar_energy = 0;
    vdw_energy = 0;
    three_body_energy = 0;

    /* get the pairwise terms necessary for the energy calculation */
    pairs(system);

    flag_all_pairs(system);

    coulombic_energy = coulombic_nopbc(system->molecules);

    /* get the polarization potential */
    if (system->polarization) {
#ifdef CUDA
        if (system->cuda)
        {
            polar_cuda(system);
            polar_energy = system->observables->polarization_energy;
        }
        else
            polar_energy = polar(system);
#else
        polar_energy = polar(system);
#endif /* CUDA */
    }

    /* get the repulsion/dispersion potential */
    if (system->sg)
        rd_energy = sg(system);
    else if (system->dreiding)
        rd_energy = dreiding(system);
    else if (system->lj_buffered_14_7)
        rd_energy = lj_buffered_14_7(system);
    else if (system->disp_expansion)
        rd_energy = disp_expansion(system);
    else if (system->rd_anharmonic)
        rd_energy = anharmonic(system);
    else if (system->cdvdw_exp_repulsion)
        rd_energy = exp_repulsion(system);
    else if (!system->gwp)
        rd_energy = lj(system);

    if (system->polarvdw) {
#ifdef CUDA
        if (system->cuda) {
            error("error: cuda polarvdw not yet implemented!\n");
            die(-1);
        } else
            vdw_energy = vdw(system);
#else
        vdw_energy = vdw(system);
#endif
    }

    if (system->axilrod_teller)
        three_body_energy = axilrod_teller(system);

    /* sum the total potential energy */
    potential_energy = rd_energy + coulombic_energy + polar_energy + vdw_energy + three_body_energy;
    return (potential_energy);
}

double calc_multi_error(system_t *system, multiConfigData_t *configs) {
    double total_error = 0.0;
    int i;
    double max_energy = system->fit_max_energy;
    double kweight = ((system->surf_weight_constant_on) ? system->surf_weight_constant : WEIGHT_CONSTANT);

    for (i = 0; i < configs->nConfigs; i++) {
        system->molecules = configs->molecules[i];
        system->pbc = configs->pbc[i];
        double model_energy = calc_multi_configurational_energy(system);
        if (configs->fitEnergy[i] == MAXVALUE)
            configs->lastFitEnergy[i] = model_energy;
        configs->fitEnergy[i] = model_energy;
        model_energy = min(model_energy, max_energy);
        double ab_initio_energy = configs->abInitioEnergy[i];
        ab_initio_energy = min(ab_initio_energy, max_energy);
        double weight = exp(kweight * (max_energy - ab_initio_energy) / max_energy);
        double error = model_energy - ab_initio_energy;
        total_error += weight * error * error;
    }

    return total_error / configs->nConfigs;
}

double calc_multi_mse(system_t *system, multiConfigData_t *configs) {
    double mse = 0.0;
    int i, count = 0;
    double max_energy = system->fit_max_energy;

    for (i = 0; i < configs->nConfigs; i++) {
        system->molecules = configs->molecules[i];
        system->pbc = configs->pbc[i];
        double current_energy = calc_multi_configurational_energy(system);
        double abInitio = configs->abInitioEnergy[i];

        if (abInitio <= max_energy) {
            mse += current_energy - abInitio;
            count += 1;
        }
    }
    return mse / count;
}

double calc_multi_mue(system_t *system, multiConfigData_t *configs) {
    double mue = 0.0;
    int i, count = 0;
    double max_energy = system->fit_max_energy;

    for (i = 0; i < configs->nConfigs; i++) {
        system->molecules = configs->molecules[i];
        system->pbc = configs->pbc[i];
        double current_energy = calc_multi_configurational_energy(system);
        double abInitio = configs->abInitioEnergy[i];

        if (abInitio <= max_energy) {
            mue += fabs(current_energy - abInitio);
            count += 1;
        }
    }
    return mue / count;
}

void output_multi_params(double temperature, double current_error, multiParamData_t *params, system_t *system, molecule_t *original_molecule, multiConfigData_t *configs) {
    printf("temperature = %f, error_func = %f\n", temperature, current_error);
    //printf("Mean Unsigned Error = %f, Mean Signed Error = %f\n", calc_multi_mue(system,configs), calc_multi_mse(system,configs));
    int i, j;
    for (i = 0; i < params->nParams; i++) {
        printf("\tatomtype %s eps = %f sig = %f pol = %f omega = %f c6 = %f c8 = %f c10 = %f\n",
               params->atomtype[i],
               params->epsilon[i],
               params->sigma[i],
               params->polarizability[i],
               params->omega[i],
               params->c6[i],
               params->c8[i],
               params->c10[i]);
    }
    // write the params file
    system->wrapall = 0;
    FILE *fp = fopen("multifit_sites.pqr", "w");
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    for (molecule_ptr = original_molecule; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            for (j = 0; j < params->nParams; j++) {
                if (strcasecmp(atom_ptr->atomtype, params->atomtype[j]) == 0) {
                    atom_ptr->c6 = params->c6[j];
                    atom_ptr->c8 = params->c8[j];
                    atom_ptr->c10 = params->c10[j];
                    atom_ptr->epsilon = params->epsilon[j];
                    atom_ptr->sigma = params->sigma[j];
                    atom_ptr->omega = params->omega[j];
                    atom_ptr->polarizability = params->polarizability[j];
                    break;
                }
            }
        }
    }
    system->molecules = original_molecule;
    system->fit_best_square_error = current_error;
    write_molecules(system, fp);
    fclose(fp);
    // write the current configs
    fp = fopen("multifit_traj.pqr", "w");
    for (i = 0; i < configs->nConfigs; i++) {
        system->molecules = configs->molecules[i];
        write_molecules(system, fp);
    }
    fclose(fp);
    fflush(stdout);
    return;
}

void output_best_config_energies(multiConfigData_t *configs) {
    int i;
    printf("\t               Ab Initio            Fit\n");
    for (i = 0; i < configs->nConfigs; i++) {
        printf("\t%8d %20.12f %20.12f\n", i, configs->abInitioEnergy[i], configs->fitEnergy[i]);
    }
    return;
}

int surface_multi_fit(system_t *system) {
    multiConfigData_t *configs = calloc(1, sizeof(multiConfigData_t));
    memnullcheck(configs, sizeof(multiConfigData_t), __LINE__ - 1, __FILE__);
    multiParamData_t *params = calloc(1, sizeof(multiParamData_t));
    memnullcheck(params, sizeof(multiParamData_t), __LINE__ - 1, __FILE__);

    molecule_t *original_molecule = system->molecules;
    pbc_t *original_pbc = system->pbc;
    load_initial_multi_params(system, params);
    read_multi_configs(system, configs, params);

    double temperature = ((system->fit_start_temp) ? system->fit_start_temp : TEMPERATURE);  // The formatting here is so pretty, thanks Brant
    double max_energy = ((system->fit_max_energy) ? system->fit_max_energy : MAX_ENERGY);
    double schedule = ((system->fit_schedule) ? system->fit_schedule : SCHEDULE);
    system->fit_max_energy = max_energy;

    int i, nSteps;

    // print some header info
    for (i = 0; i < params->nParams; i++)
        printf("SURFACE: Atom type: %d @ %s\n", i + 1, params->atomtype[i]);
    if (!system->fit_boltzmann_weight)
        printf("*** any input energy values greater than %f K will not contribute to the fit ***\n", max_energy);

    double current_error = calc_multi_error(system, configs);
    double last_error = current_error;

    // write params and current_error to stdout
    printf("*** Initial Fit: \n");
    output_multi_params(temperature, current_error, params, system, original_molecule, configs);
    printf("*****************\n");

    for (nSteps = 0; temperature > MIN_TEMPERATURE; ++nSteps) {
        perturb_multi_params(system, params);
        apply_multi_params(params, configs);
        current_error = calc_multi_error(system, configs);

        int condition = 0;
        if (system->surf_descent)
            condition = current_error < last_error;
        else
            condition = get_rand(system) < exp(-(current_error - last_error) / temperature);

        //  DO MC at this 'temperature'
        if (condition) {
            // ACCEPT
            // TODO: update last parameters in params
            output_multi_params(temperature, current_error, params, system, original_molecule, configs);
            last_error = current_error;
            accept_multi_params(params, configs);
        } else {
            // REJECT
            //output_multi_params(temperature, current_error, params, system, original_molecule, configs);
            undo_multi_params(params, configs);
            apply_multi_params(params, configs);
        }

        temperature = temperature * schedule;
    }

    // write params and current_error to stdout
    printf("*** Final Fit: \n");
    current_error = calc_multi_error(system, configs);
    output_multi_params(temperature, current_error, params, system, original_molecule, configs);
    output_best_config_energies(configs);
    printf("*****************\n");

    free(params->atomtype);
    free(params->c6);
    free(params->last_c6);
    free(params->c8);
    free(params->last_c8);
    free(params->c10);
    free(params->last_c10);
    free(params->epsilon);
    free(params->last_epsilon);
    free(params->sigma);
    free(params->last_sigma);
    free(params->omega);
    free(params->last_omega);
    free(params->polarizability);
    free(params->last_polarizability);
    free(params);


    for (i=0;i<configs->nConfigs;i++)
    {
        free(configs->pbc[i]);
        free_all_molecules(configs->molecules[i]);
    }

    free(configs->abInitioEnergy);
    free(configs->fitEnergy);
    free(configs->lastFitEnergy);
    free(configs->molecules);
    free(configs->periodic);
    free(configs->pbc);
    free(configs);

    system->molecules = original_molecule;
    system->pbc = original_pbc;
    pbc(system);
    pairs(system);

    return 0;
}
