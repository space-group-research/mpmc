#include <mc.h>
#include <surface_fit_via_arbitrary_configs.h>

// calculates the current error for a surface energy versus the ab-initio curve
double error_calc_fr_arbitrary_configurations(system_t *system, int nConfigs, configData_t *configuration, double max_energy) {
    int i;
    double abInitio, model, weight;
    double current_error = 0;

    // adjust weight constant in the weighting exponential function for calculating sq_error
    double kweight = ((system->surf_weight_constant_on) ? system->surf_weight_constant : WEIGHT_CONSTANT);

    for (i = 0; i < nConfigs; i++) {
        if (system->fit_boltzmann_weight == 0) {
            //changed so all points contribute until both curves (ab-init and fit) exceed max_energy
            abInitio = min(max_energy, configuration[i].abInitioEnergy);
            model = min(max_energy, configuration[i].currentFitEnergy);
            weight = exp(kweight * (max_energy - abInitio) / max_energy);
        } else {
            abInitio = configuration[i].abInitioEnergy;
            model = configuration[i].currentFitEnergy;
            weight = exp(kweight * (-abInitio) / max_energy);
        }

#ifdef DEBUG
        printf(
            "abInitio: %0.10lf    model: %0.10lf    weight: %0.10lf\n", abInitio, model, weight);
#endif

        current_error += weight * (abInitio - model) * (abInitio - model);
    }

    return current_error;
}

/*

void output_xyz( system_t *sys ) {
	molecule_t *molPtr;
	atom_t     *atmPtr;

	for( molPtr = sys->molecules; molPtr; molPtr = molPtr->next )
		for( atmPtr = molPtr->atoms; atmPtr; atmPtr = atmPtr->next ) {
			char out[1024];
			sprintf( out,  "H %12.5lf %12.5lf %12.5lf\n", atmPtr->pos[0], atmPtr->pos[1], atmPtr->pos[2] );
			output( out );
		}
}
*/

void dumpconfigs(system_t *system) {
    molecule_t *moleculePtr;
    atom_t *atomPtr;
    int i;

    printf(
        "SYSTEM CONFIGURATION DUMP:\n");

    for (i = 0, moleculePtr = system->molecules; moleculePtr; moleculePtr = moleculePtr->next, i++) {
        printf(
            "\tMOLECULE %d\n", i);
        for (atomPtr = moleculePtr->atoms; atomPtr; atomPtr = atomPtr->next) {
            printf(
                "\t\t%8s: %0.3lf  %0.3lf  %0.3lf   sig %0.5lf omega %0.5lf\n", atomPtr->atomtype, atomPtr->pos[0], atomPtr->pos[1], atomPtr->pos[2], atomPtr->sigma, atomPtr->omega);
        }
    }

    return;
}

// apply_config() takes the information in a configData_t record and applies
// it to the system molecules. The configData_t records have an array of molecules,
// each of which has an array of x, y and z coordinates for every site in the molecule.
// This function loops through said molecules, and for each, loops through each site/atom
// transferring the Cartesian coords in the config array to the system atoms, which are
// located in a linked list, one list for each molecule, in a linked list of molecules.

void apply_config(system_t *system, configData_t *config) {
    char out[1024];
    int i, j;

    molecule_t *moleculePtr = system->molecules;
    atom_t *atomPtr;

    // Loop through each molecule in the system. config->nMolecules should be indentic
    // with the number of molecules found in the "system->molecules" linked list

    for (i = 0; i < config->nMolecules; i++) {
        // if the molecule linked list ends while there is configuration data left to process
        // some kind of error has occurred. Likely, the user has specified a PQR file for one
        // model, and has specified config data for a different model.

        if (moleculePtr == 0) {
            if (i == 0) {
                error(
                    "PQR file contains no molecules.");
            } else {
                error(
                    "PQR file has fewer molecules than (at least) one of the configurations.\n");
                sprintf(out,
                        "Please verify that the config file(s) corresponds to the PQR file. (%d)\n", __LINE__);
                error(out);
            }
            die(1);
        }

        // Find the average coordinate, the "center" of the molecule.
        // This is the point that will form the vector along which perturbations in moveable sites will occur.
        /*
		int N=0;
		
		for( atomPtr = moleculePtr->atoms; atomPtr; atomPtr = atomPtr->next ) {
			center[0] += atomPtr->pos[0];
			center[1] += atomPtr->pos[1];
			center[2] += atomPtr->pos[2];
			N++;
		}
		center[0]/=(double)N;
		center[1]/=(double)N;
		center[2]/=(double)N;
*/

        // Retrieve the intial Center of Mass Coordinate for easy referencing
        double center[3] = {0, 0, 0};
        center[0] = moleculePtr->iCOM[0];
        center[1] = moleculePtr->iCOM[1];
        center[2] = moleculePtr->iCOM[2];
#ifdef DEBUG
        printf(
            "Center of mass of current config is:   %0.5lf   %0.5lf   %0.5lf\n", center[0], center[1], center[2]);
#endif

#ifdef DEBUG
        {
            int i;
            molecule_t *mPtr;
            for (i = 0, mPtr = system->molecules; mPtr; mPtr = mPtr->next, i++) {
                atom_t *aPtr;
                printf(
                    "\tMOLECULE %d\n", i);
                for (aPtr = mPtr->atoms; aPtr; aPtr = aPtr->next) {
                    printf(
                        "\t\t%8s: %0.3lf  %0.3lf  %0.3lf   sig %0.5lf omega %0.5lf\n", aPtr->atomtype, aPtr->pos[0], aPtr->pos[1], aPtr->pos[2], aPtr->sigma, aPtr->omega);
                }
            }
        }
#endif

        // Now loop through all the sites in each molecule. Here, each molecule's nSites variable
        // should coincide with the number of nodes in the system->molecule->atoms linked list.
        atomPtr = moleculePtr->atoms;
        for (j = 0; j < config->molecule[i].nSites; j++) {
            // If we reach the end of the "atoms" linked list, and the config data has info for more
            // sites, then some kind of error has occurred. Likely PQR/config file mismatch

            if (atomPtr == 0) {
                error(
                    "PQR file has more sites than (at least) one of the configurations.\n");
                sprintf(out,
                        "Please verify that the config file(s) corresponds to the PQR file. (%d)\n", __LINE__);
                error(out);
                die(1);
            }

            // Transfer the coords from the config record to the atom_t node

            if (config->molecule[i].type[j] == FIXED) {
                // For fixed sites, this is a direct transfer
                atomPtr->pos[0] = config->molecule[i].x[j];
                atomPtr->pos[1] = config->molecule[i].y[j];
                atomPtr->pos[2] = config->molecule[i].z[j];

            } else {
                // Calculate the length of the position vector of the moveable site relative to the center (avg coordinate) of the molecule
                double scale = sqrt((atomPtr->pos[0] - center[0]) * (atomPtr->pos[0] - center[0]) + (atomPtr->pos[1] - center[1]) * (atomPtr->pos[1] - center[1]) + (atomPtr->pos[2] - center[2]) * (atomPtr->pos[2] - center[2]));
#ifdef DEBUG
                printf(
                    "For ATOM %s--->    scale is %lf\n", atomPtr->atomtype, scale);
                printf(
                    "center is : %lf   %lf   %lf\n", center[0], center[1], center[2]);
                printf(
                    "position  : %lf   %lf   %lf\n", atomPtr->pos[0], atomPtr->pos[1], atomPtr->pos[2]);
#endif

                // Lengthen/contract the unit vector from the config file so that it will be the same length as the
                // site which moved.

                atomPtr->pos[0] = config->molecule[i].iCOM[0] + config->molecule[i].x[j] * scale;
                atomPtr->pos[1] = config->molecule[i].iCOM[1] + config->molecule[i].y[j] * scale;
                atomPtr->pos[2] = config->molecule[i].iCOM[2] + config->molecule[i].z[j] * scale;
            }

            atomPtr = atomPtr->next;

        }  // sites

        // If the config data is used up, but there are more system sites left, some kind of error
        // has occurred. Likely PQR/config file mismatch
        if (atomPtr != 0) {
            error(
                "PQR file has more sites than (at least) one of the configurations.\n");
            sprintf(out,
                    "Please verify that the config file(s) corresponds to the PQR file. (%d)\n", __LINE__);
            error(out);
            die(1);
        }

        // Store the intial Center of Mass Coordinate for this configuration into the system
        // so that the iCOM for this configuration will be there when we need to do this again.
        moleculePtr->iCOM[0] = config->molecule[i].iCOM[0];
        moleculePtr->iCOM[1] = config->molecule[i].iCOM[1];
        moleculePtr->iCOM[2] = config->molecule[i].iCOM[2];

        moleculePtr = moleculePtr->next;

    }  // molecules

    // If the config data is exhausted, but there are still molecules remaining in the system
    // then some kind of error has occurred, likely a PQR/config mismatch
    if (moleculePtr != 0) {
        error(
            "PQR file has more molecules than (at least) one of the configurations.\n");
        sprintf(out,
                "Please verify that the config file(s) corresponds to the PQR file. (%d)\n", __LINE__);
        error(out);
        die(1);
    }

#ifdef DEBUG
    dumpconfigs(system);
#endif

    //output_xyz( system ); // Testing, for making sure coords applied to system match the config data
}

void get_configuration_energies(system_t *system, configData_t *configurations, int nConfigs, param_g *params, configData_t *orig) {
    int i;

    for (i = 0; i < nConfigs; i++) {
        // Orient the system in the configuration described by the ith configData record
        apply_config(system, &configurations[i]);

        // Calculate and store the energy value for this configuration w/the current parameters
        configurations[i].currentFitEnergy = surface_energy(system, ENERGY_TOTAL);

#ifdef DEBUG
        char out[1024];
        sprintf(out,
                "energy: %15.9lf\n", configurations[i].currentFitEnergy);
        output(out);
#endif
    }

    // restore system to initial configuration (orientation)
    apply_config(system, orig);

    return;
}

void relocate_to_origin(molecule_t *molecules) {
    molecule_t *mptr;
    atom_t *aptr;
    int i;

    for (mptr = molecules; mptr; mptr = mptr->next) {
        for (aptr = mptr->atoms; aptr; aptr = aptr->next) {
            for (i = 0; i < 3; i++) {
                aptr->pos[i] -= mptr->com[i];
            }
        }
    }

    return;
}

// Frees all the dynamically allocated memory used by a configData_t record
// IF UNINITIALIZED POINTERS DO NOT HAPPEN TO BE ZERO, this routine will assume
// they are valid addresses, attempt to free(them) and seg-fault.

void delete_config(configData_t *config) {
    if (config->molecule) {
        // There is an address here, delete the information
        int i;
        for (i = 0; i < config->nMolecules; i++) {
            if (config->molecule[i].x) {
                free(config->molecule[i].x);
                config->molecule[i].x = 0;
            }
            if (config->molecule[i].y) {
                free(config->molecule[i].y);
                config->molecule[i].y = 0;
            }
            if (config->molecule[i].z) {
                free(config->molecule[i].z);
                config->molecule[i].z = 0;
            }
            if (config->molecule[i].type) {
                free(config->molecule[i].type);
                config->molecule[i].type = 0;
            }
            if (config->molecule[i].id) {
                int j;
                for (j = 0; j < config->molecule[i].nSites; j++) {
                    if (config->molecule[i].id[j]) {
                        free(config->molecule[i].id[j]);
                        config->molecule[i].id[j] = 0;
                    }
                }
                free(config->molecule[i].id);
            }
        }
    }
    free(config->molecule);
    return;
}

void store_system_config(system_t *system, configData_t *originalConfig) {
    molecule_t *moleculePtr;  // for traversing system's linked list of molecules
    atom_t *atomPtr;          // for traversing system atoms' linked lists of atoms
    int nMolecules = 0;       // a system count for the number of molecules
    int *nSites;              // an array that holds a molecule-by-molecule count of the number of atoms in each molecule
    char out[1 << 10];        // a friendly place for composing error messages

    int i;  // generic counter

    // Count the number of molecules in the system

    for (moleculePtr = system->molecules; moleculePtr; moleculePtr = moleculePtr->next) {
        nMolecules++;
    }
    if (!nMolecules) {
        error(
            "ERROR: Attempting to store system configuration, but there are no molecules in the system.");
        die(-1);
    }

    originalConfig->nMolecules = nMolecules;

    // Allocate an array to hold the molecular config records--one record per system molecule.

    nSites = calloc(nMolecules, sizeof(int));
    memnullcheck(nSites, nMolecules * sizeof(int), __LINE__, __FILE__);

    // Count the number of sites in each molecule and store that number in the corresponding entry of nSites[]
    moleculePtr = system->molecules;
    for (i = 0; i < nMolecules; i++) {
        for (atomPtr = moleculePtr->atoms; atomPtr; atomPtr = atomPtr->next)
            ++nSites[i];

        moleculePtr = moleculePtr->next;
    }

    for (i = 0; i < nMolecules; i++) {
        if (nSites[i] == 0) {
            sprintf(out,
                    "ERROR: System molecule %d has no sites or atoms.\n", i);
            error(out);
            die(-1);
        }
    }

    // Clear any old information from the originalConfig record the function was passed
    delete_config(originalConfig);

    // Allocate memory for the moleculeCoords records (meta-data and coordinate info for each molecule)
    originalConfig->molecule = calloc(nMolecules, sizeof(moleculeCoords));
    memnullcheck(originalConfig->molecule, nMolecules * sizeof(moleculeCoords), __LINE__, __FILE__);

    moleculePtr = system->molecules;
    for (i = 0; i < nMolecules; i++) {
        originalConfig->molecule[i].nSites = nSites[i];

        // Allocate space to hold the coordinate and type data for molecule i
        originalConfig->molecule[i].x = calloc(nSites[i], sizeof(double));
        originalConfig->molecule[i].y = calloc(nSites[i], sizeof(double));
        originalConfig->molecule[i].z = calloc(nSites[i], sizeof(double));
        originalConfig->molecule[i].type = calloc(nSites[i], sizeof(int));  // If calloc() here is changed, note that this needs to be initialized to 0
        memnullcheck(originalConfig->molecule[i].x, nSites[i] * sizeof(double), __LINE__ - 4, __FILE__);
        memnullcheck(originalConfig->molecule[i].y, nSites[i] * sizeof(double), __LINE__ - 4, __FILE__);
        memnullcheck(originalConfig->molecule[i].z, nSites[i] * sizeof(double), __LINE__ - 4, __FILE__);
        memnullcheck(originalConfig->molecule[i].type, nSites[i] * sizeof(int), __LINE__ - 4, __FILE__);

        int j;
        atomPtr = moleculePtr->atoms;
        for (j = 0; j < nSites[i]; j++) {
            // Transfer the atomic/site coordinates (and type code) from
            originalConfig->molecule[i].x[j] = atomPtr->pos[0];
            originalConfig->molecule[i].y[j] = atomPtr->pos[1];
            originalConfig->molecule[i].z[j] = atomPtr->pos[2];
            // The following assignment should typically be made via a bit-wise OR, to preserve other bit
            // codes, butsince we just calloc()'d the space for this item ^^, it will be initialized to 0,
            // so there are no bits to preserve.
            originalConfig->molecule[i].type[j] = determine_mobility_code_from_atom_type(atomPtr->atomtype);

            atomPtr = atomPtr->next;
        }
        moleculePtr = moleculePtr->next;
    }

    // Calculate the current center of mass for the molecule and store
    // these values as the initial Center of Mass coordinates for the config.

    update_com(system->molecules);
    int m;
    moleculePtr = system->molecules;
    for (m = 0; m < originalConfig->nMolecules; m++) {
        originalConfig->molecule[m].iCOM[0] = moleculePtr->com[0];
        originalConfig->molecule[m].iCOM[1] = moleculePtr->com[1];
        originalConfig->molecule[m].iCOM[2] = moleculePtr->com[2];

        moleculePtr = moleculePtr->next;
    }

    // Scale the moveable sites such that they become unit vectors
    // These values will be scaled and added to the molecular center to
    // produce the actual atomic (or site) coordinates

    for (m = 0; m < originalConfig->nMolecules; m++) {
        int a;
        for (a = 0; a < originalConfig->molecule[m].nSites; a++) {
            if (originalConfig->molecule[m].type[a] == MOVEABLE) {
                double dx = originalConfig->molecule[m].x[a] - originalConfig->molecule[m].iCOM[0];
                double dy = originalConfig->molecule[m].y[a] - originalConfig->molecule[m].iCOM[1];
                double dz = originalConfig->molecule[m].z[a] - originalConfig->molecule[m].iCOM[2];

                double r = sqrt(dx * dx + dy * dy + dz * dz);

                originalConfig->molecule[m].x[a] = dx / r;
                originalConfig->molecule[m].y[a] = dy / r;
                originalConfig->molecule[m].z[a] = dz / r;
            }
        }
    }

    free(nSites);
}

// Given a c-string representing an atom ID, this function will return FIXED
// if the site cannot move, or MOVEABLE if it is allowed to float.

int determine_mobility_code_from_atom_type(char *atomType) {
    char linebuf[MAXLINE];

    if (!strncasecmp(atomType,
                     "H2G", 3))
        return FIXED;
    if (!strncasecmp(atomType,
                     "N2G", 3))
        return FIXED;

    if (!strncasecmp(atomType,
                     "H2E", 3))
        return FIXED;
    if (!strncasecmp(atomType,
                     "N2E", 3))
        return FIXED;

    if (!strncasecmp(atomType,
                     "H2N", 3))
        return MOVEABLE;
    if (!strncasecmp(atomType,
                     "N2N", 3))
        return MOVEABLE;

    // Default case...
    sprintf(linebuf,
            "warning: I don't know whether atomType %.100s is MOVEABLE or FIXED. Defaulting to MOVEABLE.\n", atomType);
    error(linebuf);
    return MOVEABLE;
}

// at the end of the run, print (and write to file) the fit energy vs true energy
void dumpbestfitenergies(system_t *system, int nConfigs, configData_t *configuration) {
    char buffer[MAXLINE];
    strcpy(buffer, system->job_name);
    strcat(buffer,
           "-fit.out");
    FILE *ffitout = fopen(buffer,
                          "w");

    int j;
    fprintf(ffitout,
            "%20s\t%20s\n",
            "#abinit_energy",
            "#fit_energy");
    sprintf(buffer,
            "%20s\t%20s\n",
            "#abinit_energy",
            "#fit_energy");
    output(buffer);

    for (j = 0; j < nConfigs; j++) {
        fprintf(ffitout,
                "%20.12lf\t%20.12lf\n", configuration[j].abInitioEnergy, configuration[j].bestFitEnergy);
        sprintf(buffer,
                "%20.12lf\t%20.12lf\n", configuration[j].abInitioEnergy, configuration[j].bestFitEnergy);
        output(buffer);
    }

    fclose(ffitout);

    return;
}

// fit potential energy parameters via simulated annealing
int surface_fit_arbitrary(system_t *system) {
    // override default simulated annealing values, if they were specified in the input file
    double temperature = ((system->fit_start_temp) ? system->fit_start_temp : TEMPERATURE);
    double max_energy = ((system->fit_max_energy) ? system->fit_max_energy : MAX_ENERGY);
    double schedule = ((system->fit_schedule) ? system->fit_schedule : SCHEDULE);
    system->fit_max_energy = max_energy;

    // used only if qshift is on
    double quadrupole = 0;

    int i, j,  // generic counters
        nSteps;

    param_t *param_ptr;
    param_g *params;

    double current_error, last_error, global_minimum;
    qshiftData_t *qshiftData = NULL;  //used only with qshift

    // Store the initial configuration
    configData_t originalConfig;
    originalConfig.id = 0;        // This line and the next are there to prevent a seg-fault when
    originalConfig.molecule = 0;  // the delete_config() is called via store_system_config()...

    store_system_config(system, &originalConfig);

    // Construct the array that contains the configuration information,
    // the associated energies, and get the total number of configs
    int nConfigs = 0;
    configData_t *configuration = 0;
    read_config_fit_input_file(system, &configuration, &nConfigs);

    // Output the geometry for visual verification
    // output_pqrs ( system, nCurves, curve );

    // record parameters. we will determine which ones need to be adjusted later
    params = record_params(system);

    // Obtain initial model-calculated energy for each configuration
    get_configuration_energies(system, configuration, nConfigs, params, &originalConfig);

    // Calculate the error of all the configurations vs the ab initio energy
    current_error = error_calc_fr_arbitrary_configurations(system, nConfigs, configuration, max_energy);

#ifdef DEBUG
    printf(
        "Error: %0.5lf\n", current_error);
#endif

    // print some header info
    for (param_ptr = params->type_params; param_ptr; param_ptr = param_ptr->next)
        printf(
            "SURFACE: Atom type: %d @ %s\n", param_ptr->ntypes, param_ptr->atomtype);
    if (!system->fit_boltzmann_weight)
        printf(
            "*** any input energy values greater than %f K will not contribute to the fit ***\n", max_energy);

    // write params and current_error to stdout
    printf(
        "*** Initial Fit: \n");
    output_params(temperature, current_error, params);
    printf(
        "*****************\n");

    // set the minimum error we've found
    global_minimum = current_error;
    last_error = current_error;
    for (j = 0; j < nConfigs; j++)
        configuration[j].bestFitEnergy = configuration[j].currentFitEnergy;

    // initialize stuff for qshift
    if (system->surf_qshift_on) {
        qshiftData = malloc(sizeof(qshiftData_t));
        quadrupole = calcquadrupole(system);
    }

    // ANNEALING
    // Loop over temperature. When the temperature reaches MIN_TEMPERATURE
    // we quit. Assuming we find an error smaller than the initial error
    // we'll spit out a fit_geometry.pqr file, and you'll want to re-anneal.

    for (nSteps = 0; temperature > MIN_TEMPERATURE; ++nSteps) {
        // move molecules back to origin
        update_com(system->molecules);
        relocate_to_origin(system->molecules);

        // randomly perturb the parameters
        //surf_perturb_arbitrary_configs ( system, quadrupole, qshiftData, params );
        surf_perturb(system, quadrupole, qshiftData, params);

        // apply the trial parameters
        surface_dimer_parameters(system, params);
        get_configuration_energies(system, configuration, nConfigs, params, &originalConfig);

        // calc error
        current_error = error_calc_fr_arbitrary_configurations(system, nConfigs, configuration, max_energy);

// Remove when working properly
#ifdef DEBUG
        printf(
            "Error: %0.8lf\n", current_error);
#endif

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
                    new_global_min_arbitrary_configs(system, nConfigs, configuration);
                }
            }
        } else {
            // MOVE REJECTED ////
            revert_parameters(system, params);
        }  //END DO MONTE CARLO

        // decrement the temperature
        temperature = temperature * schedule;  // schedule
    }                                          //END ANNEALING

    dumpconfigs(system);
    dumpbestfitenergies(system, nConfigs, configuration);

    // Return all memory back to the system
    ////////////////////////////////////////////////////////////////////////////
    free(qshiftData);
    delete_config(&originalConfig);
    for (i = 0; i < nConfigs; i++) delete_config(&(configuration[i]));
    free(configuration);
    // Free the param linked list
    param_t *next;
    param_ptr = params->type_params;
    while (param_ptr) {
        next = param_ptr->next;
        free(param_ptr);
        param_ptr = next;
    }
    free(params);
    ////////////////////////////////////////////////////////////////////////////

    return (0);
}

// Reads a file of system configurations, and creates an array of configData_t records and stores the
// configuration information from the input file therein.

void read_config_fit_input_file(system_t *system, configData_t **c, int *nConfigs) {
    FILE *fp_fit;
    char linebuf[MAXLINE],  // Line buffer for reading file input
        token1[MAXLINE],    // Variables for parsing parameters in input files.
        token2[MAXLINE],
        token3[MAXLINE],
        token4[MAXLINE],
        token5[MAXLINE],
        errMsg[MAXLINE];  // Output buffer for error messages

    // Count the number of molecules, and the number of sites in each molecule
    // about which the system is aware. This number should match the data in
    // the config files, and will be used as the standard against which the
    // integrity of the files are judged.

    molecule_t *moleculePtr;
    atom_t *atomPtr;
    int nMolecules = 0;
    int *nSites = 0;
    int MAX_CONFIGS = 5000;
    int a, i, j, m;

    // Count the number of molecules in the system, and then allocate an array with this
    // many elements to hold the number of sites that each corresponding molecule has.

    for (moleculePtr = system->molecules; moleculePtr; moleculePtr = moleculePtr->next)
        ++nMolecules;
    if (nMolecules == 0) {
        error(
            "ERROR: There are no system molecules to which the config data can be applied.\n");
        sprintf(errMsg,
                "ERROR: Please check your input .PQR file. (%d)\n", __LINE__);
        error(errMsg);
        die(-1);
    }
    nSites = calloc(nMolecules, sizeof(int));
    memnullcheck(nSites, nMolecules * sizeof(int), __LINE__, __FILE__);

    // Count the number of sites in each molecule and store that number in nSites[]
    moleculePtr = system->molecules;
    for (i = 0; i < nMolecules; i++) {
        for (atomPtr = moleculePtr->atoms; atomPtr; atomPtr = atomPtr->next)
            ++nSites[i];
        moleculePtr = moleculePtr->next;
    }

    for (i = 0; i < nMolecules; i++) {
        if (nSites[i] == 0) {
            sprintf(errMsg,
                    "ERROR: System molecule %d has no sites or atoms. (%d)\n", i, __LINE__);
            error(errMsg);
            die(-1);
        }
    }

    // Allocate memory to hold the configuration data
    *c = calloc(sizeof(configData_t), MAX_CONFIGS);
    //	configData_t *config = (*c);  // for slightly more convenient referencing
    configData_t *configchk;
    int currentConfig = 0;

    // Open all the user-specified input files and extract their configuration data.
    fileNode_t *fNode = system->fit_input_list.next;
    while (fNode) {
        // Open file for reading
        ////////////////////////////

        char *filename = fNode->data.filename;
        printf(
            "INPUT: Loading %s\n", filename);
        fp_fit = fopen(filename,
                       "r");
        filecheck(fp_fit, filename, READ);
        rewind(fp_fit);

        // Parse and store data found in the input file
        //////////////////////////////////////////////////

        while (fgets(linebuf, MAXLINE, fp_fit))  // read config record header
        {
            memset(token1, 0, MAXLINE);
            memset(token2, 0, MAXLINE);
            memset(token3, 0, MAXLINE);
            memset(token4, 0, MAXLINE);
            memset(token5, 0, MAXLINE);

            if (fgets(linebuf, MAXLINE, fp_fit)) {  // read the energy
                sscanf(linebuf,
                       "%s", token1);  // parse the energy line
                if (safe_atof(token1, &((*c)[currentConfig].abInitioEnergy))) {
                    error(token1);
                    sprintf(errMsg,
                            "error reading config file (%d)\n", __LINE__);
                    error(errMsg);
                    die(-1);
                }
            } else {
                sprintf(errMsg,
                        "Error reading config file  (%d)\n", __LINE__);
                error(errMsg);
                die(-1);
            }
            (*c)[currentConfig].currentFitEnergy = 0;
            (*c)[currentConfig].bestFitEnergy = 0;
            (*c)[currentConfig].nMolecules = nMolecules;
            (*c)[currentConfig].molecule = calloc(sizeof(moleculeCoords), nMolecules);

            m = 0;
            for (moleculePtr = system->molecules; moleculePtr; moleculePtr = moleculePtr->next) {
                // initialize the configuartion record, and allocate space to hold the site types & coordinates

                (*c)[currentConfig].molecule[m].nSites = nSites[m];

                (*c)[currentConfig].molecule[m].type = calloc(sizeof(unsigned int), nSites[m]);
                memnullcheck(nSites, nMolecules * sizeof(int), __LINE__, __FILE__);

                (*c)[currentConfig].molecule[m].x = calloc(sizeof(double), nSites[m]);
                memnullcheck(nSites, nMolecules * sizeof(int), __LINE__, __FILE__);

                (*c)[currentConfig].molecule[m].y = calloc(sizeof(double), nSites[m]);
                memnullcheck(nSites, nMolecules * sizeof(int), __LINE__, __FILE__);

                (*c)[currentConfig].molecule[m].z = calloc(sizeof(double), nSites[m]);
                memnullcheck(nSites, nMolecules * sizeof(int), __LINE__, __FILE__);

                a = 0;
                double x, y, z;

                for (atomPtr = moleculePtr->atoms; atomPtr; atomPtr = atomPtr->next) {
                    if (fgets(linebuf, MAXLINE, fp_fit)) {  // read the energy

                        int moleculeID;

                        sscanf(linebuf,
                               "%s %s %s %s %s", token1, token2, token3, token4, token5);  // parse the coord line
                        if (safe_atoi(token2, &moleculeID)) {
                            sprintf(errMsg,
                                    "error reading config file (molecule ID) (%d)\n", __LINE__);
                            error(errMsg);
                            die(-1);
                        }
                        if (safe_atof(token3, &x)) {
                            sprintf(errMsg,
                                    "error reading config file (molecule x) (%d)\n", __LINE__);
                            error(errMsg);
                            die(-1);
                        }
                        if (safe_atof(token4, &y)) {
                            sprintf(errMsg,
                                    "error reading config file (molecule y) (%d)\n", __LINE__);
                            error(errMsg);
                            die(-1);
                        }
                        if (safe_atof(token5, &z)) {
                            sprintf(errMsg,
                                    "error reading config file (molecule z) (%d)\n", __LINE__);
                            error(errMsg);
                            die(-1);
                        }
                        if (moleculeID != (m + 1)) {
                            sprintf(errMsg,
                                    "Molecule ID number does not match current molecule. (%d)\n", __LINE__);
                            error(errMsg);
                            die(-1);
                        }

                    } else {
                        sprintf(errMsg,
                                "Error reading config file (Ran out of data?) (%d)\n", __LINE__);
                        error(errMsg);
                        die(-1);
                    }

                    // Store the coordinates in the configuration array
                    (*c)[currentConfig].molecule[m].x[a] = x;
                    (*c)[currentConfig].molecule[m].y[a] = y;
                    (*c)[currentConfig].molecule[m].z[a] = z;

                    // Set the moveable/fixed flag
                    (*c)[currentConfig].molecule[m].type[a] = determine_mobility_code_from_atom_type(token1);

                    a++;  // next atom...
                }
                m++;  // next molecule...

            }  // done processing this line, proceed to the next...

            ++currentConfig;
            if (currentConfig == MAX_CONFIGS) {
                MAX_CONFIGS *= 2;
                configchk = realloc(*c, MAX_CONFIGS * sizeof(configData_t));
                memnullcheck(configchk, MAX_CONFIGS * sizeof(configData_t), __LINE__ - 1, __FILE__);
                *c = configchk;
            }

        }  // while (fgets), i.e.  while(!EOF)

        fclose(fp_fit);
        fNode = fNode->next;
    }  // while(fNode) i.e., Process next file in the list...

    *nConfigs = currentConfig;

    // Shrink the array to the appropriate size
    configchk = realloc(*c, (*nConfigs) * sizeof(configData_t));
    memnullcheck(configchk, (*nConfigs) * sizeof(configData_t), __LINE__ - 1, __FILE__);
    *c = configchk;

    // Calculate initial Center of Mass (iCOM) coordinate for each configuration
    //
    // This will loop through each configuration. moleculePtr traverses the system  molecule
    // list as a counter steps through the corresponding moleculeCoord records.
    //   For each molecule, atomPtr will traverse the current molecule's atom list, as a counter
    //   steps through each corresponding element of the moleculeCoord record.
    //
    // The center of mass from the atomic linked list is used in conjunction with the coordinates
    // in the configData records to produce a center of mass. Although the center of mass may change
    // over time, this initial COM (iCOM) will be used for the life of the program, as a fixed point,
    // relative to which a site's positional perturbations can be directed.

    for (i = 0; i < (*nConfigs); i++) {
        moleculePtr = system->molecules;
        for (m = 0; m < (*c)[i].nMolecules; m++) {
            double center[3] = {0, 0, 0};
            double totMass = 0.0;

            atomPtr = moleculePtr->atoms;
            for (j = 0; j < (*c)[i].molecule[m].nSites; j++) {
                double mass = atomPtr->mass;
                totMass += mass;
                center[0] += (*c)[i].molecule[m].x[j] * mass;
                center[1] += (*c)[i].molecule[m].y[j] * mass;
                center[2] += (*c)[i].molecule[m].z[j] * mass;

                atomPtr = atomPtr->next;
            }

            (*c)[i].molecule[m].iCOM[0] = center[0] / totMass;
            (*c)[i].molecule[m].iCOM[1] = center[1] / totMass;
            (*c)[i].molecule[m].iCOM[2] = center[2] / totMass;

            moleculePtr = moleculePtr->next;
        }
    }

    // Scale the moveable sites such that they become unit vectors
    // These values will be scaled and added to the molecular center to
    // produce the actual atomic (or site) coordinates
    int cfg;
    for (cfg = 0; cfg < (*nConfigs); cfg++) {
        for (m = 0; m < (*c)[cfg].nMolecules; m++) {
            for (a = 0; a < (*c)[cfg].molecule[m].nSites; a++) {
                if ((*c)[cfg].molecule[m].type[a] == MOVEABLE) {
                    double dx = (*c)[cfg].molecule[m].x[a] - (*c)[cfg].molecule[m].iCOM[0];
                    double dy = (*c)[cfg].molecule[m].y[a] - (*c)[cfg].molecule[m].iCOM[1];
                    double dz = (*c)[cfg].molecule[m].z[a] - (*c)[cfg].molecule[m].iCOM[2];

                    double r = sqrt(dx * dx + dy * dy + dz * dz);

                    (*c)[cfg].molecule[m].x[a] = dx / r;
                    (*c)[cfg].molecule[m].y[a] = dy / r;
                    (*c)[cfg].molecule[m].z[a] = dz / r;
                }
            }
        }
    }

    /*
	// This was here for debugging purposes....
	
	printf( "\n\n\nConfiguration Data: \n\n" );
	for( int aa = 0; aa< *nConfigs; aa++ ) {
		printf( "Configuration %d\n------------------\n", aa );
		printf( "AbInit Energy: %lf\n",      config[aa].abInitioEnergy );
		printf( "Current Fit Energy: %lf\n", config[aa].currentFitEnergy );
		printf( "Best Fit Energy: %lf\n",    config[aa].bestFitEnergy );
		printf( "Number of Molecules: %d\n", config[aa].nMolecules );
		for( int bb = 0; bb < config[aa].nMolecules; bb++ ) {
			printf( "  Molecule %d (%d)\n", bb, config[aa].molecule[bb].nSites );
			for( int cc = 0; cc<config[aa].molecule[bb].nSites; cc++ ) {
				printf( "   Site %d\n", cc );
				if( config[aa].molecule[bb].type[cc] == FIXED )
					printf( "    type: FIXED\n", config[aa].molecule[bb].x[cc] );
				else if( config[aa].molecule[bb].type[cc] == MOVEABLE )
					printf( "    type: MOVEABLE\n", config[aa].molecule[bb].x[cc] );
				else
					printf( "    type: UNKNOWN\n", config[aa].molecule[bb].x[cc] );
					
				printf( "    x: %lf\n", config[aa].molecule[bb].x[cc] );
				printf( "    y: %lf\n", config[aa].molecule[bb].y[cc] );
				printf( "    z: %lf\n", config[aa].molecule[bb].z[cc] );
			}
		}
	}
	*/

    free(nSites);
    return;
}
