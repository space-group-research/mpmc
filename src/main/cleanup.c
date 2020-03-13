/* 

Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

void free_my_pairs(molecule_t *molecule) {
    int i = 0;
    pair_t **ptr_array = NULL;
    atom_t *aptr;
    pair_t *pptr;

    //build an array of pairs
    for (aptr = molecule->atoms; aptr; aptr = aptr->next) {
        for (pptr = aptr->pairs; pptr; pptr = pptr->next) {
            ptr_array = realloc(ptr_array, sizeof(pair_t *) * (i + 1));
            memnullcheck(ptr_array, sizeof(pair_t *) * (i + 1), __LINE__ - 1, __FILE__);
            ptr_array[i] = pptr;
            ++i;
        }
    }

    //free the pairs
    for (--i; i >= 0; i--) free(ptr_array[i]);

    //zero out the heads
    for (aptr = molecule->atoms; aptr; aptr = aptr->next)
        aptr->pairs = NULL;

    //free the temp array
    free(ptr_array);

    return;
}

void free_my_atoms(molecule_t *molecules) {
    int i = 0;
    atom_t *aptr;
    atom_t **aarray = NULL;

    //build an array of atoms
    for (aptr = molecules->atoms; aptr; aptr = aptr->next) {
        aarray = realloc(aarray, sizeof(atom_t *) * (i + 1));
        memnullcheck(aarray, sizeof(atom_t *) * (i + 1), __LINE__ - 1, __FILE__);
        aarray[i] = aptr;
        i++;
    }

    //free the atoms
    while (i--)
        free(aarray[i]);

    //free the temp array
    free(aarray);

    return;
}

/* free a molecule and all of it's associated stuff */
void free_molecule(system_t *system, molecule_t *molecule) {
#ifdef QM_ROTATION
    int i;
    if (system->quantum_rotation && !molecule->frozen) {
        free(molecule->quantum_rotational_energies);
        free(molecule->quantum_rotational_eigensymmetry);
        for (i = 0; i < system->quantum_rotation_level_max; i++)
            free(molecule->quantum_rotational_eigenvectors[i]);
        free(molecule->quantum_rotational_eigenvectors);
    }
#endif /* QM_ROTATION */

    //free pairs belonging to this molecule only
    free_my_pairs(molecule);
    free_my_atoms(molecule);
    free(molecule);
}

// free all pairs
void free_all_pairs(system_t *system) {
    int i, j;
    pair_t *pair_ptr;
    pair_t **ptr_array = NULL;

    j = 0;  //pair array index
    //build an array of all pair pointers
    for (i = 0; i < system->natoms; i++) {
        for (pair_ptr = system->atom_array[i]->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
            ptr_array = realloc(ptr_array, sizeof(pair_t *) * (j + 1));
            memnullcheck(ptr_array, sizeof(pair_t *) * (j + 1), __LINE__ - 1, __FILE__);
            ptr_array[j] = pair_ptr;
            j++;
        }
    }

    /* free the whole array of ptrs */
    while (j--) free(ptr_array[j]);

    /* zero out the heads */
    for (i = 0; i < system->natoms; i++)
        system->atom_array[i]->pairs = NULL;

    /* free our temporary array */
    if (ptr_array) free(ptr_array);
}

// free all molecules
void free_all_molecules(system_t *system, molecule_t *molecules) {
    int i;
    molecule_t **marray = NULL;
    molecule_t *mptr;

    /* build the ptr array */  // we can't free all of system->molecule_array since it's redundant
    for (mptr = molecules, i = 0; mptr; mptr = mptr->next) {
        free_my_atoms(mptr);
        marray = realloc(marray, sizeof(molecule_t *) * (i + 1));
        memnullcheck(marray, sizeof(molecule_t *) * (i + 1), __LINE__ - 1, __FILE__);
        marray[i] = mptr;
        i++;
    }

    /* free the whole array of ptrs */
    for (--i; i >= 0; i--) free(marray[i]);

    /* free our temporary arrays */
    free(marray);

    return;
}

/* free the polarization matrices */
void free_matrices(system_t *system) {
    int i, N;

    if (!system->A_matrix && !system->B_matrix)
        return;  //nothing to do

    if (system->checkpoint->thole_N_atom)
        N = 3 * system->checkpoint->thole_N_atom;
    else
        N = 3 * system->natoms;

    for (i = 0; i < N; i++) {
        free(system->A_matrix[i]);
        if (!system->polar_iterative)
            free(system->B_matrix[i]);
    }

    free(system->A_matrix);
    free(system->B_matrix);
    system->A_matrix = system->B_matrix = NULL;

    return;
}

/* free the cavity bias insertion grid */
void free_cavity_grid(system_t *system) {
    int i, j, N;
    N = system->cavity_grid_size;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++)
            free(system->cavity_grid[i][j]);
        free(system->cavity_grid[i]);
    }
    free(system->cavity_grid);
}

#ifdef QM_ROTATION
/* free structures associated with quantum rotations */
void free_rotational(system_t *system) {
    int i;
    molecule_t *molecule_ptr;

    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        if (!molecule_ptr->frozen) {
            free(molecule_ptr->quantum_rotational_energies);
            free(molecule_ptr->quantum_rotational_eigensymmetry);
            for (i = 0; i < system->quantum_rotation_level_max; i++)
                free(molecule_ptr->quantum_rotational_eigenvectors[i]);
            free(molecule_ptr->quantum_rotational_eigenvectors);
        }
    }
}
#endif /* QM_ROTATION */

//free vdw pointer which keeps track of e_iso energies
void free_vdw_eiso(vdw_t *vdw_eiso_info) {
    vdw_t *vp;
    vdw_t **varray = NULL;
    int i = 0;

    for (vp = vdw_eiso_info; vp; vp = vp->next) {
        varray = realloc(varray, sizeof(vdw_t *) * (i + 1));
        memnullcheck(varray, sizeof(vdw_t *) * (i + 1), __LINE__ - 1, __FILE__);
        varray[i] = vp;
        i++;
    }

    while (i > 0) {
        i--;
        free(varray[i]);
    }

    free(varray);

    return;
}

/* free all of our data structures */
void cleanup(system_t *system) {
    int i, j;

#ifdef QM_ROTATION
    if (system->quantum_rotation) free_rotational(system);
#endif /* QM_ROTATION */
    if (system->polarization && !system->cuda) free_matrices(system);

    if (system->polar_wolf_alpha_lookup)
        if (system->polar_wolf_alpha_table)
            free(system->polar_wolf_alpha_table);

    //need to rebuild atom and pair arrays so we can free everything
    system->natoms = countNatoms(system);
    rebuild_arrays(system);

    free_all_pairs(system);
    free_all_molecules(system, system->molecules);

    //free our arrays
    free(system->molecule_array);
    free(system->atom_array);

    free(system->pqr_input);
    free(system->pqr_output);
    free(system->energy_output);
    free(system->energy_output_csv);
    free(system->virial_output);

    if (system->surf_output) free(system->surf_output);
    if (system->traj_output) free(system->traj_output);
    if (system->traj_input) free(system->traj_input);
    if (system->dipole_output) free(system->dipole_output);
    if (system->field_output) free(system->field_output);
    if (system->frozen_output) free(system->frozen_output);
    if (system->surf_preserve_rotation_on) free(system->surf_preserve_rotation_on);
    if (system->cavity_bias) free_cavity_grid(system);

    if (system->surf_do_not_fit_list != NULL) {
        for (i = 0; i < 20; i++)
            free(system->surf_do_not_fit_list[i]);
        free(system->surf_do_not_fit_list);
    }

    if (system->vdw_eiso_info) free_vdw_eiso(system->vdw_eiso_info);

    // free multi sorbate related stuff
    if (system->fugacities)
        free(system->fugacities);
    if (system->insert_input) free(system->insert_input);
    //insert.pdb arrays and shit
    if (system->insertion_molecules_array) free(system->insertion_molecules_array);
    if (system->insertion_molecules) free_all_molecules(system, system->insertion_molecules);
    // free sorbate info array
    free(system->sorbateInfo);

    /* free statistics */
    free(system->nodestats);
    free(system->avg_nodestats);
    free(system->observables);
    free(system->avg_observables);
    free(system->checkpoint->observables);
    if (system->checkpoint->molecule_backup != NULL)
        free_molecule(system, system->checkpoint->molecule_backup);
    free(system->checkpoint);

    /*free histogram stuff*/
    for (i = 0; i < system->grids->histogram->x_dim; i++) {
        for (j = 0; j < system->grids->histogram->y_dim; j++) {
            free(system->grids->histogram->grid[i][j]);
        }
        free(system->grids->histogram->grid[i]);
    }
    if (!rank) {
        for (i = 0; i < system->grids->histogram->x_dim; i++) {
            for (j = 0; j < system->grids->histogram->y_dim; j++) {
                free(system->grids->avg_histogram->grid[i][j]);
            }
            free(system->grids->avg_histogram->grid[i]);
        }
        free(system->grids->avg_histogram->grid);
    }
    free(system->grids->histogram->grid);
    free(system->grids->avg_histogram);
    free(system->grids->histogram);
    free(system->grids);
    free(system->histogram_output);

    /* free fit input lists*/
    if (system->fit_input_list.data.count > 0) {
        fileNode_t *node = &(system->fit_input_list);
        fileNode_t *nextnode;
        //the first one is statically declared
        //the others are malloc'd in input.c
        if (node->next) node = node->next;
        while (node->next) {
            nextnode = node->next;
            free(node->data.filename);
            free(node);
            node = nextnode;
        }
        free(node->data.filename);
        free(node);
    }

    free(system->pqr_restart);

    free(system->pbc);

    free(system->job_name);  // (CRC)

    free(system);
}

/* on SIGTERM, cleanup and exit */
void terminate_handler(int sigtype, system_t *sys_ptr) {
    static system_t *system;

    switch (sigtype) {
        case SIGTERM:
            output(
                "CLEANUP: ************ SIGTERM received, exiting *************\n");
            close_files(sys_ptr);
            cleanup(system);
#ifdef MPI
            if (!rank) close_files(sys_ptr);
#else
            die(1);
#endif /* MPI */
            break;

#ifndef __WIN32__
        case SIGUSR1:
            output(
                "CLEANUP: ************ SIGUSR1 received, exiting *************\n");
            close_files(sys_ptr);
            cleanup(system);
#endif
#ifdef MPI
            if (!rank) close_files(sys_ptr);
#else
            die(1);
#endif /* MPI */
            break;

#ifndef __WIN32__
        case SIGUSR2:
            output(
                "CLEANUP: ************ SIGUSR2 received, exiting *************\n");
            close_files(sys_ptr);
            cleanup(system);
#endif
#ifdef MPI
            if (!rank) close_files(sys_ptr);
#else
            die(1);
#endif /* MPI */
            break;
        case -1: /* install the static ptr */
            system = sys_ptr;
            break;
    }
}
