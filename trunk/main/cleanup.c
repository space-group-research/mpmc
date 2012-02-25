/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

/* free a molecule and all of it's associated stuff */
void free_molecule(system_t *system, molecule_t *molecule) {

	int i;

	molecule->next = NULL;

#ifdef QM_ROTATION
	if(system->quantum_rotation && !molecule->frozen) {

		free(molecule->quantum_rotational_energies);
		free(molecule->quantum_rotational_eigensymmetry);
		for(i = 0; i < system->quantum_rotation_level_max; i++)
			free(molecule->quantum_rotational_eigenvectors[i]);
		free(molecule->quantum_rotational_eigenvectors);

	}
#endif /* QM_ROTATION */

	free_pairs(molecule);
	free_atoms(molecule);
	free(molecule);

}


void free_pairs(molecule_t *molecules) {

	int i;
	pair_t **ptr_array;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	pair_t *pair_ptr;

	/* build an array of ptrs to be freed */
	for(molecule_ptr = molecules, i = 0, ptr_array = NULL; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				ptr_array = realloc(ptr_array, sizeof(pair_t *)*(i + 1));
				memnullcheck(ptr_array,sizeof(pair_t *)*(i+1),34);
				ptr_array[i] = pair_ptr;
				++i;

			}
		}
	}

	/* free the whole array of ptrs */
	for(--i; i >= 0; i--) free(ptr_array[i]);

	/* zero out the heads */
	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
			atom_ptr->pairs = NULL;

	/* free our temporary array */
	if(ptr_array) free(ptr_array);

}

void free_atoms(molecule_t *molecules) {

	int i, n;
	atom_t **ptr_array;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;

	/* build the ptr array */
	for(molecule_ptr = molecules, i = 0, n = 0, ptr_array = NULL; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			ptr_array = realloc(ptr_array, sizeof(atom_t *)*(i + 1));
			memnullcheck(ptr_array,sizeof(atom_t *)*(i+1),35);
			ptr_array[i] = atom_ptr;
			++i, ++n;
		}
	}

	/* free the whole array of ptrs */
	for(i = 0; i < n; i++) free(ptr_array[i]);

	/* free our temporary array */
	free(ptr_array);

}

void free_molecules(molecule_t *molecules) {

	int i;
	molecule_t **ptr_array = NULL;
	molecule_t *molecule_ptr;

	/* build the ptr array */
	for(molecule_ptr = molecules, i = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		ptr_array = realloc(ptr_array, sizeof(molecule_t *)*(i + 1));
		memnullcheck(ptr_array,sizeof(molecule_t *)*(i + 1),36);
		ptr_array[i] = molecule_ptr;
		++i;
	}

	/* free the whole array of ptrs */
	for(--i; i >= 0; i--) free(ptr_array[i]);

	/* free our temporary array */
	free(ptr_array);

}

/* free the polarization matrices */
void free_matrices(system_t *system) {

	int i, N;

	N = 3*system->checkpoint->N_atom;

	if ( system->ensemble == ENSEMBLE_SURF_FIT ) {
		//N_atom checkpoint is not set
		//Need to manually count
		atom_t * atom_ptr;
		molecule_t * mole_ptr;
		N=0;
		for (mole_ptr = system->molecules; mole_ptr; mole_ptr = mole_ptr->next ) {
			for (atom_ptr = mole_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {
				N+=3;
			}
		}
	}

	for(i = 0; i < N; i++) {
		free(system->A_matrix[i]);
		if(!system->polar_iterative)
			free(system->B_matrix[i]);
	}

	free(system->A_matrix);
	free(system->B_matrix);

}

/* free the cavity bias insertion grid */
void free_cavity_grid(system_t *system) {

	int i, j, N;

	N = system->cavity_grid_size;

	for(i = 0; i < N; i++) {

		for(j = 0; j < N; j++)
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

	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {

		if(!molecule_ptr->frozen) {

			free(molecule_ptr->quantum_rotational_energies);
			free(molecule_ptr->quantum_rotational_eigensymmetry);
			for(i = 0; i < system->quantum_rotation_level_max; i++)
				free(molecule_ptr->quantum_rotational_eigenvectors[i]);
			free(molecule_ptr->quantum_rotational_eigenvectors);

		}

	}
}
#endif /* QM_ROTATION */

/* free all of our data structures */
void cleanup(system_t *system) {

	int i,j;

#ifdef QM_ROTATION
	if(system->quantum_rotation) free_rotational(system);
#endif /* QM_ROTATION */
	if(system->polarization && !system->cuda) free_matrices(system);

	free_pairs(system->molecules);
	free_atoms(system->molecules);
	free_molecules(system->molecules);

	free(system->pdb_input);
	free(system->pdb_output);
	free(system->energy_output);

	if(system->traj_output) free(system->traj_output);
	if(system->dipole_output) free(system->dipole_output);
	if(system->field_output) free(system->field_output);

	if(system->frozen_output) free(system->frozen_output);


	if(system->cavity_bias) free_cavity_grid(system);

	/* free statistics */
	free(system->nodestats);
	free(system->avg_nodestats);
	free(system->observables);
	free(system->avg_observables);
	free(system->checkpoint->observables);
	if ( system->checkpoint->molecule_backup != NULL ) 
		free_molecule(system, system->checkpoint->molecule_backup);
	free(system->checkpoint);
	
	/*free histogram stuff*/
	for ( i=0; i<system->grids->histogram->x_dim; i++ ) {
		for ( j=0; j<system->grids->histogram->y_dim; j++ ) {
			free(system->grids->histogram->grid[i][j]);
		}
		free(system->grids->histogram->grid[i]);
	}
	for ( i=0; i<system->grids->histogram->x_dim; i++ ) {
		for ( j=0; j<system->grids->histogram->y_dim; j++ ) {
			free(system->grids->avg_histogram->grid[i][j]);
		}
		free(system->grids->avg_histogram->grid[i]);
	}
	free(system->grids->histogram->grid);
	free(system->grids->avg_histogram->grid);
	free(system->grids->avg_histogram);
	free(system->grids->histogram);
	free(system->grids);
	free(system->histogram_output);

	/* free fit input lists*/
	if ( system->fit_input_list.data.count > 0 ) {
		fileNode_t *node = &(system->fit_input_list);
		fileNode_t *nextnode;
	//the first one is statically declared
	//the others are malloc'd in input.c
		if ( node->next ) node=node->next;
		while ( node->next ) {
			nextnode=node->next;
			free(node->data.filename);
			free(node);
			node=nextnode;
		}
		free(node->data.filename);
		free(node);
	}


	free(system->pdb_restart);

	free(system->pbc);

	free(system);

}

/* on SIGTERM, cleanup and exit */
void terminate_handler(int sigtype, system_t *sys_ptr) {

	static system_t *system;

	switch(sigtype) {

		case SIGTERM:
			output("CLEANUP: ************ SIGTERM received, exiting *************\n");
			close_files(sys_ptr);
			cleanup(system);
#ifdef MPI
			if(!rank) close_files(sys_ptr);
			MPI_Finalize();
#else
			exit(1);
#endif /* MPI */
			break;

		case SIGUSR1:
			output("CLEANUP: ************ SIGUSR1 received, exiting *************\n");
			close_files(sys_ptr);
			cleanup(system);
#ifdef MPI
			if(!rank) close_files(sys_ptr);
			MPI_Finalize();
#else
			exit(1);
#endif /* MPI */
			break;

		case SIGUSR2:
			output("CLEANUP: ************ SIGUSR2 received, exiting *************\n");
			close_files(sys_ptr);
			cleanup(system);
#ifdef MPI
			if(!rank) close_files(sys_ptr);
			MPI_Finalize();
#else
			exit(1);
#endif /* MPI */
			break;
		case -1: /* install the static ptr */
			system = sys_ptr;
			break;

	}
	

}


