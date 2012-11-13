#include <mc.h>

int setup_simulation_box(FILE * finput, system_t * system) {

	// read in input.pqr molecules
	system->molecules = read_molecules(finput,system);
	if (!system->molecules && system->ensemble == ENSEMBLE_REPLAY ) {
		output("INPUT: end of trajectory file\n");
		return 1;
	}	
	else if (!system->molecules) {
		error("INPUT: error reading in input molecules\n");
		return(-1);
	}
	else 
		output("INPUT: finished reading in molecules\n");

	//read in pqr box and calculate periodic boundary conditions if neccessary
	if (system->read_pqr_box_on) 
		read_pqr_box(finput,system);
		
	pbc(system);
	if ( system->ensemble != ENSEMBLE_SURF && system->ensemble != ENSEMBLE_SURF_FIT ) {
		if ((system->pbc->volume <= 0.0) || (system->pbc->cutoff <= 0.0) ) {
			error("INPUT: invalid simulation box dimensions.\n");
			return -1;;
		}
	}

	// read in the insertion molecules
	if (system->insert_input) {
		system->insertion_molecules = read_insertion_molecules(system);
		if (!system->insertion_molecules) {
			error("INPUT: error read in insertion molecules\n");
			return -1;
		} else
			output("INPUT: finished reading in insertion molecules\n");
	}
	else //else only 1 sorbate type
		system->sorbateCount = 1;

	// now that we've read in the sorbates, we can check that user_fugacities is properly set (if used)
	if ( system->user_fugacities ) {
		if ( system->fugacitiesCount != system->sorbateCount ) {
			error("INPUT: number of fugacities set via user_fugacities does not match the number of sorbates.\n");
			return -1;
		}
	}

	return 0;
}
