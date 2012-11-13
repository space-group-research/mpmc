#include <mc.h>

int calculate_te (system_t *system) {

	int j;
	double initial_energy, final_energy, delta_energy;
	char linebuf[MAXLINE];

	// calculate the energy of the system 
	initial_energy = energy(system);

#ifdef QM_ROTATION
	// solve for the rotational energy levels
	if(system->quantum_rotation) quantum_system_rotational_energies(system);
#endif // QM_ROTATION

	// if root, open necessary output files 
	if(!rank) {
		print_observables(system);
		if(open_files(system) < 0) {
			error("MC: could not open files\n");
			return(-1);
		}
	} 

	write_observables(system->file_pointers.fp_energy, system, system->observables);
	if ( system->file_pointers.fp_energy_csv ) {
		write_observables_csv(system->file_pointers.fp_energy_csv, system, system->observables);
	}

	// close any open files 
	if(!rank)
		close_files(system);

	return(0);
}
