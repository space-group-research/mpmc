/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>


/* get the induction energy */
double polar(system_t *system) {

	int i, num_iterations;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	double dipole_rrms, N, potential;


	/* take measures to let N fluctuate */
	if((system->ensemble == ENSEMBLE_UVT) && !system->polar_zodid)
		thole_resize_matrices(system);

	/* get the A matrix */
	if(!system->polar_zodid) {

		thole_amatrix(system);
		if(system->polarizability_tensor) {

			output("POLAR: A matrix:\n");
			print_matrix(3*((int)system->checkpoint->N_atom), system->A_matrix);

		}

	}

	/* calculate the field vectors */
	thole_field(system);


	/* find the dipoles */
	if(system->polar_iterative) {	/* solve the self-consistent field... */

		num_iterations = thole_iterative(system);
		system->nodestats->polarization_iterations = (double)num_iterations;

		/* RRMS of dipoles */
		N = 0; dipole_rrms = 0;
		for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
			for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
				if(isfinite(atom_ptr->dipole_rrms)) dipole_rrms += atom_ptr->dipole_rrms;
				N += 1.0;
			}
		}
		dipole_rrms /= N;
		system->observables->dipole_rrms = dipole_rrms;

	} else {	/* ...or do matrix inversion */

		thole_bmatrix(system);
		thole_bmatrix_dipoles(system);

		/* output the 3x3 molecular polarizability tensor */
		if(system->polarizability_tensor) {

			output("POLAR: B matrix:\n");
			print_matrix(3*((int)system->checkpoint->N_atom), system->B_matrix);

			thole_polarizability_tensor(system);
			exit(0);

		}

	}


	/* calculate the polarization energy as 1/2 mu*E */
	for(molecule_ptr = system->molecules, potential = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(i = 0; i < 3; i++) {
				potential += atom_ptr->mu[i]*atom_ptr->ef_static[i];
				if(system->polar_palmo)
					potential += atom_ptr->mu[i]*atom_ptr->ef_induced_change[i];
			}
		}

	}
	potential *= -0.5;

	return(potential);

}


