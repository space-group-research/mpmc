/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>


/* Lennard-Jones repulsion/dispersion */
double lj(system_t *system) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	pair_t *pair_ptr;
	double sigma_over_r, term12, term6;
	double first_derivative, second_derivative, third_derivative, fourth_derivative;
	double potential, potential_classical, potential_fh_second_order, potential_fh_fourth_order;
	double reduced_mass;
	double sig3, sig_cut, sig_cut3, sig_cut9;

	potential = 0;
	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if(pair_ptr->recalculate_energy) {

					pair_ptr->rd_energy = 0;

					/* make sure we're not excluded or beyond the cutoff */
					if(!((pair_ptr->rimg > system->pbc->cutoff) || pair_ptr->rd_excluded || pair_ptr->frozen)) {

						sigma_over_r = fabs(pair_ptr->sigma)/pair_ptr->rimg;

						/* the LJ potential */
						if(system->spectre) {
							term6 = 0;
							term12 = system->scale_rd*pow(sigma_over_r, 12.0);
							potential_classical = term12;
						} else {
							term6 = pow(sigma_over_r, 6.0);
							term6 *= system->scale_rd;

							if(pair_ptr->attractive_only)
								term12 = 0;
							else
								term12 = pow(term6, 2.0);

							if (system->polarvdw) //vdw calculated by vdw.c
								term6=0;

							potential_classical = 4.0*pair_ptr->epsilon*(term12 - term6);
						}

						pair_ptr->rd_energy += potential_classical;

						if(system->feynman_hibbs) {

							reduced_mass = AMU2KG*molecule_ptr->mass*pair_ptr->molecule->mass/(molecule_ptr->mass+pair_ptr->molecule->mass);

							/* FIRST DERIVATIVE */
							first_derivative = -24.0*pair_ptr->epsilon*(2.0*term12 - term6)/pair_ptr->rimg;

							/* SECOND DERIVATIVE */
							second_derivative = 24.0*pair_ptr->epsilon*(26.0*term12 - 7.0*term6)/pow(pair_ptr->rimg, 2.0);

							potential_fh_second_order = pow(METER2ANGSTROM, 2.0)*(HBAR*HBAR/(24.0*KB*system->temperature*reduced_mass))*(second_derivative + 2.0*first_derivative/pair_ptr->rimg);
							pair_ptr->rd_energy += potential_fh_second_order;

							if(system->feynman_hibbs_order >= 4) {

								/* THIRD DERIVATIVE */
								third_derivative = -1344.0*pair_ptr->epsilon*(6.0*term12 - term6)/pow(pair_ptr->rimg, 3.0);

								/* FOURTH DERIVATIVE */
								fourth_derivative = 12096.0*pair_ptr->epsilon*(10.0*term12 - term6)/pow(pair_ptr->rimg, 4.0);

								potential_fh_fourth_order = pow(METER2ANGSTROM, 4.0)*(pow(HBAR, 4.0)/(1152.0*pow(KB*system->temperature*reduced_mass, 2.0)))*(15.0*first_derivative/pow(pair_ptr->rimg, 3.0) + 4.0*third_derivative/pair_ptr->rimg + fourth_derivative);
								pair_ptr->rd_energy += potential_fh_fourth_order;

							}

						}

						/* cause an autoreject on insertions */ //cavity_autoreject_absolute is checked in energy.c
						if(system->cavity_autoreject) {
							if(pair_ptr->rimg < system->cavity_autoreject_scale*fabs(pair_ptr->sigma))
								pair_ptr->rd_energy = MAXVALUE;
						}

					}

					/* include the long-range correction */  /* I'm  not sure that I'm handling spectre pairs correctly */
					/* we can't use rd_excluded flag, since that disqualifies inter-molecular, but that DOES contribute to LRC */
					/* ALL OF THESE MUST BE TRUE TO PERFORM LRC CALCULATION */
					if( ( pair_ptr->epsilon != 0 && pair_ptr->sigma != 0 ) &&  //if these are zero, then we won't waste our time
							!( atom_ptr->spectre && pair_ptr->atom->spectre ) && //i think we want to disqualify s-s pairs 
							!( pair_ptr->frozen ) &&  //disqualify frozen pairs
							((pair_ptr->lrc == 0.0) || system->last_volume != system->pbc->volume) &&  //LRC only changes if the volume changes
							system->rd_lrc) { //flag to include LRC

						sig_cut = fabs(pair_ptr->sigma)/system->pbc->cutoff;
						sig3 = pow(fabs(pair_ptr->sigma), 3.0);
						sig_cut3 = pow(sig_cut, 3.0);
						sig_cut9 = pow(sig_cut, 9.0);

            if ( system->polarvdw ) 
              //only repulsion term, if polarvdw is on
              pair_ptr->lrc = (16.0/9.0)*M_PI*pair_ptr->epsilon*sig3*sig_cut9/system->pbc->volume;
            else 
              //if polarvdw is off, do the usual thing
              pair_ptr->lrc = ((16.0/3.0)*M_PI*pair_ptr->epsilon*sig3)*((1.0/3.0)*sig_cut9 - sig_cut3)/system->pbc->volume;
		
					}

				} /* if recalculate */

				/* sum all of the pairwise terms */
				potential += pair_ptr->rd_energy + pair_ptr->lrc;

			} /* pair */
		} /* atom */
	} /* molecule */

	/* calculate self LRC interaction */

	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

			if ( ((atom_ptr->sigma != 0)  && (atom_ptr->epsilon != 0)) && //non-zero parameters
					 !(atom_ptr->frozen) && //not frozen
					 !(atom_ptr->spectre) && //not spectre
					 (system->rd_lrc) ) { // flag for LRC calculation
		
				sig_cut = fabs(atom_ptr->sigma)/system->pbc->cutoff;
				sig3 = pow(fabs(atom_ptr->sigma), 3.0);
				sig_cut3 = pow(sig_cut, 3.0);
				sig_cut9 = pow(sig_cut, 9.0);

				if ( system->polarvdw ) 
				//only repulsion term, if polarvdw is on
					potential += (16.0/9.0)*M_PI*atom_ptr->epsilon*sig3*sig_cut9/system->pbc->volume;
				else 
				//if polarvdw is off, do the usual thing
					potential += ((16.0/3.0)*M_PI*atom_ptr->epsilon*sig3)*((1.0/3.0)*sig_cut9 - sig_cut3)/system->pbc->volume;
			}

		}
	}


	return(potential);

}



/* same as above, but no periodic boundary conditions */
double lj_nopbc(system_t * system) {

	molecule_t *molecules = system->molecules;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	pair_t *pair_ptr;
	double sigma_over_r, term12, term6;
	double potential;

	for(molecule_ptr = molecules, potential = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				/* make sure we're not excluded or beyond the cutoff */
				if(!pair_ptr->rd_excluded) {

					sigma_over_r = fabs(pair_ptr->sigma)/pair_ptr->r;

					/* the LJ potential */
					term6 = pow(sigma_over_r, 6.0);
					if(pair_ptr->attractive_only)
						term12 = 0;
					else
						term12 = pow(term6, 2.0);

					if ( system->polarvdw )  term6=0;
					potential += 4.0*pair_ptr->epsilon*(term12 - term6);

				}

			} /* pair */
		} /* atom */
	} /* molecule */


	return(potential);

}


