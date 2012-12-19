/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

double lj_fh_corr( system_t * system, molecule_t * molecule_ptr, pair_t * pair_ptr, int order, double term12, double term6 ) {
	double reduced_mass;
	double dE, d2E, d3E, d4E; //energy derivatives
	double corr;
	double ir = 1.0/pair_ptr->rimg;
	double ir2 = ir*ir;
	double ir3 = ir2*ir;
	double ir4 = ir3*ir;

	if ( (order != 2) && (order != 4) ) return NAN; //must be order 2 or 4

	reduced_mass = AMU2KG*molecule_ptr->mass*pair_ptr->molecule->mass /
		(molecule_ptr->mass+pair_ptr->molecule->mass);

	dE = -24.0*pair_ptr->epsilon*(2.0*term12 - term6) * ir;
	d2E = 24.0*pair_ptr->epsilon*(26.0*term12 - 7.0*term6) * ir2;

	//2nd order correction
	corr = M2A2 *
		(HBAR2/(24.0*KB*system->temperature*reduced_mass)) *
		(d2E + 2.0*dE/pair_ptr->rimg);

	if(order >= 4) {

		d3E = -1344.0*pair_ptr->epsilon*(6.0*term12 - term6) * ir3;
		d4E = 12096.0*pair_ptr->epsilon*(10.0*term12 - term6) * ir4;
	
		//4th order corection
		corr += M2A4 *
			(HBAR4/(1152.0*KB2*system->temperature*system->temperature*reduced_mass*reduced_mass)) *
			( 15.0*dE*ir3 +	4.0*d3E*ir + d4E );
	}

	return corr;
}

double lj_lrc_corr( system_t * system, atom_t * atom_ptr,  pair_t * pair_ptr, double cutoff ) {

	double sig_cut, sig3, sig_cut3, sig_cut9;
	double corr;

	/* include the long-range correction */  /* I'm  not sure that I'm handling spectre pairs correctly */
	/* we can't use rd_excluded flag, since that disqualifies inter-molecular, but that DOES contribute to LRC */
	/* ALL OF THESE MUST BE TRUE TO PERFORM LRC CALCULATION */
	if( ( pair_ptr->epsilon != 0 && pair_ptr->sigma != 0 ) &&  //if these are zero, then we won't waste our time
			!( atom_ptr->spectre && pair_ptr->atom->spectre ) && //i think we want to disqualify s-s pairs 
			!( pair_ptr->frozen ) &&  //disqualify frozen pairs
			((pair_ptr->lrc == 0.0) || pair_ptr->last_volume != system->pbc->volume) ) { //LRC only changes if the volume change

		pair_ptr->last_volume = system->pbc->volume;

		sig_cut = fabs(pair_ptr->sigma)/cutoff;
		sig3 = fabs(pair_ptr->sigma);
		sig3 *= sig3*sig3;
		sig_cut3 = sig_cut*sig_cut*sig_cut;
		sig_cut9 = sig_cut3*sig_cut3*sig_cut3;

		if ( system->polarvdw ) //only repulsion term, if polarvdw is on
			corr = (16.0/9.0)*M_PI*pair_ptr->epsilon*sig3*sig_cut9/system->pbc->volume;
		else //if polarvdw is off, do the usual thing
			corr = ((16.0/3.0)*M_PI*pair_ptr->epsilon*sig3)*((1.0/3.0)*sig_cut9 - sig_cut3)/system->pbc->volume;
	}
	else corr = pair_ptr->lrc; //use stored value

	return corr;
}

double lj_lrc_self ( system_t * system, atom_t * atom_ptr, double cutoff ) {
	double sig_cut, sig3, sig_cut3, sig_cut9;
	double corr = 0;

	if ( ((atom_ptr->sigma != 0)  && (atom_ptr->epsilon != 0)) && //non-zero parameters
		 !(atom_ptr->frozen) && //not frozen
		 !(atom_ptr->spectre) ) { //not spectre 
		
		sig_cut = fabs(atom_ptr->sigma)/cutoff;
		sig3 = fabs(atom_ptr->sigma);
		sig3 *= sig3*sig3;
		sig_cut3 = sig_cut*sig_cut*sig_cut;
		sig_cut9 = sig_cut3*sig_cut3*sig_cut3;

		if ( system->polarvdw ) //only repulsion term, if polarvdw is on
			corr = (16.0/9.0)*M_PI*atom_ptr->epsilon*sig3*sig_cut9/system->pbc->volume;
		else //if polarvdw is off, do the usual thing
			corr = ((16.0/3.0)*M_PI*atom_ptr->epsilon*sig3)*((1.0/3.0)*sig_cut9 - sig_cut3)/system->pbc->volume;
	}

	return corr;
}


/* Lennard-Jones repulsion/dispersion */
double lj(system_t *system) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	pair_t *pair_ptr;
	double sigma_over_r, term12, term6, sigma_over_r6, sigma_over_r12;
	double potential, potential_classical, cutoff;
	int i,j,k;

	//set the cutoff
	if ( system->rd_crystal )
		cutoff = 2.0 * system->pbc->cutoff * ((double)system->rd_crystal_order - 0.5);
	else
		cutoff = system->pbc->cutoff;

	potential = 0;
	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if(pair_ptr->recalculate_energy) {

					pair_ptr->rd_energy = 0;

					// pair LRC (lj_crystal option should check that rd_lrc is off)
					if ( system->rd_lrc ) pair_ptr->lrc = lj_lrc_corr(system,atom_ptr,pair_ptr,cutoff);

					// to include a contribution, we require
					if ( 	( pair_ptr->rimg - SMALL_dR < cutoff ) && //inside cutoff?
								! pair_ptr->rd_excluded	 && //not excluded
								! pair_ptr->frozen ) { //not frozen

						//loop over unit cells
						if ( system->rd_crystal ) {
							sigma_over_r6 = 0;
							sigma_over_r12 = 0;
							for ( i = -system->rd_crystal_order+1; i<=system->rd_crystal_order-1; i++ )
							for ( j = -system->rd_crystal_order+1; j<=system->rd_crystal_order-1; j++ )
							for ( k = -system->rd_crystal_order+1; k<=system->rd_crystal_order-1; k++ ) {
								sigma_over_r = fabs(pair_ptr->sigma)/sqrt(
									pow(atom_ptr->pos[0] - pair_ptr->atom->pos[0] + i*system->pbc->basis[0][0],2) +
									pow(atom_ptr->pos[1] - pair_ptr->atom->pos[1] + j*system->pbc->basis[1][1],2) +
									pow(atom_ptr->pos[2] - pair_ptr->atom->pos[2] + k*system->pbc->basis[2][2],2) );
								sigma_over_r6 += pow(sigma_over_r,6);
								sigma_over_r12 += pow(sigma_over_r,12);
							}
						}
						else { //otherwise, calculate as normal
							sigma_over_r = fabs(pair_ptr->sigma)/pair_ptr->rimg;
							sigma_over_r6 = sigma_over_r*sigma_over_r*sigma_over_r;
							sigma_over_r6 *= sigma_over_r6;
							sigma_over_r12 = sigma_over_r6 * sigma_over_r6;
						}

						/* the LJ potential */
						if(system->spectre) {
							term6 = 0;
							term12 = system->scale_rd*sigma_over_r12;
							potential_classical = term12;
						} 
						else {
							if ( system->polarvdw ) term6=0; //vdw calc'd by vdw.c
								else {
									term6 = sigma_over_r6;
									term6 *= system->scale_rd;
								}

							if(pair_ptr->attractive_only) term12 = 0;
								else term12 = sigma_over_r12;

							potential_classical = 4.0*pair_ptr->epsilon*(term12 - term6);
						}

						pair_ptr->rd_energy += potential_classical;

						if(system->feynman_hibbs) 
							pair_ptr->rd_energy += lj_fh_corr(system,molecule_ptr,pair_ptr,system->feynman_hibbs_order, term12, term6);

						// if cavity_autoreject is on (cavity_autoreject_absolute is performed in energy.c)
						if(system->cavity_autoreject)
							if(pair_ptr->rimg < system->cavity_autoreject_scale*fabs(pair_ptr->sigma))
								pair_ptr->rd_energy = MAXVALUE;

					}


				} /* if recalculate */

				/* sum all of the pairwise terms */
				potential += pair_ptr->rd_energy + pair_ptr->lrc;

			} /* pair */
		} /* atom */
	} /* molecule */

	/* calculate self LRC interaction */
	if ( system->rd_lrc ) 
		for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
			for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
				potential += lj_lrc_self(system,atom_ptr,cutoff);

	return(potential);

}


/* same as above, but no periodic boundary conditions */
double lj_nopbc(system_t * system) {

	molecule_t *molecules = system->molecules;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	pair_t *pair_ptr;
	double sigma_over_r, term12, term6, sigma_over_r6;
	double potential;

	for(molecule_ptr = molecules, potential = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				/* make sure we're not excluded or beyond the cutoff */
				if(!pair_ptr->rd_excluded) {

					sigma_over_r = fabs(pair_ptr->sigma)/pair_ptr->r;
					sigma_over_r6 = sigma_over_r*sigma_over_r*sigma_over_r;
					sigma_over_r6 *= sigma_over_r6;

					if ( system->polarvdw ) term6=0;
						else term6 = sigma_over_r6;

					if(pair_ptr->attractive_only) term12 = 0;
						else term12 = sigma_over_r6*sigma_over_r6;

					potential += 4.0*pair_ptr->epsilon*(term12 - term6);

				}

			} /* pair */
		} /* atom */
	} /* molecule */


	return(potential);

}


