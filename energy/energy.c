/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

/* count the number of molecules currently in the system excluding frozen, adiabatic, etc.*/
void countN ( system_t * system ) {
	molecule_t * molecule_ptr;

	system->observables->N = 0;
	system->observables->spin_ratio = 0;
	for(molecule_ptr = system->molecules, system->observables->N = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {

		if(!(molecule_ptr->frozen || molecule_ptr->adiabatic || molecule_ptr->target)) {
			/* update the molecule counter */
			system->observables->N += 1.0;
			/* update the nuclear spin ratio */
			if(molecule_ptr->nuclear_spin == NUCLEAR_SPIN_ORTHO)
				system->observables->spin_ratio += 1.0;
		}

		if(system->ensemble == ENSEMBLE_NVE) system->N = system->observables->N;
	}

	return;
}


/*check cavity_autoreject_absolute 
-- probably not the most efficient place to put this, but likely the safest*/
double cavity_absolute_check ( system_t * system ) {
	molecule_t * molecule_ptr;
	atom_t * atom_ptr;
	pair_t * pair_ptr;
	
	for ( molecule_ptr=system->molecules; molecule_ptr; molecule_ptr=molecule_ptr->next ) {
		for ( atom_ptr=molecule_ptr->atoms; atom_ptr; atom_ptr=atom_ptr->next ) {
			for ( pair_ptr=atom_ptr->pairs; pair_ptr; pair_ptr=pair_ptr->next ) {
				if ( molecule_ptr == pair_ptr->molecule ) continue; //skip if on the same molecule
				if ( pair_ptr->rimg < system->cavity_autoreject_scale ) return MAXVALUE;
			}
		}
	}
	return 0;
}

int getnatoms (system_t * system) {
	molecule_t * mptr;
	atom_t * aptr;
	int n=0;

	for (mptr = system->molecules; mptr; mptr=mptr->next ) {
		for (aptr = mptr->atoms; aptr; aptr=aptr->next ) {
			n++;
		}
	}
	return n;
}

/* returns the total potential energy for the system and updates our observables */
double energy(system_t *system) {

	molecule_t *molecule_ptr;
	double potential_energy, rd_energy, coulombic_energy, polar_energy, vdw_energy;
	double kinetic_energy;
	struct timeval old_time, new_time;
	char linebuf[MAXLINE];

	/* zero the initial values */
	potential_energy = 0;
	rd_energy = 0;
	coulombic_energy = 0;
	polar_energy = 0;
	vdw_energy = 0;

	/* get the periodic boundary conditions */
	if(system->wpi || system->fvm) pbc(system);	/* do this each time only for fvm and wpi */

	/* get the pairwise terms necessary for the energy calculation */
	pairs(system);

	/* only on the first simulation step, make sure that all recalculate flags are set */
	if(system->observables->energy == 0.0) flag_all_pairs(system);

	/* if we made a volume change (or just reverted from one) set recalculate flags OR if replaying a trajectory */
	// we set last_volume at the end of this function
	if ( system->last_volume != system->pbc->volume || system->ensemble==ENSEMBLE_REPLAY ) 
		flag_all_pairs(system);

	/* get the repulsion/dispersion potential */
	if(system->rd_anharmonic)
		rd_energy = anharmonic(system);
	else if(system->sg)
		rd_energy = sg(system);
	else if(system->dreiding)
		rd_energy = dreiding(system);
	else if(!system->gwp)
		rd_energy = lj(system);
	system->observables->rd_energy = rd_energy;

	/* get the electrostatic potential */
	if(!(system->sg || system->rd_only)) {

		if(system->spectre)
			coulombic_energy = coulombic_nopbc(system->molecules);
		else if(system->gwp) {
			coulombic_energy = coulombic_nopbc_gwp(system);
			kinetic_energy = coulombic_kinetic_gwp(system);
			system->observables->kinetic_energy = kinetic_energy;
		} else
			coulombic_energy = coulombic(system);
		system->observables->coulombic_energy = coulombic_energy;

		/* get the polarization potential */
		if(system->polarization) {

			/* get timing of polarization energy function for cuda comparison */
			gettimeofday(&old_time, NULL);

#ifdef CUDA
			if(system->cuda)
				polar_energy = (double)polar_cuda(system);
			else
				polar_energy = polar(system);
#else
			polar_energy = polar(system);
#endif /* CUDA */
			gettimeofday(&new_time, NULL);
			sprintf(linebuf, "OUTPUT: Polarization energy function took %ld us\n", (new_time.tv_sec-old_time.tv_sec)*1000000+(new_time.tv_usec-old_time.tv_usec));
			/* XXX uncomment for timing */
//			output(linebuf);

			system->observables->polarization_energy = polar_energy;

		}
		if (system->polarvdw) {
#ifdef CUDA
			if (system->cuda) {
				error("error: cuda polarvdw not yet implemented!\n");
				die(-1);
			}
			else
			vdw_energy = vdw(system);
#else
			vdw_energy = vdw(system);
#endif
			system->observables->vdw_energy = vdw_energy;
		}

	}

	/* sum the total potential energy */
	potential_energy = rd_energy + coulombic_energy + polar_energy + vdw_energy;
	/* not truly potential, but stick it there for convenience of MC */
	if(system->gwp) potential_energy += kinetic_energy;
	system->observables->energy = potential_energy;

	countN(system);
	system->observables->spin_ratio /= system->observables->N;

	/* for NVE */
	if(system->ensemble == ENSEMBLE_NVE) {
		system->observables->kinetic_energy = system->total_energy - potential_energy;
		system->observables->temperature = (2.0/3.0)*system->observables->kinetic_energy/system->observables->N;
	}

	/* SA */
	if(system->simulated_annealing)
		system->observables->temperature = system->temperature;

	/* need this for the isosteric heat */
	system->observables->NU = system->observables->N*system->observables->energy;

	/* set last known volume*/
	system->last_volume = system->pbc->volume;

	if(system->cavity_autoreject_absolute)
		potential_energy += cavity_absolute_check( system );

	return(potential_energy);

}



/// Not adjusted to work with NPT. need to run recalculates whenever volume changes

/* returns the total potential energy for the system */
/* this function is meant to be called by routines that do not */
/* require observables to be averaged in, i.e. quantum integration */
/* routines, widom insertion, etc. */
double energy_no_observables(system_t *system) {

	double potential_energy, rd_energy, coulombic_energy, polar_energy, vdw_energy;

	/* zero the initial values */
	potential_energy = 0;
	rd_energy = 0;
	coulombic_energy = 0;
	polar_energy = 0;
	vdw_energy = 0;

	/* get the periodic boundary conditions */
	if(system->wpi || system->fvm) pbc(system);	/* do this each time only for fvm and wpi */

	/* get the pairwise terms necessary for the energy calculation */
	pairs(system);

	/* get the repulsion/dispersion potential */
	if(system->sg)
		rd_energy = sg(system);
	else if(system->dreiding)
		rd_energy = dreiding(system);
	else
		rd_energy = lj(system);

	/* get the electrostatic potential */
	if(!(system->sg || system->rd_only)) {

		coulombic_energy = coulombic(system);

		/* get the polarization potential */
		if(system->polarization)
			polar_energy = polar(system);

		if(system->polarvdw)
			vdw_energy = vdw(system);

	}



	/* sum the total potential energy */
	potential_energy = rd_energy + coulombic_energy + polar_energy + vdw_energy;


	return(potential_energy);


}

