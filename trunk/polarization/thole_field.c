/* 

Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

//called from energy/polar.c
/* calculate the field with periodic boundaries */
void thole_field(system_t *system) {

	int p;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;

	/* zero the field vectors */
	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

			for(p = 0; p < 3; p++) {
				atom_ptr->ef_static[p] = 0;
				atom_ptr->ef_static_self[p] = 0;
			}

		}
	}

	/* calculate the electrostatic field */
	if(system->polar_ewald)
		ewald_estatic(system);
	else if (system->polar_wolf || system->polar_wolf_full)
		thole_field_wolf(system);
	else
		thole_field_nopbc(system);

}

/* calculate the field without ewald summation/wolf */
void thole_field_nopbc(system_t *system) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr, **atom_array;
	pair_t *pair_ptr;
	int p;
	double r;

	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if(!pair_ptr->frozen) {

					r = pair_ptr->rimg;

					//inclusive near the cutoff
					if((r - SMALL_dR < system->pbc->cutoff) && (r != 0.)) {

						if(pair_ptr->es_excluded) {
							if(system->polar_self) { //self-induction
								for(p = 0; p < 3; p++) {
									atom_ptr->ef_static_self[p] += pair_ptr->atom->charge*pair_ptr->dimg[p]/(r*r*r);
									pair_ptr->atom->ef_static_self[p] -= atom_ptr->charge*pair_ptr->dimg[p]/(r*r*r);
								}
							}
						} else {
							for(p = 0; p < 3; p++) {
								atom_ptr->ef_static[p] += pair_ptr->atom->charge*pair_ptr->dimg[p]/(r*r*r);
								pair_ptr->atom->ef_static[p] -= atom_ptr->charge*pair_ptr->dimg[p]/(r*r*r);
							}
						}

					} /* cutoff */

				} /* frozen */

			} /* pair */
		} /* atom */
	} /* molecule */

	return;
}

// calc field using wolf sum (JCP 124 234104 (2006) equation 19
void thole_field_wolf(system_t *system) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr, **atom_array;
	pair_t *pair_ptr;
	int p; //dimensionality
	double r, rr; //r and 1/r (reciprocal of r)
	double R = system->pbc->cutoff;
	double rR = 1./R;
	//used for polar_wolf_alpha (aka polar_wolf_damp)
	double a = system->polar_wolf_alpha;
	double err, erR; //complementary error functions
	erR=erfc(a*R);

	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if ( molecule_ptr == pair_ptr->molecule ) continue; //don't let molecules polarize themselves
				if ( pair_ptr->frozen ) continue; //don't let the MOF polarize itself

				r = pair_ptr->rimg;
				rr = 1./r;

				if((r - SMALL_dR < system->pbc->cutoff) && (r != 0.)) {
					for ( p=0; p<3; p++ ) { 
						//see JCP 124 (234104)
						if ( a == 0 ) {
							atom_ptr->ef_static[p] += (pair_ptr->atom->charge)*(rr*rr-rR*rR)*pair_ptr->dimg[p]*rr;
							pair_ptr->atom->ef_static[p] -= (atom_ptr->charge)*(rr*rr-rR*rR)*pair_ptr->dimg[p]*rr;
						} else {
							err=erfc(a*r);
							atom_ptr->ef_static[p] += pair_ptr->atom->charge*((err*rr*rr+2.0*a/sqrt(M_PI)*exp(-a*a*r*r)*rr) -
								(erR*rR*rR + 2.0*a/sqrt(M_PI)*exp(-a*a*R*R)*rR))*pair_ptr->dimg[p]*rr;
							pair_ptr->atom->ef_static[p] -= atom_ptr->charge*((err*rr*rr+2.0*a/sqrt(M_PI)*exp(-a*a*r*r)*rr) -
								(erR*rR*rR + 2.0*a/sqrt(M_PI)*exp(-a*a*R*R)*rR))*pair_ptr->dimg[p]*rr;
						}
						
					}
				}

			} /* pair */
		} /* atom */
	} /* molecule */

	return;
}

