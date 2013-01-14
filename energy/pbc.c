/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

/* calculate the minimum cutoff radius from the basis lattice */
/* generalized to handle skewed unit cells */
/* there is most likely a more elegant way to find this */
double pbc_cutoff(pbc_t *pbc) {

	int p, q;
	int v1, v2;				/* any two basis vectors */
	double v1_magnitude;			/* magnitude of the first basis */
	double component;			/* the terms a*a + b*a */
	double r_vector[3], r_magnitude;	/* the new vector */
	double rmin[3];				/* the radial cutoff for each basis pair */
	double cutoff;				/* the minimal cutoff */

	/* for each pair of basis vectors */
	for(p = 0; p < 3; p++) {

		/* the unique basis pairs forming the parallelogram */
		v1 = p;
		v2 = (p + 1) % 3;

		/* calculate the first basis magnitude */
		v1_magnitude = 0;
		for(q = 0; q < 3; q++) {
			v1_magnitude += pbc->basis[v1][q]*pbc->basis[v1][q];
		}
		v1_magnitude = sqrt(v1_magnitude);

		/* compute vector components of aa + ba */
		component = 0;
		for(q = 0; q < 3; q++) {
			component += pbc->basis[v1][q]*pbc->basis[v1][q];
			component += pbc->basis[v2][q]*pbc->basis[v1][q];
		}
		component /= v1_magnitude;

		/* finally, the r vector itself */
		/* r = 1/2[a + b - (aa + ba)a] */
		r_magnitude = 0;
		for(q = 0; q < 3; q++) {
			r_vector[q] = pbc->basis[v1][q] + pbc->basis[v2][q];
			r_vector[q] -= component*pbc->basis[v1][q]/v1_magnitude;
			r_vector[q] *= 0.5;
			r_magnitude += r_vector[q]*r_vector[q];
		}
		r_magnitude = sqrt(r_magnitude);

		/* store the result for this parallelogram - sort it out when done */
		rmin[p] = r_magnitude;

	}

	/* sort out the smallest radial cutoff */
	cutoff = MAXVALUE;
	for(p = 0; p < 3; p++)
		if(rmin[p] < cutoff) cutoff = rmin[p];


	return(cutoff);

}


/* take the determinant of the basis matrix */
double pbc_volume(pbc_t *pbc) {

	double volume;

	volume =  pbc->basis[0][0]*(pbc->basis[1][1]*pbc->basis[2][2] - pbc->basis[1][2]*pbc->basis[2][1]);
	volume += pbc->basis[0][1]*(pbc->basis[1][2]*pbc->basis[2][0] - pbc->basis[1][0]*pbc->basis[2][2]);
	volume += pbc->basis[0][2]*(pbc->basis[1][0]*pbc->basis[2][1] - pbc->basis[2][1]*pbc->basis[2][0]);

	return(volume);
}

/* get the reciprocal space basis */
void pbc_reciprocal(pbc_t *pbc) {

	double inverse_volume;

	inverse_volume = 1.0/pbc_volume(pbc);

	pbc->reciprocal_basis[0][0] = inverse_volume*(pbc->basis[1][1]*pbc->basis[2][2] - pbc->basis[1][2]*pbc->basis[2][1]);
	pbc->reciprocal_basis[0][1] = inverse_volume*(pbc->basis[0][2]*pbc->basis[2][1] - pbc->basis[0][1]*pbc->basis[2][2]);
	pbc->reciprocal_basis[0][2] = inverse_volume*(pbc->basis[0][1]*pbc->basis[1][2] - pbc->basis[0][2]*pbc->basis[1][1]);

	pbc->reciprocal_basis[1][0] = inverse_volume*(pbc->basis[1][2]*pbc->basis[2][0] - pbc->basis[1][0]*pbc->basis[2][2]);
	pbc->reciprocal_basis[1][1] = inverse_volume*(pbc->basis[0][0]*pbc->basis[2][2] - pbc->basis[0][2]*pbc->basis[2][0]);
	pbc->reciprocal_basis[1][2] = inverse_volume*(pbc->basis[0][2]*pbc->basis[1][2] - pbc->basis[0][0]*pbc->basis[1][2]);

	pbc->reciprocal_basis[2][0] = inverse_volume*(pbc->basis[1][0]*pbc->basis[2][1] - pbc->basis[1][1]*pbc->basis[2][0]);
	pbc->reciprocal_basis[2][1] = inverse_volume*(pbc->basis[0][1]*pbc->basis[2][0] - pbc->basis[0][0]*pbc->basis[2][1]);
	pbc->reciprocal_basis[2][2] = inverse_volume*(pbc->basis[0][0]*pbc->basis[1][1] - pbc->basis[0][1]*pbc->basis[1][0]);

}

void pbc(system_t * system) {

	pbc_t * pbc = system->pbc;

	if ( (pbc->cutoff != 0.0) && (system->checkpoint->movetype != MOVETYPE_VOLUME) && system->ensemble != ENSEMBLE_REPLAY )
		return; //nothing to do

	/* get the unit cell volume and cutoff */
	pbc->volume = pbc_volume(pbc);

	pbc->cutoff = pbc_cutoff(pbc);

	// calculate ewald_alpha and polar_ewald_alpha unless manually set
	if ( system->ewald_alpha_set != 1 ) 
		system->ewald_alpha = 3.5/system->pbc->cutoff;
	if ( system->polar_ewald_alpha_set != 1 ) 
		system->polar_ewald_alpha = 3.5/system->pbc->cutoff;

	/* get the reciprocal space lattice */
	pbc_reciprocal(pbc);

	return;
}
