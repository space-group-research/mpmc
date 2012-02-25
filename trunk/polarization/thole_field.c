/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>


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
	if(system->polar_ewald) {
		thole_field_real(system);
		thole_field_recip(system);
		thole_field_self(system);
	} else
		thole_field_nopbc(system);


}

/* calculate the field without ewald summation */
void thole_field_nopbc(system_t *system) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr, **atom_array;
	pair_t *pair_ptr;
	int p;
	double r, damping;

	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if(!pair_ptr->frozen) {

					r = pair_ptr->rimg;
					/* exponential field damping */
					if(system->field_damp == 0.0)
						damping = 1.0;
					else
						damping = 1.0 - exp(-pow(r/system->field_damp, 3.0));

					/* include self-induction if keyword polar_self is set */
					if((r < system->pbc->cutoff) && (r != 0.)) {

						if(pair_ptr->es_excluded) {

							if(system->polar_self) {

								for(p = 0; p < 3; p++) {
									atom_ptr->ef_static_self[p] += damping*pair_ptr->atom->charge*pair_ptr->dimg[p]/(r*r*r);
									pair_ptr->atom->ef_static_self[p] -= damping*atom_ptr->charge*pair_ptr->dimg[p]/(r*r*r);
								}


							}

						} else {

							for(p = 0; p < 3; p++) {
								atom_ptr->ef_static[p] += damping*pair_ptr->atom->charge*pair_ptr->dimg[p]/(r*r*r);
								pair_ptr->atom->ef_static[p] -= damping*atom_ptr->charge*pair_ptr->dimg[p]/(r*r*r);
							}

						}

					} /* cutoff */


				} /* frozen */

			} /* pair */
		} /* atom */
	} /* molecule */


}


/* real space sum */
void thole_field_real(system_t *system) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	pair_t *pair_ptr;
	int p;
	double alpha;
	double field_value;

	alpha = system->ewald_alpha;

	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if(!pair_ptr->frozen) {

					if(!pair_ptr->es_excluded) {

						field_value = erfc(alpha*pair_ptr->rimg) + 2.0*alpha*pair_ptr->rimg*exp(-alpha*alpha*pair_ptr->rimg*pair_ptr->rimg)/sqrt(M_PI);

						for(p = 0; p < 3; p++) {
							atom_ptr->ef_static[p] += pair_ptr->atom->charge*field_value*pair_ptr->dimg[p]/(pair_ptr->rimg*pair_ptr->rimg*pair_ptr->rimg);
							pair_ptr->atom->ef_static[p] -= atom_ptr->charge*field_value*pair_ptr->dimg[p]/(pair_ptr->rimg*pair_ptr->rimg*pair_ptr->rimg);
						}

					}

				} /* !frozen */

			} /* pair */
		} /* atom */
	} /* molecule */


}


/* fourier space sum */
void thole_field_recip(system_t *system) {

	int i, N, r, p, q;
	int kmax;
	double alpha, gaussian;
	double vector_product_i, vector_product_j;
	double SF_sin, SF_cos;
	int norm, l[3];
	int k_squared, k[3];
	molecule_t *molecule_ptr;
	atom_t *atom_ptr, **atom_array;
	pair_t *pair_ptr;

	/* our convergence parameters */
	alpha = system->ewald_alpha;
	kmax = system->ewald_kmax;

	/* generate an array of atom ptrs */
	for(molecule_ptr = system->molecules, N = 0, atom_array = NULL; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next, N++) {

			atom_array = realloc(atom_array, sizeof(atom_t *)*(N + 1));
			memnullcheck(atom_array,sizeof(atom_t *)*(N+1),58);
			atom_array[N] = atom_ptr;

		}
	}

	/* calculate the field vector for each nuclear coordinate */
	for(i = 0; i < N; i++) {

		/* perform the fourier sum over the reciprocal lattice */
		for(l[0] = 0; l[0] <= kmax; l[0]++) {
			for(l[1] = (!l[0] ? 0 : -kmax); l[1] <= kmax; l[1]++) {
				for(l[2] = ((!l[0] && !l[1]) ? 1 : -kmax); l[2] <= kmax; l[2]++) {

					/* k-vector norm */
					for(p = 0, norm = 0; p < 3; p++)
						norm += l[p]*l[p];

					/* get the reciprocal lattice vectors */
					for(p = 0, k_squared = 0; p < 3; p++) {
						for(q = 0, k[p] = 0; q < 3; q++)
							k[p] += system->pbc->reciprocal_basis[p][q]*2.0*M_PI*((double)l[q]);
						k_squared += k[p]*k[p];
					}

					/* compare the norm */
					if((norm <= kmax*kmax) && (k_squared > 0.0)) {

						/* our gaussian centered on the charge */
						gaussian = exp(-k_squared/(4.0*alpha*alpha))/k_squared;

						/* the structure factor */
						SF_sin = 0; SF_cos = 0;
						for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
							for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

								/* make sure it's not a frozen-frozen interaction */
								if(!(atom_ptr->frozen && atom_array[i]->frozen)) {

									/* inner product of j-th position vector and k-vector */
									for(p = 0, vector_product_j = 0; p < 3; p++)
										vector_product_j += k[p]*atom_ptr->pos[p];

									SF_sin += atom_ptr->charge*sin(vector_product_j);
									SF_cos += -atom_ptr->charge*cos(vector_product_j);

								} /* !frozen */

							} /* atom */
						} /* molecule */

						/* inner product of i-th position vector and k-vector */
						for(p = 0, vector_product_i = 0; p < 3; p++)
							vector_product_i += k[p]*atom_array[i]->pos[p];

						/* the i-th terms outside of the SF sum */
						SF_sin *= cos(vector_product_i);
						SF_cos *= sin(vector_product_i);

						/* project the field components */
						for(r = 0; r < 3; r++ )
							atom_array[i]->ef_static[r] += -4.0*M_PI*k[r]*gaussian*(SF_sin + SF_cos)/system->pbc->volume;


					} /* end if norm */

				} /* end for n */
			} /* end for m */
		} /* end for l */

	} /* end i */

	free(atom_array);

}


/* self interaction term */
void thole_field_self(system_t *system) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	pair_t *pair_ptr;
	int p;
	double alpha;
	double field_value;
	double d[3];

	alpha = system->ewald_alpha;

	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if(!pair_ptr->frozen) {

					if(pair_ptr->es_excluded) {

						field_value = erf(alpha*pair_ptr->r) - 2.0*alpha*pair_ptr->r*exp(-alpha*alpha*pair_ptr->r*pair_ptr->r)/sqrt(M_PI);

						for(p = 0; p < 3; p++) {
							d[p] = atom_ptr->pos[p] - pair_ptr->atom->pos[p];
							atom_ptr->ef_static[p] -= pair_ptr->atom->charge*field_value*d[p]/(pair_ptr->r*pair_ptr->r*pair_ptr->r);
							pair_ptr->atom->ef_static[p] += atom_ptr->charge*field_value*d[p]/(pair_ptr->r*pair_ptr->r*pair_ptr->r);
						}

					}

				} /* !frozen */

			} /* pair */
		} /* atom */
	} /* molecule */


}


