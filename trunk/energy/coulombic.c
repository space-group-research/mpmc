/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

/* total ES energy term */
double coulombic(system_t *system) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	pair_t *pair_ptr;
	double real, reciprocal;
	double potential;

	/* construct the relevant ewald terms */
	real = coulombic_real(system);
	reciprocal = coulombic_reciprocal(system);

	/* return the total electrostatic energy */
	potential = real + reciprocal;

	return(potential);

}


/* fourier space sum */
double coulombic_reciprocal(system_t *system) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	int q, p;
	int kmax;
	double alpha;
	double l[3], k[3], k_squared, norm;
	double gaussian, position_product;
	double SF_real, SF_imaginary;			/* structure factor */
	double recip_potential, self_potential, potential;

	alpha = system->ewald_alpha;
	kmax = system->ewald_kmax;

	recip_potential = 0; self_potential = 0;

	if(system->wolf) {

		for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
			for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

				if(!atom_ptr->frozen) {

					atom_ptr->es_self_point_energy = 0.5*erfc(alpha*system->pbc->cutoff)/system->pbc->cutoff*atom_ptr->charge*atom_ptr->charge;
					atom_ptr->es_self_point_energy += alpha*atom_ptr->charge*atom_ptr->charge/sqrt(M_PI);
					self_potential += atom_ptr->es_self_point_energy;

				} /* ! frozen */

			} /* for atom */
		} /* for molecule */

		potential = -self_potential;

	} else {

		/* perform the fourier sum over the reciprocal lattice for each particle */
		for(l[0] = 0; l[0] <= kmax; l[0]++) {
			for(l[1] = (!l[0] ? 0 : -kmax); l[1] <= kmax; l[1]++) {
				for(l[2] = ((!l[0] && !l[1]) ? 1 : -kmax); l[2] <= kmax; l[2]++) {

					/* compare the norm */
					for(p = 0, norm = 0; p < 3; p++)
						norm += l[p]*l[p];

					/* get the reciprocal lattice vectors */
					for(p = 0, k_squared = 0; p < 3; p++) {
						for(q = 0, k[p] = 0; q < 3; q++)
							k[p] += system->pbc->reciprocal_basis[p][q]*2.0*M_PI*((double)l[q]);
						k_squared += k[p]*k[p];
					}

					/* make sure we are within k-max */
					if((norm <= kmax*kmax) && (k_squared > 0.0)) {

						gaussian = exp(-k_squared/(4.0*alpha*alpha))/k_squared;

						/* structure factor */
						SF_real = 0; SF_imaginary = 0;
						for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
							for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

								if(!atom_ptr->frozen) {

									/* the inner product of the position vector and the k vector */
									for(p = 0, position_product = 0; p < 3; p++)
										position_product += k[p]*atom_ptr->pos[p];

									SF_real += atom_ptr->charge*cos(position_product);
									SF_imaginary += atom_ptr->charge*sin(position_product);

									/* include the self-interaction term  */
									/* this only need to happen once, on the first k-vec */
									if(!l[0] && !l[1] && (l[2] == 1)) {
										atom_ptr->es_self_point_energy = alpha*atom_ptr->charge*atom_ptr->charge/sqrt(M_PI);
										self_potential += atom_ptr->es_self_point_energy;
									}

								} /* !frozen */

							} /* atom */
						} /* molecule */
						recip_potential += gaussian*(SF_real*SF_real + SF_imaginary*SF_imaginary);

					} /* end if norm */

				} /* end for n */
			} /* end for m */
		} /* end for l */

		recip_potential *= 4.0*M_PI/system->pbc->volume;
		potential = recip_potential - self_potential;

	} /* ! wolf */

	return(potential);

}


/* real space sum */
double coulombic_real(system_t *system) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	pair_t *pair_ptr;
	double alpha, r, erfc_term, gaussian_term;
	double potential;
	double potential_classical, potential_fh_second_order, potential_fh_fourth_order;
	double first_derivative, second_derivative, third_derivative, fourth_derivative;
	double reduced_mass;

	alpha = system->ewald_alpha;

	potential = 0;
	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {


				if(pair_ptr->recalculate_energy) {

					pair_ptr->es_real_energy = 0;

					if(!pair_ptr->frozen) {

						r = pair_ptr->rimg;

						if(!((r > system->pbc->cutoff) || pair_ptr->es_excluded)) {	/* unit cell part */

							erfc_term = erfc(alpha*r);
							gaussian_term = exp(-alpha*alpha*r*r);

							potential_classical = atom_ptr->charge*pair_ptr->atom->charge*erfc_term/r;
							if(system->wolf)
								potential_classical -= atom_ptr->charge*pair_ptr->atom->charge*erfc(alpha*system->pbc->cutoff)/system->pbc->cutoff;
							pair_ptr->es_real_energy += potential_classical;

							if(system->feynman_hibbs) {

								reduced_mass = AMU2KG*molecule_ptr->mass*pair_ptr->molecule->mass/(molecule_ptr->mass+pair_ptr->molecule->mass);

								/* FIRST DERIVATIVE */
								first_derivative = -2.0*alpha*gaussian_term/(r*sqrt(M_PI)) - erfc_term/(r*r);

								/* SECOND DERIVATIVE */
								second_derivative = (4.0/sqrt(M_PI))*gaussian_term*(pow(alpha, 3.0) + pow(r, -2.0));
								second_derivative += 2.0*erfc_term/pow(r, 3.0);

								potential_fh_second_order = pow(METER2ANGSTROM, 2.0)*(HBAR*HBAR/(24.0*KB*system->temperature*reduced_mass))*(second_derivative + 2.0*first_derivative/r);
								pair_ptr->es_real_energy += potential_fh_second_order;

								if(system->feynman_hibbs_order >= 4) {

									/* THIRD DERIVATIVE */
									third_derivative = (gaussian_term/sqrt(M_PI))*(-8.0*pow(alpha, 5.0)*r - 8.0*pow(alpha, 3.0)/r - 12.0*alpha/pow(r, 3.0));
									third_derivative -= 6.0*erfc(alpha*r)/pow(r, 4.0);

									/* FOURTH DERIVATIVE */
									fourth_derivative = (gaussian_term/sqrt(M_PI))*( 8.0*pow(alpha, 5.0) + 16.0*pow(alpha, 7.0)*r*r + 32.0*pow(alpha, 3.0)/pow(r, 2.0) + 48.0/pow(r, 4.0) );
									fourth_derivative += 24.0*erfc_term/pow(r, 5.0);

									potential_fh_fourth_order = pow(METER2ANGSTROM, 4.0)*(pow(HBAR, 4.0)/(1152.0*pow(KB*system->temperature*reduced_mass, 2.0)))*(15.0*first_derivative/pow(r, 3.0) + 4.0*third_derivative/r + fourth_derivative);
									pair_ptr->es_real_energy += potential_fh_fourth_order;
								}

							}

						} else if(pair_ptr->es_excluded) /* calculate the self-intra part */
							pair_ptr->es_self_intra_energy = atom_ptr->charge*pair_ptr->atom->charge*erf(alpha*pair_ptr->r)/pair_ptr->r;

					} /* frozen */

				} /* recalculate */

				/* sum all of the pairwise terms */
				potential += pair_ptr->es_real_energy - pair_ptr->es_self_intra_energy;

			} /* pair */
		} /* atom */
	} /* molecule */


	return(potential);

}

/* allow for a neutralizing background, if necessary */
double coulombic_background(system_t *system) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	double alpha, potential;

	alpha = system->ewald_alpha;

	potential = 0;
	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
			potential += atom_ptr->charge;

	potential = M_PI*potential*potential/(2.0*alpha*alpha);

	return(potential);

}

/* no ewald summation - regular accumulation of Coulombic terms without out consideration of PBC */
/* only used by surface module */
double coulombic_nopbc(molecule_t * molecules) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	pair_t *pair_ptr;
	double pe, total_pe;

	total_pe = 0;
	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if(!pair_ptr->es_excluded) {
					pe = atom_ptr->charge*pair_ptr->atom->charge/pair_ptr->r;
					total_pe += pe;
				}


			}
		}
	}

	return(total_pe);

}


