// Keith McLaughlin
// University of South Florida
// 7 June 2012

#include <mc.h>
#define OneOverSqrtPi 0.5641895835477562869480794515607725858440506293289988

void zero_out ( molecule_t * m ) {
	molecule_t * mptr;
	atom_t * aptr;
	int p;

	//zero out the electric field for each site
	for ( mptr = m; mptr; mptr=mptr->next ) {
		for ( aptr = mptr->atoms; aptr; aptr=aptr->next ) {
			for ( p=0; p<3; p++ ) {
				aptr->ef_static[p] = 0.0;
				aptr->ef_static_self[p] = 0.0;
			}
		}
	}

	return;
}


//damping term (see e.g. Souaille et al., Comp. Phys. Comm. 180 276-301) below eq (9).
//signs are intentionally reversed (different convention)
double damp_factor ( double t, int i ) {
	double poo;

	poo = 1.0 + t + 0.5*t*t;
	if ( i == 3 ) poo += t*t*t/6.0;

	return poo*exp(-t);
}

void real_term ( system_t * system ) {
	molecule_t * mptr;
	atom_t * aptr;
	pair_t * pptr;
	int p;
	double r, r2, factor, ea;
	ea = system->ewald_alpha; //some ambiguity between ea and ea^2 across the literature

	for (mptr = system->molecules; mptr; mptr=mptr->next ) {
		for (aptr = mptr->atoms; aptr; aptr=aptr->next ) {
			for (pptr = aptr->pairs; pptr; pptr=pptr->next ) { //for each pair
				if (pptr->frozen) continue; //if the pair is frozen (i.e. MOF-MOF interaction) it doesn't contribute to polar
				r = pptr->rimg;
				if ( (r > system->pbc->cutoff) || (r == 0.0) ) continue; //if outside cutoff sphere (not sure why r==0 ever) -> skip
				r2 = r*r; 
				if (pptr->es_excluded) {
					//need to subtract self-term (interaction between a site and a neighbor's screening charge (on the same molecule)
					factor = (2.0*ea*OneOverSqrtPi*exp(-ea*ea*r2)*r - erf(ea*r))/(r*r2);
					for ( p=0; p<3; p++ ) {
						aptr->ef_static_self[p] += factor*pptr->atom->charge * pptr->dimg[p];
						pptr->atom->ef_static_self[p] -= factor*aptr->charge * pptr->dimg[p];
					}
				} //excluded
				else { //not excluded
					r2 = r*r;
					factor = (2.0*ea*OneOverSqrtPi*exp(-ea*ea*r2)*r + erfc(ea*r))/(r2*r);
					for ( p=0; p<3; p++ ) { // for each dim, add e-field contribution for the pair
						aptr->ef_static[p] += factor*pptr->atom->charge * pptr->dimg[p];
						pptr->atom->ef_static[p] -= factor*aptr->charge * pptr->dimg[p];
					}
				} //excluded else
			} //ptr
		} //aptr
	} //mptr

	return;
}

//we deviate from drexel's treatment, and instead do a trig identity to get from a pairwise sum to two atomwise rums
//or ignore drexel, and derive this term from eq (29) in nymand and linse
void recip_term ( system_t * system ) {
	molecule_t * mptr;
	atom_t * aptr;
	pair_t * pptr;
	int p, q, l[3], kmax;
	double r, ea, k[3], k2, kdotr, kweight[3], float1, float2;
	ea = system->ewald_alpha; //actually sqrt(ea)
	kmax = system->ewald_kmax;

	//k-space sum (symmetry for k -> -k, so we sum over hemisphere, avoiding double-counting on the face)
	for (l[0] = 0; l[0] <= kmax; l[0]++) {
		for (l[1] = (!l[0] ? 0 : -kmax); l[1] <= kmax; l[1]++) {
			for (l[2] = ((!l[0] && !l[1]) ? 1 : -kmax); l[2] <= kmax; l[2]++) {

				// if |l|^2 > kmax^2, then it doesn't contribute (outside the k-space cutoff)
				if ( iidotprod(l,l) > kmax*kmax ) continue; 

				for ( p=0; p<3; p++ ) //make recip lattice vector
					k[p] = 2.0*M_PI*didotprod(system->pbc->reciprocal_basis[p],l);

				k2 = dddotprod(k,k); // |k|^2

				kweight[0] = k[0]/k2 * exp(-k2/(4.0*ea*ea));
				kweight[1] = k[1]/k2 * exp(-k2/(4.0*ea*ea));
				kweight[2] = k[2]/k2 * exp(-k2/(4.0*ea*ea));

				float1 = float2 = 0;
				for (mptr = system->molecules; mptr; mptr=mptr->next )
					for (aptr = mptr->atoms; aptr; aptr=aptr->next ) {
						float1 += aptr->charge * cos(dddotprod(k,aptr->pos));
						float2 += aptr->charge * sin(dddotprod(k,aptr->pos));
					}

				for (mptr = system->molecules; mptr; mptr=mptr->next )
					for (aptr = mptr->atoms; aptr; aptr=aptr->next ) {
						for ( p=0; p<3; p++ ) {
							aptr->ef_static[p] += kweight[p] * sin(dddotprod(k,aptr->pos)) * float1;
							aptr->ef_static[p] -= kweight[p] * cos(dddotprod(k,aptr->pos)) * float2;
						}
					}

			} //l2
		} //l1
	} //l0

	for ( mptr=system->molecules; mptr; mptr=mptr->next ) {
		for ( aptr=mptr->atoms; aptr; aptr=aptr->next ) {
			for ( p=0; p<3; p++ ) {
				//factor of 2 more, since we only summed over hemisphere
				aptr->ef_static[p] *= 8.0*M_PI/system->pbc->volume; 
			}
		}
	}

	return;
}

//set zeroth iteration dipoles
void init_dipoles_ewald( system_t * system ) {
	molecule_t * mptr;
	atom_t * aptr;
	int p;

	for ( mptr=system->molecules; mptr; mptr=mptr->next ) {
		for ( aptr=mptr->atoms; aptr; aptr=aptr->next ) {
			for ( p=0; p<3; p++ ) {
				aptr->mu[p] = aptr->polarizability*aptr->ef_static[p];
				aptr->ef_induced[p] = 0;
			}
		}
	}
	return;
}

//only calculate the static e-field via ewald
//see http://www.pages.drexel.edu/~cfa22/msim/node50.html
void ewald_estatic ( system_t * system ) {

	//calculate static e-field
	zero_out(system->molecules);
	recip_term(system);
	real_term(system);

	return;
}



void induced_real_term(system_t * system) {
	molecule_t * mptr;
	atom_t * aptr;
	pair_t * pptr;
	double erfcar, expa2r2, r, ir, ir3, ir5, explr;
	double T; //dipole-interaction tensor component
	int p, q; //dimensions
	double a = system->ewald_alpha; //ewald damping
	double l = system->polar_damp; //polar damping

	//zero out this term for each atom ptr
	for ( mptr = system->molecules; mptr; mptr=mptr->next ) 
		for ( aptr = mptr->atoms; aptr; aptr=aptr->next ) 
			for ( p=0; p<3; p++ )  
				aptr->ef_induced[p] = 0;

	for ( mptr = system->molecules; mptr; mptr=mptr->next ) {
		for ( aptr = mptr->atoms; aptr; aptr=aptr->next ) {
			for ( pptr = aptr->pairs; pptr; pptr=pptr->next ) { 
				if ( aptr->polarizability == 0 || pptr->atom->polarizability == 0 ) continue; //don't waste CPU time

				//some things we'll need
				r=pptr->rimg;
				ir=1.0/r; ir3=ir*ir*ir; ir5=ir*ir*ir3;
				erfcar=erfc(system->ewald_alpha*r);
				expa2r2=exp(-system->ewald_alpha*system->ewald_alpha*r*r);

				//E_static_realspace_i = sum(i!=j) d_xi d_xj erfc(a*r)/r u_j 
				for ( p=0; p<3; p++ ) {
					for ( q=p; q<3; q++ ) {  //it's symmetric!
						T = pptr->dimg[p] * pptr->dimg[q] * ir5; //common term
						T *= 3*erfcar + (4*a*a*r*r+6)*a*r*OneOverSqrtPi*expa2r2 //common term (cont.)
								-3*damp_factor(l*r,3); //common damp term
						if ( p==q ) {
							T -= (erfcar+2*a*r*OneOverSqrtPi*expa2r2)*ir3; //kronecker delta term
							T += ir3*damp_factor(l*r,2); //kronecker damp term
						}

						aptr->ef_induced[p] += T*pptr->atom->mu[q];
						pptr->atom->ef_induced[p] += T*aptr->mu[q];

						if ( p!=q ) {
							aptr->ef_induced[q] += T*pptr->atom->mu[p];
							pptr->atom->ef_induced[q] += T*aptr->mu[p];
						}

					} //loop over q dim
				} //loop over p dim

			} //pptr loop
		} //aptr loop
	} //mptr loop
						
	return;
}

void induced_recip_term(system_t * system) {
	molecule_t * mptr;
	atom_t * aptr;
	int p, q, l[3], kmax;
	double r, ea, k[3], k2, kdotr, kweight[3], float1;
	double sinsum, cossum;
	ea = system->ewald_alpha; //actually sqrt(ea)
	kmax = system->ewald_kmax;

	//k-space sum (symmetry for k -> -k, so we sum over hemisphere, avoiding double-counting on the face)
	for (l[0] = 0; l[0] <= kmax; l[0]++) {
		for (l[1] = (!l[0] ? 0 : -kmax); l[1] <= kmax; l[1]++) {
			for (l[2] = ((!l[0] && !l[1]) ? 1 : -kmax); l[2] <= kmax; l[2]++) {

				// if |l|^2 > kmax^2, then it doesn't contribute (outside the k-space cutoff)
				if ( iidotprod(l,l) > kmax*kmax ) continue; 

				for ( p=0; p<3; p++ ) //make recip lattice vector
					k[p] = 2.0*M_PI*didotprod(system->pbc->reciprocal_basis[p],l);

				k2 = dddotprod(k,k); // |k|^2
				kweight[0] = k[0]/k2 * exp(-k2/(4.0*ea*ea));
				kweight[1] = k[1]/k2 * exp(-k2/(4.0*ea*ea));
				kweight[2] = k[2]/k2 * exp(-k2/(4.0*ea*ea));

				//start calculating first term
				cossum=0; sinsum=0;
				for ( mptr = system->molecules; mptr; mptr=mptr->next ) 
					for ( aptr = mptr->atoms; aptr; aptr=aptr->next ) {
						float1 = dddotprod(aptr->mu,k);
						cossum += float1*cos(dddotprod(k,aptr->pos));
						sinsum += float1*sin(dddotprod(k,aptr->pos));
					}
				for ( mptr = system->molecules; mptr; mptr=mptr->next ) 
					for ( aptr = mptr->atoms; aptr; aptr=aptr->next ) 
						for ( p=0; p<3; p++ ) {
							aptr->ef_induced[p] -= 8*M_PI/system->pbc->volume * kweight[p] * cos(dddotprod(k,aptr->pos)) * cossum;
							aptr->ef_induced[p] -= 8*M_PI/system->pbc->volume * kweight[p] * sin(dddotprod(k,aptr->pos)) * sinsum;
						}

				//second term
				for ( mptr = system->molecules; mptr; mptr=mptr->next ) 
					for ( aptr = mptr->atoms; aptr; aptr=aptr->next ) 
						for ( p=0; p<3; p++ )  //extra factor of two from symmetry
							aptr->ef_induced[p] += (8.0*M_PI/system->pbc->volume)*kweight[p]*dddotprod(aptr->mu,k);

			} //l2
		} //l1
	} //l0

	//other term
	for ( mptr = system->molecules; mptr; mptr=mptr->next ) 
		for ( aptr = mptr->atoms; aptr; aptr=aptr->next )
			for ( p=0; p<3; p++ ) 
				aptr->ef_induced[p] += 4.0*M_PI/(3.0*system->pbc->volume)*aptr->mu[p];

	return;
}

void induced_self_term(system_t * system) {
	molecule_t * mptr;
	atom_t * aptr;
	pair_t * pptr;
	double T; //dipole-interaction tensor for self terms
	double erfar, expa2r2, r, ir, ir2, ir3, ir4, ir5;
	double a = system->ewald_alpha;
	int p, q; //dimensions

	for ( mptr = system->molecules; mptr; mptr=mptr->next ) {
		for ( aptr = mptr->atoms; aptr; aptr=aptr->next ) {
			for ( pptr = aptr->pairs; pptr; pptr=pptr->next ) { 
				if ( !pptr->es_excluded ) continue; //same molecule only
				if ( aptr->polarizability == 0 || pptr->atom->polarizability == 0 ) continue; //don't waste CPU time
				//some things we'll need
				r=pptr->rimg;
				ir=1.0/r; ir2=ir*ir; ir3=ir2*ir; ir4=ir3*ir; ir5=ir4*ir;
				erfar=erf(system->ewald_alpha*r);
				expa2r2=exp(-system->ewald_alpha*system->ewald_alpha*r*r);

				for ( p=0; p<3; p++ ) {
					for ( q=p; q<3; q++ ) {  //it's symmetric!
						T = pptr->dimg[p]*pptr->dimg[q]*ir5;  //common term
						T *= 3*erfar - (2*a*a*r*r+3)*(2*a*r*expa2r2*OneOverSqrtPi);
						if ( p==q ) 
							T += -erfar*ir3 + 2*a*expa2r2*OneOverSqrtPi*ir2; //kronecker term

						aptr->ef_induced[p] -= T*pptr->atom->mu[q];
						pptr->atom->ef_induced[p] -= T*aptr->mu[q];
						if ( p!=q ) {
							aptr->ef_induced[q] -= T*pptr->atom->mu[p];
							pptr->atom->ef_induced[q] -= T*aptr->mu[p];
						}

					} //loop over q dim
				} //loop over p dim

			} //pptr loop
		} //aptr loop
	} //mptr loop
								
	return;
}

void new_dipoles(system_t * system) {
	molecule_t * mptr;
	atom_t * aptr;
	int p;

	for ( mptr = system->molecules; mptr; mptr=mptr->next ) {
		for ( aptr = mptr->atoms; aptr; aptr=aptr->next ) {
			for ( p=0; p<3; p++ ) {
				aptr->mu[p] = aptr->polarizability*(aptr->ef_static[p]+aptr->ef_induced[p]);
				aptr->ef_induced[p] = 0;
			}
		}
	}

	return;
}

//do full polarization calculation using ewald
//see nymand and linse jcp 112 6152 (2000)
void ewald_full ( system_t * system ) {

	int max_iter=system->polar_max_iter;
	int i;

	//calculate static e-field
	zero_out(system->molecules);
	recip_term(system);
	real_term(system);

	//calculate induced e-field
	init_dipoles_ewald(system);	

	for ( i=0; i<max_iter; i++ ) {
		induced_real_term(system);
		induced_recip_term(system);
		induced_self_term(system);
		new_dipoles(system);
	}

	return;
}

