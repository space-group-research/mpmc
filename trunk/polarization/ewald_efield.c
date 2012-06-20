// Keith McLaughlin
// University of South Florida
// 7 June 2012

#include <mc.h>
#include <fenv.h>

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

				//now that we have the K-vector, we need to calculate the corresponding contribution to each pair			
				for (mptr = system->molecules; mptr; mptr=mptr->next ) {
					for (aptr = mptr->atoms; aptr; aptr=aptr->next ) {
						for (pptr = aptr->pairs; pptr; pptr=pptr->next ) { //for each pair
							if ( pptr->frozen ) continue; //if the pair is frozen (i.e. MOF-MOF) it doens't contribute

							kdotr = dddotprod(k,pptr->dimg);
						
							//we're going to loop through here N*N*K times (atoms*atoms*kpts), let's calc sin as few times as possible
							float1 = sin(kdotr); 
							float2 = float1 * aptr->charge;
							float1 *= pptr->atom->charge;

							//assign e-fields (short by a factor of 8Pi/V) 
							for ( p=0; p<3; p++ ) {
								aptr->ef_static[p] += float1 * kweight[p];
								pptr->atom->ef_static[p] -= float2 * kweight[p];
							}

						} //pptr
					}//aptr
				}//mptr
				
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
void init_dipoles( system_t * system ) {
	molecule_t * mptr;
	atom_t * aptr;
	int p;

	for ( mptr=system->molecules; mptr; mptr=mptr->next ) {
		for ( aptr=mptr->atoms; aptr; aptr=aptr->next ) {
			for ( p=0; p<3; p++ )
				aptr->mu[p] = aptr->polarizability*aptr->ef_static[p];
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

	return;
}

void induced_recip_term(system_t * system) {

	return;
}

void induced_self_term(system_t * system) {

	return;
}

void new_dipoles(system_t * system) {


	return;
}

//do full polarization calculation using ewald
//see nymand and linse jcp 112 6152 (2000)
void ewald_full ( system_t * system ) {

	int max_iter=10; //temporary
	int i;

	//calculate static e-field
	zero_out(system->molecules);
	recip_term(system);
	real_term(system);

	//calculate induced e-field
	init_dipoles(system);	

	for ( i=0; i<max_iter; i++ ) {
		induced_real_term(system);
		induced_recip_term(system);
		induced_self_term(system);
		new_dipoles(system);
	}

	return;
}

