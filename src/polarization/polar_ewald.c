// Keith McLaughlin
// University of South Florida
// 7 June 2012

#include <mc.h>
#define OneOverSqrtPi 0.5641895835477562869480794515607725858440506293289988
#define SqrtPi 1.77245385091

void zero_out(molecule_t *m) {
    molecule_t *mptr;
    atom_t *aptr;
    int p;

    //zero out the electric field for each site
    for (mptr = m; mptr; mptr = mptr->next) {
        for (aptr = mptr->atoms; aptr; aptr = aptr->next) {
            for (p = 0; p < 3; p++) {
                aptr->ef_static[p] = 0.0;
                aptr->ef_static_self[p] = 0.0;
            }
        }
    }

    return;
}

//damping term (see e.g. Souaille et al., Comp. Phys. Comm. 180 276-301) below eq (9).
//signs are intentionally reversed (different convention)
double damp_factor(double t, int i) {
    double poo;

    poo = 1.0 + t + 0.5 * t * t;
    if (i == 3) poo += t * t * t / 6.0;

    return poo * exp(-t);
}

void real_term(system_t *system) {
    molecule_t *mptr;
    atom_t *aptr;
    pair_t *pptr;
    int p;
    double r, r2, factor, a;
    a = system->polar_ewald_alpha;  //some ambiguity between ea and ea^2 across the literature

    for (mptr = system->molecules; mptr; mptr = mptr->next) {
        for (aptr = mptr->atoms; aptr; aptr = aptr->next) {
            for (pptr = aptr->pairs; pptr; pptr = pptr->next) {  //for each pair
                if (pptr->frozen) continue;                      //if the pair is frozen (i.e. MOF-MOF interaction) it doesn't contribute to polar
                r = pptr->rimg;
                if ((r > system->pbc->cutoff) || (r == 0.0)) continue;  //if outside cutoff sphere (not sure why r==0 ever) -> skip
                r2 = r * r;
                if (pptr->es_excluded) {
                    //need to subtract self-term (interaction between a site and a neighbor's screening charge (on the same molecule)
                    factor = (2.0 * a * OneOverSqrtPi * exp(-a * a * r2) * r - erf(a * r)) / (r * r2);
                    for (p = 0; p < 3; p++) {
                        aptr->ef_static[p] += factor * pptr->atom->charge * pptr->dimg[p];
                        pptr->atom->ef_static[p] -= factor * aptr->charge * pptr->dimg[p];
                    }
                }       //excluded
                else {  //not excluded

                    factor = (2.0 * a * OneOverSqrtPi * exp(-a * a * r2) * r + erfc(a * r)) / (r2 * r);
                    for (p = 0; p < 3; p++) {  // for each dim, add e-field contribution for the pair
                        aptr->ef_static[p] += factor * pptr->atom->charge * pptr->dimg[p];
                        pptr->atom->ef_static[p] -= factor * aptr->charge * pptr->dimg[p];
                    }
                }  //excluded else
            }      //ptr
        }          //aptr
    }              //mptr

    return;
}

//we deviate from drexel's treatment, and instead do a trig identity to get from a pairwise sum to two atomwise rums
//or ignore drexel, and derive this term from eq (29) in nymand and linse
void recip_term(system_t *system) {
    molecule_t *mptr;
    atom_t *aptr;
    int p, q, l[3], kmax;
    double ea, k[3], k2, kweight[3], float1, float2;
    ea = system->polar_ewald_alpha;  //actually sqrt(ea)
    kmax = system->ewald_kmax;

    //k-space sum (symmetry for k -> -k, so we sum over hemisphere, avoiding double-counting on the face)
    for (l[0] = 0; l[0] <= kmax; l[0]++) {
        for (l[1] = (!l[0] ? 0 : -kmax); l[1] <= kmax; l[1]++) {
            for (l[2] = ((!l[0] && !l[1]) ? 1 : -kmax); l[2] <= kmax; l[2]++) {
                // if |l|^2 > kmax^2, then it doesn't contribute (outside the k-space cutoff)
                if (iidotprod(l, l) > kmax * kmax) continue;

                for (p = 0; p < 3; p++) {
                    for (q = 0, k[p] = 0; q < 3; q++)
                        k[p] += 2.0 * M_PI * system->pbc->reciprocal_basis[p][q] * l[q];
                }
                k2 = dddotprod(k, k);

                kweight[0] = k[0] / k2 * exp(-k2 / (4.0 * ea * ea));
                kweight[1] = k[1] / k2 * exp(-k2 / (4.0 * ea * ea));
                kweight[2] = k[2] / k2 * exp(-k2 / (4.0 * ea * ea));

                float1 = float2 = 0;
                for (mptr = system->molecules; mptr; mptr = mptr->next)
                    for (aptr = mptr->atoms; aptr; aptr = aptr->next) {
                        float1 += aptr->charge * cos(dddotprod(k, aptr->pos));
                        float2 += aptr->charge * sin(dddotprod(k, aptr->pos));
                    }

                for (mptr = system->molecules; mptr; mptr = mptr->next)
                    for (aptr = mptr->atoms; aptr; aptr = aptr->next) {
                        for (p = 0; p < 3; p++) {
                            aptr->ef_static[p] += kweight[p] * sin(dddotprod(k, aptr->pos)) * float1;
                            aptr->ef_static[p] -= kweight[p] * cos(dddotprod(k, aptr->pos)) * float2;
                        }
                    }

            }  //l2
        }      //l1
    }          //l0

    for (mptr = system->molecules; mptr; mptr = mptr->next) {
        for (aptr = mptr->atoms; aptr; aptr = aptr->next) {
            for (p = 0; p < 3; p++) {
                //factor of 2 more, since we only summed over hemisphere
                aptr->ef_static[p] *= 8.0 * M_PI / system->pbc->volume;
            }
        }
    }

    return;
}

//set zeroth iteration dipoles
void init_dipoles_ewald(system_t *system) {
    molecule_t *mptr;
    atom_t *aptr;
    int p;

    for (mptr = system->molecules; mptr; mptr = mptr->next)
        for (aptr = mptr->atoms; aptr; aptr = aptr->next)
            for (p = 0; p < 3; p++) {
                aptr->old_mu[p] = 0;
                aptr->new_mu[p] = aptr->mu[p] = aptr->polarizability * aptr->ef_static[p];
            }

    return;
}

//reset the ef_induced values
void clear_ef_induced(system_t *system) {
    molecule_t *mptr;
    atom_t *aptr;
    int p;

    for (mptr = system->molecules; mptr; mptr = mptr->next)
        for (aptr = mptr->atoms; aptr; aptr = aptr->next)
            for (p = 0; p < 3; p++)
                aptr->ef_induced[p] = 0;

    return;
}

//only calculate the static e-field via ewald
//see http://www.pages.drexel.edu/~cfa22/msim/node50.html
void ewald_estatic(system_t *system) {
    //calculate static e-field
    // no need to zero out dipoles; this is done in polar.c

    recip_term(system);
    real_term(system);

    return;
}

void induced_real_term(system_t *system) {
    molecule_t *mptr;
    atom_t *aptr;
    pair_t *pptr;
    double erfcar, expa2r2, r, ir, ir3, ir5;
    double T;                              //dipole-interaction tensor component
    double s1, s2;                         //common term (s_2 via eq 10. JCP 133 243101)
    int p, q;                              //dimensions
    double a = system->polar_ewald_alpha;  //ewald damping
    double l = system->polar_damp;         //polar damping

    for (mptr = system->molecules; mptr; mptr = mptr->next) {
        for (aptr = mptr->atoms; aptr; aptr = aptr->next) {
            for (pptr = aptr->pairs; pptr; pptr = pptr->next) {
                if (aptr->polarizability == 0 || pptr->atom->polarizability == 0) continue;  //don't waste CPU time
                if (pptr->rimg > system->pbc->cutoff) continue;                              //if outside cutoff sphere skip
                //some things we'll need
                r = pptr->rimg;
                ir = 1.0 / r;
                ir3 = ir * ir * ir;
                ir5 = ir * ir * ir3;
                erfcar = erfc(a * r);
                expa2r2 = exp(-a * a * r * r);

                //E_static_realspace_i = sum(i!=j) d_xi d_xj erfc(a*r)/r u_j
                s2 = erfcar + 2.0 * a * r * OneOverSqrtPi * expa2r2 + 4.0 * a * a * a * r * r * r / 3.0 * OneOverSqrtPi * expa2r2 - damp_factor(l * r, 3);

                for (p = 0; p < 3; p++) {
                    for (q = p; q < 3; q++) {  //it's symmetric!

                        if (p == q)
                            s1 = erfcar + 2.0 * a * r * OneOverSqrtPi * expa2r2 - damp_factor(l * r, 2);
                        else
                            s1 = 0;

                        //real-space dipole interaction tensor
                        T = 3.0 * pptr->dimg[p] * pptr->dimg[q] * s2 * ir5 - s1 * ir3;

                        aptr->ef_induced[p] += T * pptr->atom->mu[q];
                        pptr->atom->ef_induced[p] += T * aptr->mu[q];

                        if (p != q) {
                            aptr->ef_induced[q] += T * pptr->atom->mu[p];
                            pptr->atom->ef_induced[q] += T * aptr->mu[p];
                        }

                    }  //loop over q dim
                }      //loop over p dim

            }  //pptr loop
        }      //aptr loop
    }          //mptr loop

    return;
}

void induced_recip_term(system_t *system) {
    molecule_t *mptr;
    atom_t *aptr;
    atom_t **aarray = NULL;
    int i, j, N, l[3];
    int p, q;
    double Psin, Pcos, kweight;
    double a = system->polar_ewald_alpha;
    double kmax = system->ewald_kmax;
    double dotprod1, dotprod2, k[3], k2;

    //make atom array
    N = 0;
    for (mptr = system->molecules; mptr; mptr = mptr->next) {
        for (aptr = mptr->atoms; aptr; aptr = aptr->next) {
            aarray = realloc(aarray, sizeof(atom_t *) * (N + 1));
            memnullcheck(aarray, sizeof(atom_t *) * (N + 1), __LINE__ - 1, __FILE__);
            aarray[N] = aptr;
            N++;
        }
    }

    //k-space sum (symmetry for k -> -k, so we sum over hemisphere, avoiding double-counting on the face)
    for (l[0] = 0; l[0] <= kmax; l[0]++)
        for (l[1] = (!l[0] ? 0 : -kmax); l[1] <= kmax; l[1]++)
            for (l[2] = ((!l[0] && !l[1]) ? 1 : -kmax); l[2] <= kmax; l[2]++) {
                // if |l|^2 > kmax^2, then it doesn't contribute (outside the k-space cutoff)
                if (iidotprod(l, l) > kmax * kmax) continue;

                for (p = 0; p < 3; p++) {
                    for (q = 0, k[p] = 0; q < 3; q++)
                        k[p] += 2.0 * M_PI * system->pbc->reciprocal_basis[p][q] * l[q];
                }
                k2 = dddotprod(k, k);

                for (p = 0; p < 3; p++)
                    kweight = 8.0 * M_PI / system->pbc->volume * exp(-k2 / (4.0 * a * a)) / k2 * k[p];

                //calculate Pcos, Psin for this k-point
                Pcos = Psin = 0;
                for (j = 0; j < N; j++) {
                    dotprod1 = dddotprod(k, aarray[j]->mu);
                    dotprod2 = dddotprod(k, aarray[j]->pos);
                    Pcos += dotprod1 * cos(dotprod2);
                    Psin += dotprod1 * sin(dotprod2);
                }

                //calculate ef_induced over atom array
                for (i = 0; i < N; i++) {
                    //for each cartesian dimension
                    for (p = 0; p < 3; p++) {
                        dotprod1 = dddotprod(k, aarray[i]->pos);
                        aarray[i]->ef_induced[p] += kweight * (-sin(dotprod1) * Psin - cos(dotprod1) * Pcos);

                    }  //dim
                }      //ef_incuded over atom array

            }  //kspace

    free(aarray);

    return;
}

void induced_corr_term(system_t *system) {
    molecule_t *mptr;
    atom_t *aptr;
    int p;
    double a = system->polar_ewald_alpha;
    double totalmu[3];

    for (p = 0; p < 3; p++) totalmu[p] = 0;
    for (mptr = system->molecules; mptr; mptr = mptr->next)
        for (aptr = mptr->atoms; aptr; aptr = aptr->next)
            for (p = 0; p < 3; p++)
                totalmu[p] += aptr->mu[p];

    //other term
    for (mptr = system->molecules; mptr; mptr = mptr->next)
        for (aptr = mptr->atoms; aptr; aptr = aptr->next)
            for (p = 0; p < 3; p++)

                aptr->ef_induced[p] += -4.0 * M_PI / (3.0 * system->pbc->volume) * totalmu[p] + 4.0 * a * a * a / (3.0 * SqrtPi) * aptr->mu[p];

    return;
}

void new_dipoles(system_t *system, int count) {
    molecule_t *mptr;
    atom_t *aptr;
    int p;

    for (mptr = system->molecules; mptr; mptr = mptr->next) {
        for (aptr = mptr->atoms; aptr; aptr = aptr->next) {
            for (p = 0; p < 3; p++) {
                //set dipoles
                aptr->old_mu[p] = aptr->mu[p];
                if (system->polar_sor) {
                    aptr->new_mu[p] = aptr->polarizability * (aptr->ef_static[p] + aptr->ef_induced[p]);
                    aptr->new_mu[p] = aptr->mu[p] = system->polar_gamma * aptr->new_mu[p] + (1.0 - system->polar_gamma) * aptr->old_mu[p];
                } else if (system->polar_esor) {
                    aptr->new_mu[p] = aptr->polarizability * (aptr->ef_static[p] + aptr->ef_induced[p]);
                    aptr->new_mu[p] = aptr->mu[p] = (1.0 - exp(-system->polar_gamma * (count + 1))) * aptr->new_mu[p] +
                                                    exp(-system->polar_gamma * (count + 1)) * aptr->old_mu[p];
                } else {
                    //if no sor, still need new_mu for polar_palmo
                    aptr->mu[p] = aptr->new_mu[p] = aptr->polarizability * (aptr->ef_static[p] + aptr->ef_induced[p]);
                }
            }
        }
    }

    return;
}

void ewald_palmo_contraction(system_t *system) {
    molecule_t *mptr;
    atom_t *aptr;
    int p;

    //set induced field to zero
    clear_ef_induced(system);

    //calculate induced field
    induced_real_term(system);
    induced_recip_term(system);
    induced_corr_term(system);

    for (mptr = system->molecules; mptr; mptr = mptr->next) {
        for (aptr = mptr->atoms; aptr; aptr = aptr->next) {
            if (aptr->polarizability == 0) continue;
            for (p = 0; p < 3; p++) {
                aptr->ef_induced_change[p] =  //current induced - last induced (backed out from dipole values)
                    aptr->ef_induced[p] - (aptr->new_mu[p] / aptr->polarizability - aptr->ef_static[p]);
            }
        }
    }

    return;
}

//do full polarization calculation using ewald
//see nymand and linse jcp 112 6152 (2000)
void ewald_full(system_t *system) {
    //int max_iter=system->polar_max_iter;  (unused variable)
    int keep_iterating, iteration_counter;

    //calculate static e-field
    zero_out(system->molecules);
    recip_term(system);
    real_term(system);

    //calculate induced e-field
    init_dipoles_ewald(system);

    keep_iterating = 1;
    iteration_counter = 0;
    while (keep_iterating) {
        if (iteration_counter >= MAX_ITERATION_COUNT && system->polar_precision) {
            system->iter_success = 1;
            return;
        }

        //set induced field to zero
        clear_ef_induced(system);

        //calculate induced field
        induced_real_term(system);
        induced_recip_term(system);
        induced_corr_term(system);

        if (system->polar_rrms || system->polar_precision > 0)
            calc_dipole_rrms(system);

        //recalculate dipoles using new induced field
        new_dipoles(system, iteration_counter);

        keep_iterating = are_we_done_yet(system, iteration_counter);

        if (system->polar_palmo && !keep_iterating)  //if last iteration
            ewald_palmo_contraction(system);

        iteration_counter++;
    }

    return;
}
