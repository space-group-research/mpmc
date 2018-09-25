/* 

@2010, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

/* dreiding potential */
double dreiding(system_t *system) {
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr;
    double gamma, r_over_sigma, termexp, term6, potential, potential_classical;
#ifdef XXX
    double first_derivative, second_derivative, third_derivative, fourth_derivative;
    double potential, potential_classical, potential_fh_second_order, potential_fh_fourth_order;
    double reduced_mass;
    double sig3, sig_cut, sig_cut3, sig_cut9;
#endif /* XXX */

    gamma = DREIDING_GAMMA;

    potential = 0;
    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
                if (pair_ptr->recalculate_energy) {
                    pair_ptr->rd_energy = 0;

                    /* make sure we're not excluded or beyond the cutoff */
                    if (!((pair_ptr->rimg > system->pbc->cutoff) || pair_ptr->rd_excluded || pair_ptr->frozen)) {
                        r_over_sigma = pair_ptr->rimg / pair_ptr->sigma;

                        /* the DREIDING potential */
                        term6 = pow(r_over_sigma, -6);
                        term6 *= gamma / (gamma - 6.0);

                        if (pair_ptr->attractive_only)
                            termexp = 0;
                        else {
                            if (pair_ptr->rimg < 0.4 * pair_ptr->sigma)
                                termexp = MAXVALUE;
                            else {
                                termexp = exp(gamma * (1.0 - r_over_sigma));
                                termexp *= (6.0 / (gamma - 6.0));
                            }
                        }
                        potential_classical = pair_ptr->epsilon * (termexp - term6);

                        pair_ptr->rd_energy += potential_classical;

/* XXX */
/* need to do fh for dreiding */
#ifdef XXX
                        if (system->feynman_hibbs) {
                            reduced_mass = AMU2KG * molecule_ptr->mass * pair_ptr->molecule->mass / (molecule_ptr->mass + pair_ptr->molecule->mass);

                            /* FIRST DERIVATIVE */
                            first_derivative = -24.0 * pair_ptr->epsilon * (2.0 * term12 - term6) / pair_ptr->rimg;

                            /* SECOND DERIVATIVE */
                            second_derivative = 24.0 * pair_ptr->epsilon * (26.0 * term12 - 7.0 * term6) / pow(pair_ptr->rimg, 2);

                            potential_fh_second_order = pow(METER2ANGSTROM, 2) * (HBAR * HBAR / (24.0 * KB * system->temperature * reduced_mass)) * (second_derivative + 2.0 * first_derivative / pair_ptr->rimg);
                            pair_ptr->rd_energy += potential_fh_second_order;

                            if (system->feynman_hibbs_order >= 4) {
                                /* THIRD DERIVATIVE */
                                third_derivative = -1344.0 * pair_ptr->epsilon * (6.0 * term12 - term6) / pow(pair_ptr->rimg, 3);

                                /* FOURTH DERIVATIVE */
                                fourth_derivative = 12096.0 * pair_ptr->epsilon * (10.0 * term12 - term6) / pow(pair_ptr->rimg, 4);

                                potential_fh_fourth_order = pow(METER2ANGSTROM, 4) * (pow(HBAR, 4) / (1152.0 * pow(KB * system->temperature * reduced_mass, 2))) * (15.0 * first_derivative / pow(pair_ptr->rimg, 3) + 4.0 * third_derivative / pair_ptr->rimg + fourth_derivative);
                                pair_ptr->rd_energy += potential_fh_fourth_order;
                            }
                        }

#endif /* XXX */
                    }

/* XXX need to derive lrc for dreiding */
#ifdef XXX
                    /* include the long-range correction */
                    if (!(pair_ptr->rd_excluded || pair_ptr->frozen) && (pair_ptr->lrc == 0.0) && system->rd_lrc) {
                        sig_cut = fabs(pair_ptr->sigma) / system->pbc->cutoff;
                        sig3 = pow(fabs(pair_ptr->sigma), 3);
                        sig_cut3 = pow(sig_cut, 3);
                        sig_cut9 = pow(sig_cut, 9);
                        pair_ptr->lrc = ((-8.0 / 3.0) * M_PI * pair_ptr->epsilon * sig3) * sig_cut3 / system->pbc->volume;
                    }
#endif /* XXX */

                } /* if recalculate */

                /* sum all of the pairwise terms */
                potential += pair_ptr->rd_energy + pair_ptr->lrc;

            } /* pair */
        }     /* atom */
    }         /* molecule */

    return (potential);
}

/* same as above, but no periodic boundary conditions */
double dreiding_nopbc(molecule_t *molecules) {
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr;
    double gamma, r_over_sigma, termexp, term6;
    double potential;

    gamma = DREIDING_GAMMA;

    for (molecule_ptr = molecules, potential = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
                /* make sure we're not excluded or beyond the cutoff */
                if (!pair_ptr->rd_excluded) {
                    r_over_sigma = pair_ptr->r / pair_ptr->sigma;

                    /* the DREIDING potential */
                    term6 = pow(r_over_sigma, -6);
                    term6 *= gamma / (gamma - 6.0);

                    if (pair_ptr->attractive_only)
                        termexp = 0;
                    else {
                        if (pair_ptr->rimg < 0.35 * pair_ptr->sigma)
                            termexp = MAXVALUE;
                        else {
                            termexp = exp(gamma * (1.0 - r_over_sigma));
                            termexp *= (6.0 / (gamma - 6.0));
                        }
                    }
                    potential += pair_ptr->epsilon * (termexp - term6);
                }

            } /* pair */
        }     /* atom */
    }         /* molecule */

    return (potential);
}
