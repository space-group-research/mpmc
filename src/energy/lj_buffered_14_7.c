#include <mc.h>

//Copyright 2013 Adam Hogan
//TO DO: long range correction

double lj_buffered_14_7(system_t *system) {
    double potential = 0.0, potential_classical;
    double r_over_sigma, first_term, second_term;

    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr;

    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
                if (pair_ptr->recalculate_energy) {
                    pair_ptr->rd_energy = 0;

                    /* make sure we're not excluded or beyond the cutoff */
                    if (!((pair_ptr->rimg > system->pbc->cutoff) || pair_ptr->rd_excluded || pair_ptr->frozen)) {
                        r_over_sigma = pair_ptr->rimg / pair_ptr->sigma;
                        first_term = pow(1.07 / (r_over_sigma + 0.07), 7);
                        second_term = (1.12 / (pow(r_over_sigma, 7) + 0.12) - 2);
                        potential_classical = pair_ptr->epsilon * first_term * second_term;
                        pair_ptr->rd_energy += potential_classical;
                    }
                }
                potential += pair_ptr->rd_energy;
            }
        }
    }

    return potential;
}

double lj_buffered_14_7_nopbc(system_t *system) {
    double potential = 0.0, potential_classical;
    double r_over_sigma, first_term, second_term;

    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr;

    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
                if (pair_ptr->recalculate_energy) {
                    pair_ptr->rd_energy = 0;

                    /* make sure we're not excluded or beyond the cutoff */
                    if (!((pair_ptr->rimg > system->pbc->cutoff) || pair_ptr->rd_excluded || pair_ptr->frozen)) {
                        r_over_sigma = pair_ptr->rimg / pair_ptr->sigma;
                        first_term = pow(1.07 / (r_over_sigma + 0.07), 7);
                        second_term = (1.12 / (pow(r_over_sigma, 7) + 0.12) - 2);
                        potential_classical = pair_ptr->epsilon * first_term * second_term;
                        pair_ptr->rd_energy += potential_classical;
                    }
                }
                potential += pair_ptr->rd_energy;
            }
        }
    }

    return potential;
}
