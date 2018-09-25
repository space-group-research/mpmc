/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

/* Silvera-Goldman H2 potential */
double sg(system_t *system) {
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr;
    double rimg, r6, r8, r9, r10, r_rm;
    double repulsive_term, multipole_term, exponential_term;
    double first_r_diff_term, second_r_diff_term;
    double first_derivative, second_derivative;
    double potential_classical, potential_fh_second_order;
    double potential;
    double temperature;

    temperature = system->temperature;

    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
                if (pair_ptr->recalculate_energy) {
                    pair_ptr->rd_energy = 0;
                    rimg = pair_ptr->rimg;

                    if (rimg < system->pbc->cutoff) {
                        /* convert units to Bohr radii */
                        rimg /= AU2ANGSTROM;

                        /* classical pairwise part */
                        repulsive_term = exp(ALPHA - BETA * rimg - GAMMA * rimg * rimg);

                        r6 = pow(rimg, 6);
                        r8 = pow(rimg, 8);
                        r9 = pow(rimg, 9);
                        r10 = pow(rimg, 10);
                        multipole_term = C6 / r6 + C8 / r8 + C10 / r10 - C9 / r9;

                        r_rm = RM / rimg;
                        if (rimg < RM)
                            exponential_term = exp(-pow((r_rm - 1.0), 2));
                        else
                            exponential_term = 1.0;

                        potential_classical = (repulsive_term - multipole_term * exponential_term);
                        pair_ptr->rd_energy += potential_classical;

                        if (system->feynman_hibbs) {
                            /* FIRST DERIVATIVE */
                            first_derivative = (-BETA - 2.0 * GAMMA * rimg) * repulsive_term;
                            first_derivative += (6.0 * C6 / pow(rimg, 7) + 8.0 * C8 / pow(rimg, 9) - 9.0 * C9 / pow(rimg, 10) + 10.0 * C10 / pow(rimg, 11)) * exponential_term;
                            first_r_diff_term = (r_rm * r_rm - r_rm) / rimg;
                            first_derivative += -2.0 * multipole_term * exponential_term * first_r_diff_term;

                            /* SECOND DERIVATIVE */
                            second_derivative = (pow((BETA + 2.0 * GAMMA * rimg), 2) - 2.0 * GAMMA) * repulsive_term;
                            second_derivative += (-exponential_term) * (42.0 * C6 / pow(rimg, 8) + 72.0 * C8 / pow(rimg, 10) - 90.0 * C9 / pow(rimg, 11) + 110.0 * C10 / pow(rimg, 10));
                            second_derivative += exponential_term * first_r_diff_term * (12.0 * C6 / pow(rimg, 7) + 16.0 * C8 / pow(rimg, 9) - 18.0 * C9 / pow(rimg, 10) + 20.0 * C10 / pow(rimg, 11));
                            second_derivative += exponential_term * pow(first_r_diff_term, 2) * 4.0 * multipole_term;
                            second_r_diff_term = (3.0 * r_rm * r_rm - 2.0 * r_rm) / (rimg * rimg);
                            second_derivative += exponential_term * second_r_diff_term * 2.0 * multipole_term;

                            potential_fh_second_order = pow(METER2ANGSTROM, 2) * (HBAR * HBAR / (24.0 * KB * temperature * (AMU2KG * molecule_ptr->mass))) * (second_derivative + 2.0 * first_derivative / rimg);
                            pair_ptr->rd_energy += potential_fh_second_order;
                        }

                        /* convert units from Hartrees back to Kelvin */
                        pair_ptr->rd_energy *= HARTREE2KELVIN;
                    }

                } /* recalculate */
            }     /* pair */
        }         /* atom */
    }             /* molecule */

    potential = 0;
    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
            for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next)
                potential += pair_ptr->rd_energy;

    return (potential);
}

/* same as above, but no periodic boundary conditions */
double sg_nopbc(molecule_t *molecules) {
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr;
    double r, r6, r8, r10, r9;
    double r_rm, r_rm_2, r_exp;
    double potential, multipole_term;
    double result, exp_result;

    for (molecule_ptr = molecules, potential = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
                if (pair_ptr->recalculate_energy) {
                    r = pair_ptr->r / AU2ANGSTROM;

                    r6 = pow(r, 6);
                    r8 = pow(r, 8);
                    r10 = pow(r, 10);
                    r9 = pow(r, 9);

                    multipole_term = C6 / r6 + C8 / r8 + C10 / r10 - C9 / r9;

                    if (r < RM) {
                        r_rm = RM / r;
                        r_rm -= 1.0;
                        r_rm_2 = pow(r_rm, 2);
                        r_rm_2 *= -1.0;
                        r_exp = exp(r_rm_2);

                        multipole_term *= r_exp;
                    }

                    result = ALPHA - BETA * r - GAMMA * r * r;

                    exp_result = exp(result);
                    exp_result -= multipole_term;
                    pair_ptr->rd_energy = HARTREE2KELVIN * exp_result;

                } /* recalculate */
            }     /* pair */
        }         /* atom */
    }             /* molecule */

    potential = 0;
    for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
            for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next)
                potential += pair_ptr->rd_energy;

    return (potential);
}
