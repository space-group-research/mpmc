/* 

Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

double exp_fh_corr(system_t *system, molecule_t *molecule_ptr, pair_t *pair_ptr, int order, double pot) {
    double reduced_mass;
    double dE, d2E, d3E, d4E;  //energy derivatives
    double corr;
    double ir = 1.0 / pair_ptr->rimg;
    double ir2 = ir * ir;
    double ir3 = ir2 * ir;
    // double ir4 = ir3*ir;   (unused variable)

    if ((order != 2) && (order != 4)) return NAN;  //must be order 2 or 4

    reduced_mass = AMU2KG * molecule_ptr->mass * pair_ptr->molecule->mass /
                   (molecule_ptr->mass + pair_ptr->molecule->mass);

    dE = -pot / (2.0 * pair_ptr->epsilon);
    d2E = dE / (2.0 * pair_ptr->epsilon);

    //2nd order correction
    corr = M2A2 *
           (HBAR2 / (24.0 * KB * system->temperature * reduced_mass)) *
           (d2E + 2.0 * dE / pair_ptr->rimg);

    if (order >= 4) {
        d3E = -d2E / (2.0 * pair_ptr->epsilon);
        d4E = d3E / (2.0 * pair_ptr->epsilon);

        //4th order corection
        corr += M2A4 *
                (HBAR4 / (1152.0 * KB2 * system->temperature * system->temperature * reduced_mass * reduced_mass)) *
                (15.0 * dE * ir3 + 4.0 * d3E * ir + d4E);
    }

    return corr;
}

double exp_lrc_corr(system_t *system, atom_t *atom_ptr, pair_t *pair_ptr, double cutoff) {
    double eps = pair_ptr->epsilon;
    double rover2e = cutoff / (2.0 * eps);

    /* include the long-range correction */ /* I'm  not sure that I'm handling spectre pairs correctly */
    /* we can't use rd_excluded flag, since that disqualifies inter-molecular, but that DOES contribute to LRC */
    /* ALL OF THESE MUST BE TRUE TO PERFORM LRC CALCULATION */
    if ((pair_ptr->epsilon != 0 && pair_ptr->sigma != 0) &&                          //if these are zero, then we won't waste our time
        !(atom_ptr->spectre && pair_ptr->atom->spectre) &&                           //i think we want to disqualify s-s pairs
        !(pair_ptr->frozen) &&                                                       //disqualify frozen pairs
        ((pair_ptr->lrc == 0.0) || pair_ptr->last_volume != system->pbc->volume)) {  //LRC only changes if the volume change

        pair_ptr->last_volume = system->pbc->volume;

        return (8.0 * M_PI) * exp(1. - rover2e) * (cutoff * cutoff + 4.0 * eps * cutoff + 8.0 * eps * eps) * pair_ptr->sigma / system->pbc->volume;

    } else
        return pair_ptr->lrc;  //use stored value
}

double exp_lrc_self(system_t *system, atom_t *atom_ptr, double cutoff) {
    double eps = atom_ptr->epsilon;
    double rover2e = cutoff / (2.0 * eps);

    if (((atom_ptr->sigma != 0) && (atom_ptr->epsilon != 0)) &&  //non-zero parameters
        !(atom_ptr->frozen) &&                                   //not frozen
        !(atom_ptr->spectre)) {                                  //not spectre

        return (8.0 * M_PI) * exp(1. - rover2e) * (cutoff * cutoff + 4.0 * eps * cutoff + 8.0 * eps * eps) * atom_ptr->sigma / system->pbc->volume;
    }
    return 0;
}

double exp_crystal_self(system_t *system, atom_t *aptr, double cutoff) {
    double term = 0;
    double r;
    double eps = aptr->epsilon;
    int i[3], p, q;
    double a[3];
    // double prefactor;   (unused variable)

    if (aptr->sigma == 0 || aptr->epsilon == 0) return 0;  //skip if no LJ interaction

    for (i[0] = -(system->rd_crystal_order); i[0] <= system->rd_crystal_order; i[0]++)
        for (i[1] = -(system->rd_crystal_order); i[1] <= system->rd_crystal_order; i[1]++)
            for (i[2] = -(system->rd_crystal_order); i[2] <= system->rd_crystal_order; i[2]++) {
                if (!i[0] && !i[1] && !i[2]) continue;  //no (0,0,0)
                //calculate pair separation (atom with it's image)
                for (p = 0; p < 3; p++) {
                    a[p] = 0;
                    for (q = 0; q < 3; q++)
                        a[p] += system->pbc->basis[q][p] * i[q];
                }
                r = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

                if (r > cutoff) continue;             //too far away! will be included in LRC if enabled
                term += 0.5 * exp(-r / (2.0 * eps));  //multiply by 0.5 to get counting correct
            }

    return aptr->sigma * term;
}

double exp_repulsion(system_t *system) {
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr;
    double r, term;
    double potential, potential_classical, cutoff;
    int i[3], p, q;
    double a[3];

    //set the cutoff
    if (system->rd_crystal)
        cutoff = 2.0 * system->pbc->cutoff * ((double)system->rd_crystal_order - 0.5);
    else
        cutoff = system->pbc->cutoff;

    potential = 0;
    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
                if (pair_ptr->recalculate_energy) {
                    pair_ptr->rd_energy = 0;

                    // pair LRC
                    if (system->rd_lrc) pair_ptr->lrc = exp_lrc_corr(system, atom_ptr, pair_ptr, cutoff);

                    // to include a contribution, we require
                    if ((pair_ptr->rimg - SMALL_dR < cutoff)               //inside cutoff?
                        && (!pair_ptr->rd_excluded || system->rd_crystal)  //either not excluded OR rd_crystal is ON
                        && !pair_ptr->frozen)                              //not frozen
                    {
                        //loop over unit cells
                        if (system->rd_crystal) {
                            term = 0;
                            for (i[0] = -(system->rd_crystal_order); i[0] <= system->rd_crystal_order; i[0]++)
                                for (i[1] = -(system->rd_crystal_order); i[1] <= system->rd_crystal_order; i[1]++)
                                    for (i[2] = -(system->rd_crystal_order); i[2] <= system->rd_crystal_order; i[2]++) {
                                        if (!i[0] && !i[1] && !i[2] && pair_ptr->rd_excluded) continue;  //no i=j=k=0 for excluded pairs (intra-molecular)
                                        //calculate pair separation (atom with it's image)
                                        for (p = 0; p < 3; p++) {
                                            a[p] = 0;
                                            for (q = 0; q < 3; q++)
                                                a[p] += system->pbc->basis[q][p] * i[q];
                                            a[p] += atom_ptr->pos[p] - pair_ptr->atom->pos[p];
                                        }
                                        r = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

                                        if (r + SMALL_dR > cutoff) continue;
                                        term += exp(-r / (2.0 * pair_ptr->epsilon));
                                    }
                        } else  //otherwise, calculate as normal
                            term = exp(-pair_ptr->rimg / (2.0 * pair_ptr->epsilon));

                        potential_classical = pair_ptr->sigma * term;
                        pair_ptr->rd_energy += potential_classical;

                        if (system->feynman_hibbs)
                            pair_ptr->rd_energy +=
                                exp_fh_corr(system, molecule_ptr, pair_ptr, system->feynman_hibbs_order, potential_classical);

                    }  //count contributions

                } /* if recalculate */

                /* sum all of the pairwise terms */
                potential += pair_ptr->rd_energy + pair_ptr->lrc;

            } /* pair */
        }     /* atom */
    }         /* molecule */

    /* molecule self-energy for rd_crystal -> energy of molecule interacting with its periodic neighbors */

    if (system->rd_crystal)
        for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
            for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
                potential += exp_crystal_self(system, atom_ptr, cutoff);

    /* calculate self LRC interaction */
    if (system->rd_lrc)
        for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
            for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
                potential += exp_lrc_self(system, atom_ptr, cutoff);

    return (potential);
}

/* same as above, but no periodic boundary conditions */
double exp_repulsion_nopbc(system_t *system) {
    molecule_t *molecules = system->molecules;
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr;
    double potential;

    for (molecule_ptr = molecules, potential = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
                /* make sure we're not excluded or beyond the cutoff */
                if (!pair_ptr->rd_excluded)
                    potential += pair_ptr->sigma * exp(-pair_ptr->rimg / (2.0 * pair_ptr->epsilon));

            } /* pair */
        }     /* atom */
    }         /* molecule */

    return (potential);
}

#ifdef DEBUG
void test_exp(system_t *system) {
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr;
    char poo[MAXLINE];
    sprintf(poo,
            "%d.lj", system->step);
    FILE *fp = fopen(poo,
                     "w");

    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
            for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next)
                fprintf(fp,
                        "DEBUG_LJ: m_id %d %d a_id %d %d rimg %.3lf\n", molecule_ptr->id, pair_ptr->molecule->id, atom_ptr->id, pair_ptr->atom->id, pair_ptr->rimg);

    fclose(fp);

    return;
}
#endif
