/* 

Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

double lj_fh_corr(system_t *system, molecule_t *molecule_ptr, pair_t *pair_ptr, int order, double term12, double term6) {
    double reduced_mass;
    double dE, d2E, d3E, d4E;  //energy derivatives
    double corr;
    double ir = 1.0 / pair_ptr->rimg;
    double ir2 = ir * ir;
    double ir3 = ir2 * ir;
    double ir4 = ir3 * ir;

    if ((order != 2) && (order != 4)) return NAN;  //must be order 2 or 4

    reduced_mass = AMU2KG * molecule_ptr->mass * pair_ptr->molecule->mass /
                   (molecule_ptr->mass + pair_ptr->molecule->mass);

    if (system->cdvdw_sig_repulsion) {
        dE = -6.0 * pair_ptr->sigrep * (2.0 * term12 - term6) * ir;
        d2E = 6.0 * pair_ptr->sigrep * (26.0 * term12 - 7.0 * term6) * ir2;
    } else {
        dE = -24.0 * pair_ptr->epsilon * (2.0 * term12 - term6) * ir;
        d2E = 24.0 * pair_ptr->epsilon * (26.0 * term12 - 7.0 * term6) * ir2;
    }

    //2nd order correction
    corr = M2A2 *
           (HBAR2 / (24.0 * KB * system->temperature * reduced_mass)) *
           (d2E + 2.0 * dE / pair_ptr->rimg);

    if (order >= 4) {
        if (system->cdvdw_sig_repulsion) {
            d3E = -336.0 * pair_ptr->sigrep * (6.0 * term12 - term6) * ir3;
            d4E = 3024.0 * pair_ptr->sigrep * (10.0 * term12 - term6) * ir4;
        } else {
            d3E = -1344.0 * pair_ptr->epsilon * (6.0 * term12 - term6) * ir3;
            d4E = 12096.0 * pair_ptr->epsilon * (10.0 * term12 - term6) * ir4;
        }

        //4th order corection
        corr += M2A4 *
                (HBAR4 / (1152.0 * KB2 * system->temperature * system->temperature * reduced_mass * reduced_mass)) *
                (15.0 * dE * ir3 + 4.0 * d3E * ir + d4E);
    }

    return corr;
}

double lj_lrc_corr(system_t *system, atom_t *atom_ptr, pair_t *pair_ptr, double cutoff) {
    double sig_cut, sig3, sig_cut3, sig_cut9;

    /* include the long-range correction */ /* I'm  not sure that I'm handling spectre pairs correctly */
    /* we can't use rd_excluded flag, since that disqualifies inter-molecular, but that DOES contribute to LRC */
    /* ALL OF THESE MUST BE TRUE TO PERFORM LRC CALCULATION */
    if ((pair_ptr->epsilon != 0 && pair_ptr->sigma != 0) &&                          //if these are zero, then we won't waste our time
        !(atom_ptr->spectre && pair_ptr->atom->spectre) &&                           //i think we want to disqualify s-s pairs
        !(pair_ptr->frozen) &&                                                       //disqualify frozen pairs
        ((pair_ptr->lrc == 0.0) || pair_ptr->last_volume != system->pbc->volume)) {  //LRC only changes if the volume change

        pair_ptr->last_volume = system->pbc->volume;

        sig_cut = fabs(pair_ptr->sigma) / cutoff;
        sig3 = fabs(pair_ptr->sigma);
        sig3 *= sig3 * sig3;
        sig_cut3 = sig_cut * sig_cut * sig_cut;
        sig_cut9 = sig_cut3 * sig_cut3 * sig_cut3;

        if (system->cdvdw_sig_repulsion)
            return (4.0 / 9.0) * M_PI * pair_ptr->sigrep * sig3 * sig_cut9 / system->pbc->volume;
        else if (system->polarvdw)  //only repulsion term, if polarvdw is on
            return (16.0 / 9.0) * M_PI * pair_ptr->epsilon * sig3 * sig_cut9 / system->pbc->volume;
        else  //if polarvdw is off, do the usual thing
            return ((16.0 / 3.0) * M_PI * pair_ptr->epsilon * sig3) * ((1.0 / 3.0) * sig_cut9 - sig_cut3) / system->pbc->volume;
    } else
        return pair_ptr->lrc;  //use stored value
}

double lj_lrc_self(system_t *system, atom_t *atom_ptr, double cutoff) {
    double sig_cut, sig3, sig_cut3, sig_cut9;

    if (((atom_ptr->sigma != 0) && (atom_ptr->epsilon != 0)) &&  //non-zero parameters
        !(atom_ptr->frozen) &&                                   //not frozen
        !(atom_ptr->spectre)) {                                  //not spectre

        sig_cut = fabs(atom_ptr->sigma) / cutoff;
        sig3 = fabs(atom_ptr->sigma);
        sig3 *= sig3 * sig3;
        sig_cut3 = sig_cut * sig_cut * sig_cut;
        sig_cut9 = sig_cut3 * sig_cut3 * sig_cut3;

        if (system->cdvdw_sig_repulsion)
            return (1.0 / 3.0) * M_PI * HBAR / KB * au2invseconds * atom_ptr->omega * atom_ptr->polarizability * atom_ptr->polarizability / sig3 * sig_cut9 / system->pbc->volume;
        else if (system->polarvdw)  //only repulsion term, if polarvdw is on
            return (16.0 / 9.0) * M_PI * atom_ptr->epsilon * sig3 * sig_cut9 / system->pbc->volume;
        else  //if polarvdw is off, do the usual thing
            return ((16.0 / 3.0) * M_PI * atom_ptr->epsilon * sig3) * ((1.0 / 3.0) * sig_cut9 - sig_cut3) / system->pbc->volume;
    }

    return 0;
}

double rd_crystal_self(system_t *system, atom_t *aptr, double cutoff) {
    double curr_pot, term12, term6;
    double sigma_over_r6, sigma_over_r12, sigma_over_r, r;
    int i[3], p, q;
    double a[3];
    curr_pot = 0;

    if (aptr->sigma == 0 && aptr->epsilon == 0) return 0;  //skip if no LJ interaction

    sigma_over_r6 = 0;
    sigma_over_r12 = 0;  //need to init these guys

    for (i[0] = -(system->rd_crystal_order - 1); i[0] <= system->rd_crystal_order - 1; i[0]++)
        for (i[1] = -(system->rd_crystal_order - 1); i[1] <= system->rd_crystal_order - 1; i[1]++)
            for (i[2] = -(system->rd_crystal_order - 1); i[2] <= system->rd_crystal_order - 1; i[2]++) {
                if (!i[0] && !i[1] && !i[2]) continue;  //no (0,0,0)
                //calculate pair separation (atom with it's image)
                for (p = 0; p < 3; p++) {
                    a[p] = 0;
                    for (q = 0; q < 3; q++)
                        a[p] += system->pbc->basis[q][p] * i[q];
                }
                r = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

                if (r > cutoff) continue;  //too far away! will be included in LRC if enabled
                sigma_over_r = fabs(aptr->sigma) / r;
                sigma_over_r6 += 0.5 * pow(sigma_over_r, 6);  //multiply by 0.5 to get counting correct
                sigma_over_r12 += 0.5 * pow(sigma_over_r, 12);
            }

    if (system->spectre) {
        term6 = 0;
        curr_pot = term12 = sigma_over_r12;
    } else {
        if (system->polarvdw)
            term6 = 0;  //vdw calc'd by vdw.c
        else
            term6 = sigma_over_r6;

        if (aptr->sigma < 0.0)
            term12 = 0;  //attractive only
        else
            term12 = sigma_over_r12;

        if (system->cdvdw_sig_repulsion)
            curr_pot = 0.75 * HBAR / KB * au2invseconds * aptr->omega * aptr->polarizability * aptr->polarizability /
                       pow(aptr->sigma, 6) * term12;  //C6*sig^6/r^12
        else if (system->polarvdw)
            curr_pot = 4.0 * aptr->epsilon * term12;
        else
            curr_pot = 4.0 * aptr->epsilon * (term12 - term6);
    }
    return curr_pot;
}

/* Lennard-Jones repulsion/dispersion */
double lj(system_t *system) {
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr;
    double sigma_over_r, term12, term6, sigma_over_r6, sigma_over_r12, r;  // , sigma6;   (unused variable)
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
                    if (system->rd_lrc) pair_ptr->lrc = lj_lrc_corr(system, atom_ptr, pair_ptr, cutoff);

                    // to include a contribution, we require
                    if ((pair_ptr->rimg - SMALL_dR < cutoff) &&            //inside cutoff?
                        (!pair_ptr->rd_excluded || system->rd_crystal) &&  //either not excluded OR rd_crystal is ON
                        !pair_ptr->frozen) {                               //not frozen

                        //loop over unit cells
                        if (system->rd_crystal) {
                            sigma_over_r6 = 0;
                            sigma_over_r12 = 0;
                            for (i[0] = -(system->rd_crystal_order - 1); i[0] <= system->rd_crystal_order - 1; i[0]++)
                                for (i[1] = -(system->rd_crystal_order - 1); i[1] <= system->rd_crystal_order - 1; i[1]++)
                                    for (i[2] = -(system->rd_crystal_order - 1); i[2] <= system->rd_crystal_order - 1; i[2]++) {
                                        if (!i[0] && !i[1] && !i[2] && pair_ptr->rd_excluded) continue;  //no i=j=k=0 for excluded pairs (intra-molecular)
                                        //calculate pair separation (atom with it's image)
                                        for (p = 0; p < 3; p++) {
                                            a[p] = 0;
                                            for (q = 0; q < 3; q++)
                                                a[p] += system->pbc->basis[q][p] * i[q];
                                            a[p] += atom_ptr->pos[p] - pair_ptr->atom->pos[p];
                                        }
                                        r = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

                                        if (r > cutoff) continue;
                                        sigma_over_r = fabs(pair_ptr->sigma) / r;
                                        sigma_over_r6 += pow(sigma_over_r, 6);
                                        sigma_over_r12 += pow(sigma_over_r, 12);
                                    }
                        } else {  //otherwise, calculate as normal
                            sigma_over_r = fabs(pair_ptr->sigma) / pair_ptr->rimg;
                            sigma_over_r6 = sigma_over_r * sigma_over_r * sigma_over_r;
                            sigma_over_r6 *= sigma_over_r6;
                            sigma_over_r12 = sigma_over_r6 * sigma_over_r6;
                        }

                        /* the LJ potential */
                        if (system->spectre) {
                            term6 = 0;
                            term12 = sigma_over_r12;
                            potential_classical = term12;
                        } else {
                            if (system->polarvdw)
                                term6 = 0;  //vdw calc'd by vdw.c
                            else
                                term6 = sigma_over_r6;

                            if (pair_ptr->attractive_only)
                                term12 = 0;
                            else
                                term12 = sigma_over_r12;

                            if (system->cdvdw_sig_repulsion)
                                potential_classical = pair_ptr->sigrep * term12;  //C6*sig^6/r^12
                            else
                                potential_classical = 4.0 * pair_ptr->epsilon * (term12 - term6);
                        }

                        pair_ptr->rd_energy += potential_classical;

                        if (system->feynman_hibbs)
                            pair_ptr->rd_energy += lj_fh_corr(system, molecule_ptr, pair_ptr, system->feynman_hibbs_order, term12, term6);

                    }  // if qualified contributions

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
                potential += rd_crystal_self(system, atom_ptr, cutoff);

    /* calculate self LRC interaction */
    if (system->rd_lrc)
        for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
            for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
                potential += lj_lrc_self(system, atom_ptr, cutoff);

    return (potential);
}

/* same as above, but no periodic boundary conditions */
double lj_nopbc(system_t *system) {
    molecule_t *molecules = system->molecules;
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr;
    double sigma_over_r, term12, term6, sigma_over_r6;
    double potential;

    for (molecule_ptr = molecules, potential = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
                /* make sure we're not excluded or beyond the cutoff */
                if (!pair_ptr->rd_excluded) {
                    sigma_over_r = fabs(pair_ptr->sigma) / pair_ptr->r;
                    sigma_over_r6 = sigma_over_r * sigma_over_r * sigma_over_r;
                    sigma_over_r6 *= sigma_over_r6;

                    if (system->polarvdw)
                        term6 = 0;
                    else
                        term6 = sigma_over_r6;

                    if (pair_ptr->attractive_only)
                        term12 = 0;
                    else
                        term12 = sigma_over_r6 * sigma_over_r6;

                    if (system->cdvdw_sig_repulsion)
                        potential += pair_ptr->sigrep * term12;  //C6*sig^6/r^12
                    else
                        potential += 4.0 * pair_ptr->epsilon * (term12 - term6);
                }

            } /* pair */
        }     /* atom */
    }         /* molecule */

    return (potential);
}

#ifdef DEBUG
void test_lj(system_t *system) {
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
