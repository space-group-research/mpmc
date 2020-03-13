/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

void rebuild_arrays(system_t *system) {
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    int n = system->natoms;

    free(system->atom_array);
    free(system->molecule_array);

    //allocate the arrays
    system->molecule_array = malloc(n * sizeof(molecule_t *));
    memnullcheck(system->molecule_array, n * sizeof(molecule_t *), __LINE__ - 1, __FILE__);
    system->atom_array = malloc(n * sizeof(atom_t *));
    memnullcheck(system->atom_array, n * sizeof(atom_t *), __LINE__ - 1, __FILE__);

    n = 0;
    //build the arrays
    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            system->molecule_array[n] = molecule_ptr;
            system->atom_array[n] = atom_ptr;
            n++;
        }
    }
    system->natoms = n;

    return;
}

/* flag all pairs to have their energy calculated */
/* needs to be called at simulation start, or can */
/* be called to periodically keep the total energy */
/* from drifting */
void flag_all_pairs(system_t *system) {
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr;

    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
            for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next)
                pair_ptr->recalculate_energy = 1;
}

/* set the exclusions and LJ mixing for relevant pairs */
void pair_exclusions(system_t *system, molecule_t *molecule_i, molecule_t *molecule_j, atom_t *atom_i, atom_t *atom_j, pair_t *pair_ptr) {
    double si3, sj3, si6, sj6;
    double repul1, repul2, repulmix;

    /* recalculate exclusions */
    if ((molecule_i == molecule_j) && !system->gwp) { /* if both on same molecule, exclude all interactions */

        pair_ptr->rd_excluded = 1;
        pair_ptr->es_excluded = 1;

    } else {
        /* exclude null repulsion/dispersion interations */
        if (((atom_i->epsilon == 0.0) || (atom_i->sigma == 0.0) || (atom_j->epsilon == 0.0) || (atom_j->sigma == 0.0)) && (atom_i->c6 == 0.0 && atom_i->c8 == 0.0 && atom_i->c10 == 0.0 && atom_j->c6 == 0.0 && atom_j->c8 == 0.0 && atom_j->c10 == 0.0))
            pair_ptr->rd_excluded = 1;
        else
            pair_ptr->rd_excluded = 0;

        /* exclude null electrostatic interactions */
        if ((atom_i->charge == 0.0) || (atom_j->charge == 0.0))
            pair_ptr->es_excluded = 1;
        else
            pair_ptr->es_excluded = 0;
    }

    /* get the frozen interactions */
    pair_ptr->frozen = atom_i->frozen && atom_j->frozen;

    /* get the mixed LJ parameters */
    if (!system->sg) {
        if (system->waldmanhagler && !system->cdvdw_sig_repulsion) {  //wh mixing rule
            si3 = atom_i->sigma;
            si3 *= si3 * si3;
            si6 = si3 * si3;

            sj3 = atom_j->sigma;
            sj3 *= sj3 * sj3;
            sj6 = sj3 * sj3;

            if ((atom_i->sigma < 0.0) || (atom_j->sigma < 0.0)) {
                pair_ptr->attractive_only = 1;
                pair_ptr->sigma = pow(0.5 * (si6 + sj6), 1. / 6.);
            } else if ((atom_i->sigma == 0 || atom_j->sigma == 0)) {
                pair_ptr->sigma = 0;
                pair_ptr->epsilon = sqrt(atom_i->epsilon * atom_j->epsilon);  //can't use sigma weights -> div by 0
            } else {
                pair_ptr->sigma = pow(0.5 * (si6 + sj6), 1. / 6.);
                pair_ptr->epsilon = sqrt(atom_i->epsilon * atom_j->epsilon) * 2.0 * si3 * sj3 / (si6 + sj6);
            }
        } else if (system->halgren_mixing) {  //halgren mixing rules
            if (atom_i->sigma > 0.0 && atom_j->sigma > 0.0) {
                pair_ptr->sigma = (atom_i->sigma * atom_i->sigma * atom_i->sigma + atom_j->sigma * atom_j->sigma * atom_j->sigma) / (atom_i->sigma * atom_i->sigma + atom_j->sigma * atom_j->sigma);
            } else {
                pair_ptr->sigma = 0;
            }
            if (atom_i->epsilon > 0.0 && atom_j->epsilon > 0.0) {
                pair_ptr->epsilon = 4 * atom_i->epsilon * atom_j->epsilon / pow(sqrt(atom_i->epsilon) + sqrt(atom_j->epsilon), 2);
            } else {
                pair_ptr->epsilon = 0;
            }
        } else if (system->cdvdw_9th_repulsion) {  //9th power mixing for repulsion
            si3 = atom_i->sigma;
            si3 *= si3 * si3;
            si6 = si3 * si3;
            sj3 = atom_j->sigma;
            sj3 *= sj3 * sj3;
            sj6 = sj3 * sj3;
            repul1 = 4.0 * si6 * si6 * atom_i->epsilon;
            repul2 = 4.0 * sj6 * sj6 * atom_j->epsilon;
            repulmix = pow(0.5 * (pow(repul1, 1. / 9.) + pow(repul2, 1. / 9.)), 9);
            pair_ptr->sigma = 1.0;
            pair_ptr->epsilon = repulmix / 4.0;
        } else if (system->cdvdw_sig_repulsion) {  //sigma repulsion for coupled-dipole vdw
            si3 = atom_i->sigma;
            si3 *= si3 * si3;
            si6 = si3 * si3;
            sj3 = atom_j->sigma;
            sj3 *= sj3 * sj3;
            sj6 = sj3 * sj3;
            pair_ptr->sigma = pow(0.5 * (si6 + sj6), 1. / 6.);
            pair_ptr->sigrep = 1.5 * HBAR / KB * au2invseconds * atom_i->omega * atom_j->omega *
                               atom_i->polarizability * atom_j->polarizability / (atom_i->omega + atom_j->omega) / pow(pair_ptr->sigma, 6);
        } else if ((system->polarvdw && system->cdvdw_exp_repulsion)) {  // mix for buckingham repulsion
            // sigma == C, epsilon == rho
            // U = C exp(-R/(2*rho))
            pair_ptr->sigma = pow(pow(atom_i->sigma, atom_i->epsilon) * pow(atom_j->sigma, atom_j->epsilon), 1.0 / ((atom_i->epsilon + atom_j->epsilon)));
            pair_ptr->epsilon = 0.5 * (atom_i->epsilon + atom_j->epsilon);
        } else if (system->disp_expansion) {
            // https://journals.aps.org/pra/pdf/10.1103/PhysRevA.5.1708
            pair_ptr->sigma = 0.5 * (atom_i->sigma + atom_j->sigma);
            pair_ptr->epsilon = 2.0 * atom_i->epsilon * atom_j->epsilon / (atom_i->epsilon + atom_j->epsilon);

            // get mixed dispersion coefficients
            pair_ptr->c6 = sqrt(atom_i->c6 * atom_j->c6) * 0.021958709 / (3.166811429 * 0.000001);   // Convert H*Bohr^6 to K*Angstrom^6, etc
            pair_ptr->c8 = sqrt(atom_i->c8 * atom_j->c8) * 0.0061490647 / (3.166811429 * 0.000001);  // Dispersion coeffs should be inputed in a.u.

            if (system->extrapolate_disp_coeffs && pair_ptr->c6 != 0.0 && pair_ptr->c8 != 0.0)
                pair_ptr->c10 = 49.0 / 40.0 * pair_ptr->c8 * pair_ptr->c8 / pair_ptr->c6;
            else if (system->extrapolate_disp_coeffs)
                pair_ptr->c10 = 0.0;  //either c6 or c8 is zero so lets set c10 to zero too
            else
                pair_ptr->c10 = sqrt(atom_i->c10 * atom_j->c10) * 0.0017219135 / (3.166811429 * 0.000001);
        } else if (system->c6_mixing) {
            pair_ptr->sigma = 0.5 * (atom_i->sigma + atom_j->sigma);
            if (pair_ptr->sigma != 0.0)
                pair_ptr->epsilon = 64.0 * sqrt(atom_i->epsilon * atom_j->epsilon) * pow(atom_i->sigma, 3.0) * pow(atom_j->sigma, 3.0) / pow(atom_i->sigma + atom_j->sigma, 6.0);
            else
                pair_ptr->epsilon = 0.0;
        } else { /* lorentz-berthelot */
            if ((atom_i->sigma < 0.0) || (atom_j->sigma < 0.0)) {
                pair_ptr->attractive_only = 1;
                pair_ptr->sigma = 0.5 * (fabs(atom_i->sigma) + fabs(atom_j->sigma));
            } else if ((atom_i->sigma == 0 || atom_j->sigma == 0)) {
                pair_ptr->sigma = 0;
                pair_ptr->epsilon = sqrt(atom_i->epsilon * atom_j->epsilon);
            } else {
                pair_ptr->sigma = 0.5 * (atom_i->sigma + atom_j->sigma);
                pair_ptr->epsilon = sqrt(atom_i->epsilon * atom_j->epsilon);
            }
        }
    } /*!sg*/

    /* ensure that no ES calc for S-S pairs, and ONLY ES for S-everythingelse */
    if (system->spectre) {
        if (atom_i->spectre && atom_j->spectre) { /* case for both spectre */

            pair_ptr->rd_excluded = 0;
            pair_ptr->es_excluded = 1;

        } else if (atom_i->spectre || atom_j->spectre) { /* case for one spectre */

            pair_ptr->rd_excluded = 1;
            pair_ptr->es_excluded = 0;
        }
    }
}

/* perform the modulo minimum image for displacements */
void minimum_image(system_t *system, atom_t *atom_i, atom_t *atom_j, pair_t *pair_ptr) {
    int p, q;
    double img[3];
    double d[3], r, r2;
    double di[3], ri, ri2;
    // double rnd;  (unused variable)

    /* get the real displacement */
    pair_ptr->recalculate_energy = 0; /* reset the recalculate flag */
    for (p = 0; p < 3; p++) {
        d[p] = atom_i->pos[p] - atom_j->pos[p];

        /* this pair changed and so it will have it's energy recalculated */
        if (d[p] != pair_ptr->d_prev[p]) {
            pair_ptr->recalculate_energy = 1;
            pair_ptr->d_prev[p] = d[p]; /* reset */
        }
    }

    //relative position didn't change. nothing to do here.
    if (pair_ptr->recalculate_energy == 0) return;

    for (p = 0; p < 3; p++) {
        for (q = 0, img[p] = 0; q < 3; q++) {
            img[p] += system->pbc->reciprocal_basis[q][p] * d[q];
        }
        img[p] = rint(img[p]);
    }

    /* matrix multiply to project back into our basis */
    for (p = 0; p < 3; p++)
        for (q = 0, di[p] = 0; q < 3; q++)
            di[p] += system->pbc->basis[q][p] * img[q];

    /* now correct the displacement */
    for (p = 0; p < 3; p++)
        di[p] = d[p] - di[p];

    /* pythagorean terms */
    for (p = 0, r2 = 0, ri2 = 0; p < 3; p++) {
        r2 += d[p] * d[p];
        ri2 += di[p] * di[p];
    }
    r = sqrt(r2);
    ri = sqrt(ri2);

    /* store the results for this pair */
    pair_ptr->r = r;

    if (isnan(ri) != 0) {
        pair_ptr->rimg = r;
        for (p = 0; p < 3; p++)
            pair_ptr->dimg[p] = d[p];
    } else {
        pair_ptr->rimg = ri;
        for (p = 0; p < 3; p++)
            pair_ptr->dimg[p] = di[p];
    }

    return;
}

/* update everything necessary to describe the complete pairwise system */
void pairs(system_t *system) {
    int i, j, n;
    // molecule_t *molecule_ptr;     (unused variable)
    // atom_t *atom_ptr;    (unused variable)
    pair_t *pair_ptr;
    molecule_t **molecule_array;
    atom_t **atom_array;
    /* needed for GS ranking metric */
    // int p;   (unused variable)
    // double r;    (unused variable)
    double rmin;

    // get array of atom ptrs
    rebuild_arrays(system);
    atom_array = system->atom_array;
    molecule_array = system->molecule_array;
    n = system->natoms;

    /* loop over all atoms and pair */
    for (i = 0; i < (n - 1); i++) {
        for (j = (i + 1), pair_ptr = atom_array[i]->pairs; j < n; j++, pair_ptr = pair_ptr->next) {
            /* set the link */
            pair_ptr->atom = atom_array[j];
            pair_ptr->molecule = molecule_array[j];

            //this is dangerous and has already been responsible for numerous bugs, most recently
            //in UVT runs. after and insert/remove move there is no guarantee that pair_ptr->rd_excluded is properly set
            //if ( !pair_ptr->frozen && !(pair_ptr->rd_excluded && pair_ptr->es_excluded) )
            pair_exclusions(system, molecule_array[i], molecule_array[j], atom_array[i], atom_array[j], pair_ptr);

            /* recalc min image */
            if (!pair_ptr->frozen || system->polarization)  //need induced-induced interaction for frozen atoms
                minimum_image(system, atom_array[i], atom_array[j], pair_ptr);

        } /* for j */
    }     /* for i */

    /* update the com of each molecule */
    update_com(system->molecules);

    /* store wrapped coords */
    wrapall(system->molecules, system->pbc);

    /* rank metric */
    if (system->polar_iterative && system->polar_gs_ranked) {
        /* determine smallest polarizable separation */
        rmin = MAXVALUE;
        for (i = 0; i < n; i++) {
            if (atom_array[i]->polarizability == 0.0) continue;
            for (pair_ptr = atom_array[i]->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
                if (pair_ptr->atom->polarizability == 0.0) continue;
                if (pair_ptr->rimg < rmin) rmin = pair_ptr->rimg;
            }
        }
        //calculate rank shits
        for (i = 0; i < n; i++)
            atom_array[i]->rank_metric = 0;
        for (i = 0; i < n; i++) {
            if (atom_array[i]->polarizability == 0.0) continue;
            for (pair_ptr = atom_array[i]->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
                if (pair_ptr->atom->polarizability == 0.0) continue;
                if (pair_ptr->r <= rmin * 1.5) {
                    atom_array[i]->rank_metric += 1.0;
                    pair_ptr->atom->rank_metric += 1.0;
                }
            }
        }
    }
}

/* molecular center of mass */
void update_com(molecule_t *molecules) {
    int i;
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;

    for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (i = 0; i < 3; i++)
            molecule_ptr->com[i] = 0;

        if (!(molecule_ptr->spectre || molecule_ptr->target)) {
            for (atom_ptr = molecule_ptr->atoms, molecule_ptr->mass = 0; atom_ptr; atom_ptr = atom_ptr->next) {
                molecule_ptr->mass += atom_ptr->mass;

                for (i = 0; i < 3; i++)
                    molecule_ptr->com[i] += atom_ptr->mass * atom_ptr->pos[i];
            }

            for (i = 0; i < 3; i++)
                molecule_ptr->com[i] /= molecule_ptr->mass;
        }
    }
}

/* add new pairs for when a new molecule is created */
void update_pairs_insert(system_t *system) {
    int i, n;
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr;

    /* count the number of atoms per molecule */
    for (atom_ptr = system->checkpoint->molecule_altered->atoms, n = 0; atom_ptr; atom_ptr = atom_ptr->next, n++)
        ;

    /* add n number of pairs to altered and all molecules ahead of it in the list */
    for (molecule_ptr = system->molecules; molecule_ptr != system->checkpoint->tail; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            /* go to the end of the pair list */
            if (atom_ptr->pairs) {
                /* go to the end of the pair list */
                for (pair_ptr = atom_ptr->pairs; pair_ptr->next; pair_ptr = pair_ptr->next)
                    ;

                /* tag on the extra pairs */
                for (i = 0; i < n; i++) {
                    pair_ptr->next = calloc(1, sizeof(pair_t));
                    memnullcheck(pair_ptr->next, sizeof(pair_t), __LINE__ - 1, __FILE__);
                    pair_ptr = pair_ptr->next;
                }

            } else {
                /* needs a new list */
                atom_ptr->pairs = calloc(1, sizeof(pair_t));
                memnullcheck(atom_ptr->pairs, sizeof(pair_t), __LINE__ - 1, __FILE__);
                pair_ptr = atom_ptr->pairs;
                for (i = 0; i < (n - 1); i++) {
                    pair_ptr->next = calloc(1, sizeof(pair_t));
                    memnullcheck(pair_ptr->next, sizeof(pair_t), __LINE__ - 1, __FILE__);
                    pair_ptr = pair_ptr->next;
                }
            }

        } /* for atom */
    }     /* for molecule */
}

/* remove pairs when a molecule is deleted */
void update_pairs_remove(system_t *system) {
    int i, n, m;
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr, **pair_array;

    /* count the number of atoms per molecule */
    for (atom_ptr = system->checkpoint->molecule_backup->atoms, n = 0; atom_ptr; atom_ptr = atom_ptr->next, n++)
        ;

    /* remove n number of pairs for all molecules ahead of the removal point */
    pair_array = calloc(1, sizeof(pair_t *));
    memnullcheck(pair_array, sizeof(pair_t *), __LINE__ - 1, __FILE__);
    for (molecule_ptr = system->molecules; molecule_ptr != system->checkpoint->tail; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            /* build the pair pointer array */
            for (pair_ptr = atom_ptr->pairs, m = 0; pair_ptr; pair_ptr = pair_ptr->next, m++) {
                pair_array = realloc(pair_array, sizeof(pair_t *) * (m + 1));
                memnullcheck(pair_array, sizeof(pair_t *) * (m + 1), __LINE__ - 1, __FILE__);
                pair_array[m] = pair_ptr;
            }

            for (i = (m - n); i < m; i++)
                free(pair_array[i]);

            /* handle the end of the list */
            if ((m - n) > 0)
                pair_array[(m - n - 1)]->next = NULL;
            else
                atom_ptr->pairs = NULL;
        }
    }

    /* free our temporary array */
    free(pair_array);
}

/* if an insert move is rejected, remove the pairs that were previously added */
void unupdate_pairs_insert(system_t *system) {
    int i, n, m;
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr, **pair_array;

    /* count the number of atoms per molecule */
    for (atom_ptr = system->checkpoint->molecule_altered->atoms, n = 0; atom_ptr; atom_ptr = atom_ptr->next)
        ++n;

    /* remove n number of pairs for all molecules ahead of the removal point */
    pair_array = calloc(1, sizeof(pair_t *));
    memnullcheck(pair_array, sizeof(pair_t *), __LINE__ - 1, __FILE__);
    for (molecule_ptr = system->molecules; molecule_ptr != system->checkpoint->tail; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            /* build the pair pointer array */
            for (pair_ptr = atom_ptr->pairs, m = 0; pair_ptr; pair_ptr = pair_ptr->next, m++) {
                pair_array = realloc(pair_array, sizeof(pair_t *) * (m + 1));
                memnullcheck(pair_array, sizeof(pair_t *) * (m + 1), __LINE__ - 1, __FILE__);
                pair_array[m] = pair_ptr;
            }

            for (i = (m - n); i < m; i++)
                free(pair_array[i]);

            /* handle the end of the list */
            if ((m - n) > 0)
                pair_array[(m - n - 1)]->next = NULL;
            else
                atom_ptr->pairs = NULL;
        }
    }

    /* free our temporary array */
    free(pair_array);
}

/* if a remove is rejected, then add back the pairs that were previously deleted */
void unupdate_pairs_remove(system_t *system) {
    int i, n;
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr;

    /* count the number of atoms per molecule */
    for (atom_ptr = system->checkpoint->molecule_backup->atoms, n = 0; atom_ptr; atom_ptr = atom_ptr->next)
        ++n;

    /* add n number of pairs to altered and all molecules ahead of it in the list */
    for (molecule_ptr = system->molecules; molecule_ptr != system->checkpoint->molecule_backup; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            /* go to the end of the pair list */
            if (atom_ptr->pairs) {
                /* go to the end of the pair list */
                for (pair_ptr = atom_ptr->pairs; pair_ptr->next; pair_ptr = pair_ptr->next)
                    ;

                /* tag on the extra pairs */
                for (i = 0; i < n; i++) {
                    pair_ptr->next = calloc(1, sizeof(pair_t));
                    memnullcheck(pair_ptr->next, sizeof(pair_t), __LINE__ - 1, __FILE__);
                    pair_ptr = pair_ptr->next;
                }

            } else {
                /* needs a new list */
                atom_ptr->pairs = calloc(1, sizeof(pair_t));
                memnullcheck(atom_ptr->pairs, sizeof(pair_t), __LINE__ - 1, __FILE__);
                pair_ptr = atom_ptr->pairs;
                for (i = 0; i < (n - 1); i++) {
                    pair_ptr->next = calloc(1, sizeof(pair_t));
                    memnullcheck(pair_ptr->next, sizeof(pair_t), __LINE__ - 1, __FILE__);
                    pair_ptr = pair_ptr->next;
                }
            }

        } /* for atom */
    }     /* for molecule */
}

/* allocate the pair lists */
void setup_pairs(system_t *system) {
    int i, j, n;
    // molecule_t *molecule_ptr;  (unused variable)
    atom_t **atom_array;  //, *atom_ptr;   (unused variable)
    pair_t *pair_ptr, *prev_pair_ptr;

    system->natoms = countNatoms(system);

    //build atom and molecule arrays
    rebuild_arrays(system);
    atom_array = system->atom_array;
    n = system->natoms;

    /* setup the pairs, lower triangular */
    for (i = 0; i < (n - 1); i++) {
        atom_array[i]->pairs = calloc(1, sizeof(pair_t));
        memnullcheck(atom_array[i]->pairs, sizeof(pair_t), __LINE__ - 1, __FILE__);
        pair_ptr = atom_array[i]->pairs;
        prev_pair_ptr = pair_ptr;

        for (j = (i + 1); j < n; j++) {
            pair_ptr->next = calloc(1, sizeof(pair_t));
            memnullcheck(pair_ptr->next, sizeof(pair_t), __LINE__ - 1, __FILE__);
            prev_pair_ptr = pair_ptr;
            pair_ptr = pair_ptr->next;
        }

        prev_pair_ptr->next = NULL;
        free(pair_ptr);
    }
}

#ifdef DEBUG
void test_pairs(molecule_t *molecules) {
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr;

    for (molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
                if (!(pair_ptr->frozen || pair_ptr->rd_excluded || pair_ptr->es_excluded))
                    printf(
                        "DEBUG_PAIRS: atomid %d charge = %f, epsilon = %f, sigma = %f, r = %f, rimg = %f\n",
                        atom_ptr->id, pair_ptr->atom->charge, pair_ptr->epsilon,
                        pair_ptr->sigma, pair_ptr->r, pair_ptr->rimg);
                fflush(stdout);
            }
        }
    }
}
#endif
