/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

//set them to alpha*E_static
void init_dipoles(system_t *system) {
    int i, p;
    atom_t **aa = system->atom_array;

    for (i = 0; i < system->natoms; i++) {
        for (p = 0; p < 3; p++) {
            aa[i]->mu[p] = aa[i]->polarizability * (aa[i]->ef_static[p] + aa[i]->ef_static_self[p]);
            // should improve convergence since mu's typically grow as induced fields are added in
            if (!system->polar_sor && !system->polar_esor) aa[i]->mu[p] *= system->polar_gamma;
        }
    }
    return;
}

void contract_dipoles(system_t *system, int *ranked_array) {
    int i, j, ii, jj, p, index;
    atom_t **aa = system->atom_array;

    for (i = 0; i < system->natoms; i++) {
        index = ranked_array[i];  //do them in the order of the ranked index
        ii = index * 3;
        if (aa[index]->polarizability == 0) {  //if not polar
            //aa[index]->ef_induced[p] is already 0
            aa[index]->new_mu[0] = aa[index]->new_mu[1] = aa[index]->new_mu[2] = 0;  //might be redundant?
            aa[index]->mu[0] = aa[index]->mu[1] = aa[index]->mu[2] = 0;              //might be redundant?
            continue;
        }
        for (j = 0; j < system->natoms; j++) {
            jj = j * 3;
            if (index != j)
                for (p = 0; p < 3; p++)
                    aa[index]->ef_induced[p] -= dddotprod((system->A_matrix[ii + p] + jj), aa[j]->mu);
        } /* end j */

        /* dipole is the sum of the static and induced parts */
        for (p = 0; p < 3; p++) {
            aa[index]->new_mu[p] = aa[index]->polarizability * (aa[index]->ef_static[p] + aa[index]->ef_static_self[p] + aa[index]->ef_induced[p]);

            /* Gauss-Seidel */
            if (system->polar_gs || system->polar_gs_ranked)
                aa[index]->mu[p] = aa[index]->new_mu[p];
        }

    } /* end matrix multiply */

    return;
}

void calc_dipole_rrms(system_t *system) {
    int i, p;
    double carry;
    atom_t **aa = system->atom_array;

    /* get the dipole RRMS */
    for (i = 0; i < system->natoms; i++) {
        /* mean square difference */
        aa[i]->dipole_rrms = 0;
        for (p = 0; p < 3; p++) {
            carry = aa[i]->new_mu[p] - aa[i]->old_mu[p];
            aa[i]->dipole_rrms += carry * carry;
        }

        /* normalize */
        aa[i]->dipole_rrms /= dddotprod(aa[i]->new_mu, aa[i]->new_mu);
        aa[i]->dipole_rrms = sqrt(aa[i]->dipole_rrms);
        if (!isfinite(aa[i]->dipole_rrms)) aa[i]->dipole_rrms = 0;
    }

#ifdef DEBUG
/*
double totalrms = 0;
for ( i=0; i< system->natoms; i++ ) {
		totalrms +=  aa[i]->dipole_rrms;
}
fprintf(stderr,"TOTAL DIPOLE RRMS %lf\n", totalrms);
*/
#endif

    return;
}

int are_we_done_yet(system_t *system, int iteration_counter) {
    int i, p;
    int N = system->natoms;
    atom_t **aa = system->atom_array;
    double allowed_sqerr, error;

    if (system->polar_precision == 0.0) { /* ... by fixed iteration ... */
        if (iteration_counter != system->polar_max_iter)
            return 1;
    }

    else { /* ... or by dipole precision */
        allowed_sqerr = system->polar_precision * system->polar_precision * DEBYE2SKA * DEBYE2SKA;
        for (i = 0; i < N; i++) {  //check the change in each dipole component
            for (p = 0; p < 3; p++) {
                error = aa[i]->new_mu[p] - aa[i]->old_mu[p];
                if (error * error > allowed_sqerr)
                    return 1;  //we broke tolerance
            }
        }
    }

    return 0;
}

void palmo_contraction(system_t *system, int *ranked_array) {
    int i, j, ii, jj, index, p;
    int N = system->natoms;
    atom_t **aa = system->atom_array;

    /* calculate change in induced field due to this iteration */
    for (i = 0; i < N; i++) {
        index = ranked_array[i];
        ii = index * 3;

        for (p = 0; p < 3; p++)
            aa[index]->ef_induced_change[p] = -aa[index]->ef_induced[p];

        for (j = 0; j < N; j++) {
            jj = j * 3;
            if (index != j)
                for (p = 0; p < 3; p++)
                    aa[index]->ef_induced_change[p] -= dddotprod(system->A_matrix[ii + p] + jj, aa[j]->mu);
        }
    }

    return;
}

void update_ranking(system_t *system, int *ranked_array) {
    int i, j, sorted, tmp;
    int N = system->natoms;
    atom_t **aa = system->atom_array;

    /* rank the dipoles by bubble sort */
    if (system->polar_gs_ranked) {
        for (i = 0; i < N; i++) {
            for (j = 0, sorted = 1; j < (N - 1); j++) {
                if (aa[ranked_array[j]]->rank_metric < aa[ranked_array[j + 1]]->rank_metric) {
                    sorted = 0;
                    tmp = ranked_array[j];
                    ranked_array[j] = ranked_array[j + 1];
                    ranked_array[j + 1] = tmp;
                }
            }
            if (sorted) break;
        }
    }

    return;
}

/* iterative solver of the dipole field tensor */
/* returns the number of iterations required */
int thole_iterative(system_t *system) {
    int i, N, p;
    int iteration_counter, keep_iterating;
    atom_t **aa;  //atom array
    int *ranked_array;

    aa = system->atom_array;
    N = system->natoms;

    /* array for ranking */
    ranked_array = calloc(N, sizeof(int));
    memnullcheck(ranked_array, N * sizeof(int), __LINE__ - 1, __FILE__);
    for (i = 0; i < N; i++) ranked_array[i] = i;

    //set all dipoles to alpha*E_static * polar_gamma
    init_dipoles(system);

    /* if ZODID is enabled, then stop here and just return the alpha*E dipoles */
    if (system->polar_zodid) {
        free(ranked_array);
        return (0);
    }

    /* iterative solver of the dipole field equations */
    keep_iterating = 1;
    iteration_counter = 0;
    while (keep_iterating) {
        iteration_counter++;

        /* divergence detection */
        /* if we fail to converge, then return dipoles as alpha*E */
        if (iteration_counter >= MAX_ITERATION_COUNT && system->polar_precision) {
            for (i = 0; i < N; i++)
                for (p = 0; p < 3; p++) {
                    aa[i]->mu[p] = aa[i]->polarizability * (aa[i]->ef_static[p] + aa[i]->ef_static_self[p]);
                    aa[i]->ef_induced_change[p] = 0.0;  //so we don't break palmo
                }
            //set convergence failure flag
            system->iter_success = 1;

            free(ranked_array);
            return (iteration_counter);
        }

        //zero out induced e-field
        for (i = 0; i < system->natoms; i++)
            for (p = 0; p < 3; p++)
                aa[i]->ef_induced[p] = 0;

        //save the current dipole information if we want to calculate precision (or if needed for relaxation)
        if (system->polar_rrms || system->polar_precision > 0 || system->polar_sor || system->polar_esor) {
            for (i = 0; i < N; i++)
                for (p = 0; p < 3; p++)
                    aa[i]->old_mu[p] = aa[i]->mu[p];
        }

        // contract the dipoles with the field tensor (gauss-seidel/gs-ranked optional)
        contract_dipoles(system, ranked_array);

        if (system->polar_rrms || system->polar_precision > 0)
            calc_dipole_rrms(system);

        /* determine if we are done... */
        keep_iterating = are_we_done_yet(system, iteration_counter);

        // if we would be finished, we want to contract once more to get the next induced field for palmo
        if (system->polar_palmo && !keep_iterating)
            palmo_contraction(system, ranked_array);

        //new gs_ranking if needed
        if (system->polar_gs_ranked && keep_iterating)
            update_ranking(system, ranked_array);

        /* save the dipoles for the next pass */
        for (i = 0; i < N; i++) {
            for (p = 0; p < 3; p++) {
                /* allow for different successive over-relaxation schemes */
                if (system->polar_sor)
                    aa[i]->mu[p] = system->polar_gamma * aa[i]->new_mu[p] + (1.0 - system->polar_gamma) * aa[i]->old_mu[p];
                else if (system->polar_esor)
                    aa[i]->mu[p] = (1.0 - exp(-system->polar_gamma * iteration_counter)) * aa[i]->new_mu[p] + exp(-system->polar_gamma * iteration_counter) * aa[i]->old_mu[p];
                else
                    aa[i]->mu[p] = aa[i]->new_mu[p];
            }
        }

    }  //end iterate
    free(ranked_array);

    /* return the iteration count */
    return (iteration_counter);
}
