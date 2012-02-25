/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>


/* iterative solver of the dipole field tensor */
/* returns the number of iterations required */
int thole_iterative(system_t *system) {

	int i, j, ii, jj, N, p, q, sorted;
	int iteration_counter, keep_iterating;
	double error;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr, **atom_array;
	int *ranked_array, ranked, tmp, index;


	/* generate an array of atom ptrs */
	for(molecule_ptr = system->molecules, N = 0, atom_array = NULL; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next, N++) {

			atom_array = realloc(atom_array, sizeof(atom_t *)*(N + 1));
			memnullcheck(atom_array,sizeof(atom_t *)*(N+1),80);
			atom_array[N] = atom_ptr;

		}
	}

	/* array for ranking */
	ranked_array = calloc(N, sizeof(int));
	memnullcheck(ranked_array,N*sizeof(int),81);
	for(i = 0; i < N; i++) ranked_array[i] = i;

	/* initialize things */
	for(i = 0; i < N; i++) {
		for(p = 0; p < 3; p++) {

			/* set the first guess to alpha*E */
			atom_array[i]->mu[p] = atom_array[i]->polarizability*atom_array[i]->ef_static[p];
			if(!system->polar_sor) atom_array[i]->mu[p] *= system->polar_gamma;

		}
	}

	/* if ZODID is enabled, then stop here and just return the alpha*E dipoles */
	if(system->polar_zodid) return(0);

	/* iterative solver of the dipole field equations */
	for(iteration_counter = 0, keep_iterating = 1; keep_iterating; iteration_counter++) {

		/* divergence detection */
		/* if we fail to converge, then return dipoles as alpha*E */
		if(iteration_counter >= MAX_ITERATION_COUNT) {

			for(i = 0; i < N; i++)
				for(p = 0; p < 3; p++)
					atom_array[i]->mu[p] = atom_array[i]->polarizability*atom_array[i]->ef_static[p];

			free(atom_array);
			return(iteration_counter);

		}


		/* save the current dipole set and clear the induced field vectors */
		for(i = 0; i < N; i++) {
			for(p = 0; p < 3; p++) {

				atom_array[i]->old_mu[p] = atom_array[i]->mu[p];
				atom_array[i]->ef_induced[p] = 0;
				atom_array[i]->ef_induced_change[p] = 0;

			}
		}


		/* contract the dipoles with the field tensor */
		for(i = 0; i < N; i++) {
			index = ranked_array[i];
			ii = index*3;

			for(j = 0; j < N; j++) {
				jj = j*3;

				if(index != j) {

					for(p = 0; p < 3; p++)
						for(q = 0; q < 3; q++)
							atom_array[index]->ef_induced[p] -= system->A_matrix[ii+p][jj+q]*atom_array[j]->mu[q];

				}

			} /* end j */


			/* dipole is the sum of the static and induced parts */
			for(p = 0; p < 3; p++) {

				atom_array[index]->new_mu[p] = atom_array[index]->polarizability*(atom_array[index]->ef_static[p] + atom_array[index]->ef_static_self[p] + atom_array[index]->ef_induced[p]);

				/* Gauss-Seidel */
				if(system->polar_gs || system->polar_gs_ranked)
					atom_array[index]->mu[p] = atom_array[index]->new_mu[p];

			}

		} /* end i */

		/* get the dipole RRMS */
		for(i = 0; i < N; i++) {

			/* mean square difference */
			atom_array[i]->dipole_rrms = 0;
			for(p = 0; p < 3; p++)
				atom_array[i]->dipole_rrms += pow((atom_array[i]->new_mu[p] - atom_array[i]->old_mu[p]), 2.0);

			/* normalize */
			atom_array[i]->dipole_rrms /=  pow(atom_array[i]->new_mu[0], 2.0) + pow(atom_array[i]->new_mu[1], 2.0) + pow(atom_array[i]->new_mu[2], 2.0);
			atom_array[i]->dipole_rrms = sqrt(atom_array[i]->dipole_rrms);

		}

		/* determine if we are done... */
		if(system->polar_precision == 0.0) {	/* ... by fixed iteration ... */

			if(iteration_counter == system->polar_max_iter)
				keep_iterating = 0;
			else
				keep_iterating = 1;

		} else { /* ... or by dipole precision */

			/* immediately reiterate if any component broke tolerance, otherwise we are done */
			for(i = 0, keep_iterating = 0; i < N; i++) {

				/* get the change of dipole between iterations */
				for(p = 0; p < 3; p++) {
					error = pow((atom_array[i]->new_mu[p] - atom_array[i]->old_mu[p]), 2.0);
					if(error > pow(system->polar_precision*DEBYE2SKA, 2.0)) keep_iterating = 1;
				}

			}
		}

		/* contract once more for the next induced field */
		if(system->polar_palmo) {

			/* calculate change in induced field due to this iteration */
			for(i = 0; i < N; i++) {
				index = ranked_array[i];
				ii = index*3;

	                        for(j = 0; j < N; j++) {
					jj = j*3;

					if(index != j) {

						for(p = 0; p < 3; p++)
							for(q = 0; q < 3; q++)
								atom_array[index]->ef_induced_change[p] -= system->A_matrix[ii+p][jj+q]*atom_array[j]->mu[q];
					}
				}
				for(p = 0; p < 3; p++)
					atom_array[index]->ef_induced_change[p] -= atom_array[index]->ef_induced[p];
			} /* end i */
		} /* palmo */

		/* rank the dipoles by bubble sort */
		if(system->polar_gs_ranked) {
			for(i = 0; i < N; i++) {
				for(j = 0, sorted = 1; j < (N-1); j++) {
					if(atom_array[ranked_array[j]]->rank_metric < atom_array[ranked_array[j+1]]->rank_metric) {
						sorted = 0;
						tmp = ranked_array[j];
						ranked_array[j] = ranked_array[j+1];
						ranked_array[j+1] = tmp;
					}
				}
				if(sorted) break;
			}
		}

		/* save the dipoles for the next pass */
		for(i = 0; i < N; i++) {
			for(p = 0; p < 3; p++) {

				/* allow for different successive over-relaxation schemes */
				if(system->polar_sor)
					atom_array[i]->mu[p] = system->polar_gamma*atom_array[i]->new_mu[p] + (1.0 - system->polar_gamma)*atom_array[i]->old_mu[p];
				else if(system->polar_esor)
					atom_array[i]->mu[p] = (1.0 - exp(-system->polar_gamma*iteration_counter))*atom_array[i]->new_mu[p] + exp(-system->polar_gamma*iteration_counter)*atom_array[i]->old_mu[p];
				else
					atom_array[i]->mu[p] = atom_array[i]->new_mu[p];

			}
		}

	} /* end keep iterating */

	free(atom_array);
	free(ranked_array);

	/* return the iteration count */
	return(iteration_counter);

}


