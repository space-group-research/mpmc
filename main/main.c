/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/
#define VERSION "0.11\n\tMC insertions selected from a separate linked list."


int rank, size;

#include <mc.h>


void usage(char *progname) {

	if(!rank) fprintf(stderr, "usage: %s <config>\n", progname);

#ifdef MPI
	MPI_Finalize();
#else
	exit(1);
#endif /* MPI */

}

void seed_rng(long int seed) {
	rule30_rng(seed);
}

double get_rand(void) {
	return(rule30_rng(0));
}

int main(int argc, char **argv) {
        printf( "MPMC Version id: %s\n", VERSION );
        int i, j, N;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	char linebuf[MAXLINE];
	char input_file[MAXLINE];
	system_t *system;
	


	/* set the default rank */
	rank = 0; size = 1;

	/* check args */
	if(argc < 2) usage(argv[0]);

	if(!argv[1]) {
		error("MAIN: invalid config file specified");
		exit(1);
	}

	/* start up the MPI chain */
#ifdef MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	sprintf(linebuf, "MAIN: processes started on %d cores\n", size);
	output(linebuf);
#endif /* MPI */

	/* get the config file arg */
	strcpy(input_file, argv[1]);
	sprintf(linebuf, "MAIN: running parameters found in %s\n", input_file);
	output(linebuf);

	// moved this allocation, because we'd like to set system->obsevables->volume during setup_system

	/* read the input files and setup the simulation */
	system = setup_system(input_file);
	if(!system) {
		error("MAIN: could not initialize the simulation\n");
#ifdef MPI
		MPI_Finalize();
#else
		exit(1);
#endif /* MPI */
	} else {
		output("MAIN: the simulation has been initialized\n");
	}
	
	/* install the signal handler to catch SIGTERM cleanly */
	terminate_handler(-1, system);
	signal(SIGTERM, ((void *)(terminate_handler)));
	signal(SIGUSR1, ((void *)(terminate_handler)));
	signal(SIGUSR2, ((void *)(terminate_handler)));
	output("MAIN: signal handler installed\n");

	/* allocate space for the statistics */
	system->nodestats = calloc(1, sizeof(nodestats_t));
	memnullcheck(system->nodestats,sizeof(nodestats_t),21);
	system->avg_nodestats = calloc(1, sizeof(avg_nodestats_t));
	memnullcheck(system->avg_nodestats,sizeof(avg_nodestats_t), 22);
	system->observables = calloc(1, sizeof(observables_t));
	memnullcheck(system->observables,sizeof(observables_t),23);
	system->avg_observables = calloc(1, sizeof(avg_observables_t));
	memnullcheck(system->avg_observables,sizeof(avg_observables_t),24);
	system->checkpoint = calloc(1, sizeof(checkpoint_t));
	memnullcheck(system->checkpoint,sizeof(checkpoint_t),25);
	system->checkpoint->observables = calloc(1, sizeof(observables_t));
	memnullcheck(system->checkpoint->observables,sizeof(observables_t),26);
	system->grids = calloc(1,sizeof(grid_t));
	memnullcheck(system->grids,sizeof(grid_t),27);
	system->grids->histogram = calloc(1,sizeof(histogram_t));
	memnullcheck(system->grids->histogram,sizeof(histogram_t),28);
	system->grids->avg_histogram = calloc(1,sizeof(histogram_t));
	memnullcheck(system->grids->avg_histogram,sizeof(histogram_t),29);

	/* if polarization active, allocate the necessary matrices */
	if(system->polarization && !system->cuda) {

		/* count the number of atoms initially in the system */
		for(molecule_ptr = system->molecules, N = 0; molecule_ptr; molecule_ptr = molecule_ptr->next)
			for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
				++N;

		system->A_matrix = calloc(3*N, sizeof(double *));
		memnullcheck(system->A_matrix,3*N*sizeof(double *),30);
		for(i = 0; i < 3*N; i++) {
			system->A_matrix[i] = calloc(3*N, sizeof(double));
			memnullcheck(system->A_matrix[i],3*N*sizeof(double), 31);
		}

		if(!system->polar_iterative) {
			system->B_matrix = calloc(3*N, sizeof(double *));
			memnullcheck(system->B_matrix,3*N*sizeof(double *),32);
			for(i = 0; i < 3*N; i++) {
				system->B_matrix[i] = calloc(3*N, sizeof(double));
				memnullcheck(system->B_matrix[i],3*N*sizeof(double),33);
			}
		}


	}

	/* if histogram calculation flag is set, allocate grid */
	if(system->calc_hist){
		setup_histogram(system);
		allocate_histogram_grid(system);
	}

	/* seed the rng */
	seed_rng(system->seed + rank);
	output("MAIN: the random seeds on the compute cores are all set\n");

#ifdef MPI
	MPI_Barrier(MPI_COMM_WORLD);
	sprintf(linebuf, "MAIN: all %d cores are in sync\n", size);
	output(linebuf);
#endif /* MPI */



	///////////////////////////////////////////////////////////////////
	// Testing Purposes only, for unique tagging of molecules erase
	// when no longer needed.
	{
		molecule_t *mp;
		int nMolecules = 0;
		for( mp = system->molecules; mp; mp = mp->next )
			nMolecules++;
		system->tag = nMolecules + 1;
	}
	///////////////////////////////////////////////////////////////////




	/* start the MC simulation */
	if(system->ensemble == ENSEMBLE_UVT) {
		output("MAIN: *******************************************************\n");
		output("MAIN: *** starting Grand Canonical Monte Carlo simulation ***\n");
		output("MAIN: *******************************************************\n");
	} else if(system->ensemble == ENSEMBLE_NVT) {
		output("MAIN: *************************************************\n");
		output("MAIN: *** starting Canonical Monte Carlo simulation ***\n");
		output("MAIN: *************************************************\n");
	} else if(system->ensemble == ENSEMBLE_NVE) {
		output("MAIN: ******************************************************\n");
		output("MAIN: *** starting Microcanonical Monte Carlo simulation ***\n");
		output("MAIN: ******************************************************\n");
	} else if(system->ensemble == ENSEMBLE_SURF) {	/* surface run */
		output("MAIN: *****************************************************\n");
		output("MAIN: *** starting potential energy surface calculation ***\n");
		output("MAIN: *****************************************************\n");
	} else if(system->ensemble == ENSEMBLE_SURF_FIT) {	/* surface fitting run */
		output("MAIN: *************************************************************\n");
		output("MAIN: *** starting potential energy surface fitting calculation ***\n");
		output("MAIN: *************************************************************\n");
	} else if(system->ensemble == ENSEMBLE_TE) {	/* surface fitting run */
		output("MAIN: *************************************************\n");
		output("MAIN: *** starting single-point energy calculation  ***\n");
		output("MAIN: *************************************************\n");
	}

	if(!((system->ensemble == ENSEMBLE_SURF) || (system->ensemble == ENSEMBLE_SURF_FIT))) {

		if(mc(system) < 0) {
			error("MAIN: MC failed on error, exiting\n");
#ifdef MPI
			MPI_Finalize();
#else
			exit(1);
#endif /* MPI */
		} else {
			if(system->ensemble == ENSEMBLE_UVT) {
				output("MAIN: ********************************************************\n");
				output("MAIN: *** finishing Grand Canonical Monte Carlo simulation ***\n");
				output("MAIN: ********************************************************\n\n");
			} else if(system->ensemble == ENSEMBLE_NVT) {
				output("MAIN: **************************************************\n");
				output("MAIN: *** finishing Canonical Monte Carlo simulation ***\n");
				output("MAIN: **************************************************\n\n");
			} else if(system->ensemble == ENSEMBLE_NVE) {
				output("MAIN: *******************************************************\n");
				output("MAIN: *** finishing Microcanonical Monte Carlo simulation ***\n");
				output("MAIN: *******************************************************\n\n");
			}
		}

	} else if(system->ensemble == ENSEMBLE_SURF) { /* surface */

		if(surface(system) < 0) {
			error("MAIN: surface module failed on error, exiting\n");
#ifdef MPI
			MPI_Finalize();
#else
			exit(1);
#endif /* MPI */
		} else {
			output("MAIN: ******************************************************\n");
			output("MAIN: *** finishing potential energy surface calculation ***\n");
			output("MAIN: ******************************************************\n");
		}

	} else if(system->ensemble == ENSEMBLE_SURF_FIT) { /* surface fitting */

		if(surface_fit(system) < 0) {
			error("MAIN: surface fitting module failed on error, exiting\n");
#ifdef MPI
			MPI_Finalize();
#else
			exit(1);
#endif /* MPI */
		} else {
			output("MAIN: **************************************************************\n");
			output("MAIN: *** finishing potential energy surface fitting calculation ***\n");
			output("MAIN: **************************************************************\n");
		}

	}



	/* cleanup */
	output("MAIN: freeing all data structures....");
	cleanup(system);
	output("...done\n");

	output("MAIN: simulation exiting successfully\n\n");
#ifdef MPI
	MPI_Finalize();
#else
	exit(0);
#endif /* MPI */

}

