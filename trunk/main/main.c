/* 

Space Research Group
Department of Chemistry
University of South Florida

*/

int rank, size;

#include <mc.h>

//kill MPI before quitting, when neccessary
void die( int code ){
#ifdef MPI
	MPI_Finalize();
#endif
	exit(code);
}

void usage(char *progname) {

	if(!rank) fprintf(stderr, "usage: %s <config>\n", progname);
	exit(1);
}

int main(int argc, char **argv) {

	output("MPMC (Massively Parallel Monte Carlo) 2012 GNU Public License\n");
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

	/* start up the MPI chain */
#ifdef MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	output(linebuf);
#endif /* MPI */

	output("For version info, use \"svn info\"\n" );
	sprintf(linebuf, "MAIN: processes started on %d cores\n", size);
	output(linebuf);

	/* get the config file arg */
	strcpy(input_file, argv[1]);
	sprintf(linebuf, "MAIN: running parameters found in %s\n", input_file);
	output(linebuf);

	/* output warning about PDB --> PQR change */
	output("MAIN: *** PLEASE NOTE THAT THE PDB FILE FORMAT IS NO LONGER SUPPORTED IN MPMC ***\n" );

	/* read the input files and setup the simulation */
	system = setup_system(input_file);
	if(!system) {
		error("MAIN: could not initialize the simulation\n");
		die(1);
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
	memnullcheck(system->nodestats,sizeof(nodestats_t), __LINE__-1, __FILE__);
	system->avg_nodestats = calloc(1, sizeof(avg_nodestats_t));
	memnullcheck(system->avg_nodestats,sizeof(avg_nodestats_t), __LINE__-1, __FILE__);
	system->observables = calloc(1, sizeof(observables_t));
	memnullcheck(system->observables,sizeof(observables_t), __LINE__-1, __FILE__);
	system->avg_observables = calloc(1, sizeof(avg_observables_t));
	memnullcheck(system->avg_observables,sizeof(avg_observables_t), __LINE__-1, __FILE__);
	system->checkpoint = calloc(1, sizeof(checkpoint_t));
	memnullcheck(system->checkpoint,sizeof(checkpoint_t), __LINE__-1, __FILE__);
	system->checkpoint->observables = calloc(1, sizeof(observables_t));
	memnullcheck(system->checkpoint->observables,sizeof(observables_t), __LINE__-1, __FILE__);
	system->grids = calloc(1,sizeof(grid_t));
	memnullcheck(system->grids,sizeof(grid_t), __LINE__-1, __FILE__);
	system->grids->histogram = calloc(1,sizeof(histogram_t));
	memnullcheck(system->grids->histogram,sizeof(histogram_t), __LINE__-1, __FILE__);
	system->grids->avg_histogram = calloc(1,sizeof(histogram_t));
	memnullcheck(system->grids->avg_histogram,sizeof(histogram_t), __LINE__-1, __FILE__);

	/* if polarization active, allocate the necessary matrices */
	if(system->polarization && !system->cuda)
		allocate_thole_matrices(system);

	/* if histogram calculation flag is set, allocate grid */
	if(system->calc_hist){
		setup_histogram(system);
		allocate_histogram_grid(system);
	}

	/* seed the rng if neccessary */
	if ( system->ensemble != ENSEMBLE_TE && system->ensemble != ENSEMBLE_REPLAY )
		seed_rng(system, rank);

#ifdef MPI
	MPI_Barrier(MPI_COMM_WORLD);
	sprintf(linebuf, "MAIN: all %d cores are in sync\n", size);
	output(linebuf);
#endif /* MPI */

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
	} else if(system->ensemble == ENSEMBLE_REPLAY) {	/* surface fitting run */
		output("MAIN: **********************************\n");
		output("MAIN: *** starting trajectory replay ***\n");
		output("MAIN: **********************************\n");
	} else if(system->ensemble == ENSEMBLE_SURF_FIT) {	/* surface fitting run */
		output("MAIN: *************************************************************\n");
		output("MAIN: *** starting potential energy surface fitting calculation ***\n");
		output("MAIN: *************************************************************\n");
	} else if(system->ensemble == ENSEMBLE_TE) {	/* surface fitting run */
		output("MAIN: *************************************************\n");
		output("MAIN: *** starting single-point energy calculation  ***\n");
		output("MAIN: *************************************************\n");
	}

	if(system->ensemble == ENSEMBLE_SURF) { /* surface */
		if(surface(system) < 0) {
			error("MAIN: surface module failed on error, exiting\n");
			die(1);
		}
	}
	else if(system->ensemble == ENSEMBLE_SURF_FIT) { /* surface fitting */
		if(surface_fit(system) < 0) {
			error("MAIN: surface fitting module failed on error, exiting\n");
			die(1);
		}
	}
	else if(system->ensemble == ENSEMBLE_REPLAY) { /* replay trajectory and recalc energies, etc. */
		if(replay_trajectory(system) < 0) {
			error("MAIN: trajectory replay failed, exiting\n");
			die(1);
		}
	}
	else if(system->ensemble == ENSEMBLE_TE) {
		if(calculate_te(system) < 0) {
			error("MAIN: single-point energy calculation failed, exiting\n");
			die(1);
		}
	}
	else { //else run monte carlo
		if(mc(system) < 0) {
			error("MAIN: MC failed on error, exiting\n");
			die(1);
		}
	}

	/* cleanup */
	output("MAIN: freeing all data structures....");
	cleanup(system);
	output("...done\n");

	output("MAIN: simulation exiting successfully\n\n");
	die(0);

}
