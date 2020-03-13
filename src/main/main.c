/* 

Space Research Group
Department of Chemistry
University of South Florida

*/
int rank, size;

#include <mc.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "surface_fit_arbitrary.h"
#include "surface_multi_fit.h"

//kill MPI before quitting, when neccessary
void die(int code) {
#ifdef MPI
    MPI_Finalize();
#endif
    exit(code);
}

void usage(char *progname) {
    if (!rank) {
        fprintf(stderr,
                "usage: %s <config>\n", progname);
        fprintf(stderr,
                "See: https://github.com/mpmccode/mpmc\n");
    }
    exit(1);
}

int main(int argc, char **argv) {
    char linebuf[MAXLINE];
    char input_file[MAXLINE];
    system_t *system;
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);

    /* set the default rank */
    rank = 0;
    size = 1;

    /* check args */
    if (argc < 2) usage(argv[0]);

/* start up the MPI chain */
#ifdef MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif /* MPI */

    if (!rank) {
#ifndef __WIN32__
        sprintf(linebuf,
                "MPMC (Massively Parallel Monte Carlo) r%d - 2012-2019 GNU Public License\n", VERSION);
        output(linebuf);
#endif
        sprintf(linebuf,
                "MAIN: processes started on %d cores @ %d-%d-%d %d:%d:%d\n", size, tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
        output(linebuf);
    }

#ifndef MPI
    FILE *procfile;
    FILE *host;
    char nodename[MAXLINE];
    char cpu[MAXLINE];
    struct stat info;
    char *get_nodename_err;
    char *get_cpu_err;

    // These system calls were causing a fork() that does not play nice with some
    // MPI implementations (causing some processes to never end... )

    /* Collect and output hostname & processor info */
    if (access(
            "/proc/cpuinfo", R_OK) != -1) {  // Linux
        host = popen(
            "hostname",
            "r");
        get_nodename_err = fgets(nodename, MAXLINE, host);
        sprintf(linebuf,
                "MAIN: Job running on node -> %.400s", nodename);
        output(linebuf);
        pclose(host);

        procfile = fopen(
            "/proc/cpuinfo",
            "r");
        while (!feof(procfile)) {
            get_cpu_err = fgets(cpu, MAXLINE, procfile);
            if (strncasecmp(cpu,
                            "model name", 10) == 0) {
                sprintf(linebuf,
                        "MAIN: CPU -> %.400s", cpu);
                output(linebuf);
                break;
            }
        }
        fclose(procfile);
    } else if (stat(
                   "/Applications", &info) == 0) {  // Mac OS
        output(
            "MAIN: Mac OS detected\n");
        host = popen(
            "hostname",
            "r");
        get_nodename_err = fgets(nodename, MAXLINE, host);
        sprintf(linebuf,
                "MAIN: Job running on node -> %.400s", nodename);
        output(linebuf);
        pclose(host);

        procfile = popen(
            "sysctl -n machdep.cpu.brand_string",
            "r");
        get_cpu_err = fgets(cpu, MAXLINE, procfile);
        sprintf(linebuf,
                "MAIN: CPU -> %.400s", cpu);
        output(linebuf);
        pclose(procfile);
    }
#endif  // !MPI

    /* get the config file arg */
    strcpy(input_file, argv[1]);
    sprintf(linebuf,
            "MAIN: running parameters found in %.400s\n", input_file);
    output(linebuf);

    /* read the input files and setup the simulation */
    system = setup_system(input_file);
    if (!system) {
        error(
            "MAIN: could not initialize the simulation\n");
        die(1);
    } else {
        output(
            "MAIN: the simulation has been initialized\n");
    }

    /* install the signal handler to catch SIGTERM cleanly */
    terminate_handler(-1, system);
    signal(SIGTERM, ((void *)(terminate_handler)));
#ifdef __linux__
    signal(SIGUSR1, ((void *)(terminate_handler)));
    signal(SIGUSR2, ((void *)(terminate_handler)));
#endif
    output(
        "MAIN: signal handler installed\n");

    /* allocate space for the statistics */
    system->nodestats = calloc(1, sizeof(nodestats_t));
    memnullcheck(system->nodestats, sizeof(nodestats_t), __LINE__ - 1, __FILE__);
    system->avg_nodestats = calloc(1, sizeof(avg_nodestats_t));
    memnullcheck(system->avg_nodestats, sizeof(avg_nodestats_t), __LINE__ - 1, __FILE__);
    system->observables = calloc(1, sizeof(observables_t));
    memnullcheck(system->observables, sizeof(observables_t), __LINE__ - 1, __FILE__);
    system->avg_observables = calloc(1, sizeof(avg_observables_t));
    memnullcheck(system->avg_observables, sizeof(avg_observables_t), __LINE__ - 1, __FILE__);
    system->checkpoint = calloc(1, sizeof(checkpoint_t));
    memnullcheck(system->checkpoint, sizeof(checkpoint_t), __LINE__ - 1, __FILE__);
    system->checkpoint->observables = calloc(1, sizeof(observables_t));
    memnullcheck(system->checkpoint->observables, sizeof(observables_t), __LINE__ - 1, __FILE__);
    system->grids = calloc(1, sizeof(grid_t));
    memnullcheck(system->grids, sizeof(grid_t), __LINE__ - 1, __FILE__);
    system->grids->histogram = calloc(1, sizeof(histogram_t));
    memnullcheck(system->grids->histogram, sizeof(histogram_t), __LINE__ - 1, __FILE__);
    system->grids->avg_histogram = calloc(1, sizeof(histogram_t));
    memnullcheck(system->grids->avg_histogram, sizeof(histogram_t), __LINE__ - 1, __FILE__);

    /* if polarization active, allocate the necessary matrices */
    if (system->polarization && !system->cuda && !system->polar_zodid)
        thole_resize_matrices(system);

    /* if histogram calculation flag is set, allocate grid */
    if (system->calc_hist) {
        setup_histogram(system);
        allocate_histogram_grid(system);
    }

#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
    sprintf(linebuf,
            "MAIN: all %d cores are in sync\n", size);
    output(linebuf);
#endif /* MPI */

    /* start the MC simulation */
    if (system->ensemble == ENSEMBLE_UVT) {
        output(
            "MAIN: *******************************************************\n");
        output(
            "MAIN: *** starting Grand Canonical Monte Carlo simulation ***\n");
        output(
            "MAIN: *******************************************************\n");
    } else if (system->ensemble == ENSEMBLE_NVT) {
        output(
            "MAIN: *************************************************\n");
        output(
            "MAIN: *** starting Canonical Monte Carlo simulation ***\n");
        output(
            "MAIN: *************************************************\n");
    } else if (system->ensemble == ENSEMBLE_NVE) {
        output(
            "MAIN: ******************************************************\n");
        output(
            "MAIN: *** starting Microcanonical Monte Carlo simulation ***\n");
        output(
            "MAIN: ******************************************************\n");
    } else if (system->ensemble == ENSEMBLE_SURF) { /* surface run */
        output(
            "MAIN: *****************************************************\n");
        output(
            "MAIN: *** starting potential energy surface calculation ***\n");
        output(
            "MAIN: *****************************************************\n");
    } else if (system->ensemble == ENSEMBLE_REPLAY) { /* surface fitting run */
        output(
            "MAIN: **********************************\n");
        output(
            "MAIN: *** starting trajectory replay ***\n");
        output(
            "MAIN: **********************************\n");
    } else if (system->ensemble == ENSEMBLE_SURF_FIT) { /* surface fitting run */
        output(
            "MAIN: *************************************************************\n");
        output(
            "MAIN: *** starting potential energy surface fitting calculation ***\n");
        output(
            "MAIN: *************************************************************\n");
    } else if (system->ensemble == ENSEMBLE_TE) { /* surface fitting run */
        output(
            "MAIN: *************************************************\n");
        output(
            "MAIN: *** starting single-point energy calculation  ***\n");
        output(
            "MAIN: *************************************************\n");
    }

    if (system->ensemble == ENSEMBLE_SURF) { /* surface */
        if (surface(system) < 0) {
            error(
                "MAIN: surface module failed on error, exiting\n");
            die(1);
        }
    }

    else if (system->ensemble == ENSEMBLE_SURF_FIT) { /* surface fitting */

        if (system->surf_fit_multi_configs) {
            if (surface_multi_fit(system) < 0) {
                error("MAIN: surface fitting module (for multiple configurations) failed on error, exiting\n");
                die(1);
            }
        } else if (system->surf_fit_arbitrary_configs) {
            if (surface_fit_arbitrary(system) < 0) {
                error(
                    "MAIN: surface fitting module (for arbitrary configurations) failed on error, exiting\n");
                die(1);
            }
        } else if (surface_fit(system) < 0) {
            error(
                "MAIN: surface fitting module failed on error, exiting\n");
            die(1);
        }
    }

    else if (system->ensemble == ENSEMBLE_REPLAY) { /* replay trajectory and recalc energies, etc. */
        if (replay_trajectory(system) < 0) {
            error(
                "MAIN: trajectory replay failed, exiting\n");
            die(1);
        }
    }

    else if (system->ensemble == ENSEMBLE_TE) {
        if (calculate_te(system) < 0) {
            error(
                "MAIN: single-point energy calculation failed, exiting\n");
            die(1);
        }
    }

    else {  //else run monte carlo
        if (mc(system) < 0) {
            error(
                "MAIN: MC failed on error, exiting\n");
            die(1);
        }
    }

    /* cleanup */
    output(
        "MAIN: freeing all data structures....");
    cleanup(system);
    output(
        "...done\n");

    output(
        "MAIN: simulation exiting successfully\n\n");
    die(0);

    return 0;
}
