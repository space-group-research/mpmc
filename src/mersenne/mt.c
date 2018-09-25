#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dSFMT.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include <mc.h>
#ifdef MPI
#include <mpi.h>
#endif

#define RANDARRAYSIZE 1000

//global variables (sorry bout that)
int rand_array_used, rand_array_size;
double* rand_array;
dsfmt_t rand_state;

uint32_t getseed(void) {
    uint32_t urand;
    uint32_t rval;

    urand = open(
        "/dev/urandom", O_RDONLY);
    if (urand <= 0) {
        error(
            "MT: failed to open /dev/urandom\n");
        die(-1);
    }
    if (!read(urand, &rval, sizeof(rval))) {
        error(
            "MT: failed to read from /dev/urandom\n");
        die(-1);
    }
    close(urand);

    return rval;
}

void seed_mt_rand(system_t* system, dsfmt_t* state, int rank) {
    static uint32_t* seeds = NULL;
    int i;

    seeds = malloc(4 * sizeof(uint32_t));
    memnullcheck(seeds, 4 * sizeof(uint32_t), __LINE__ - 1, __FILE__);

    if (system->preset_seeds_on) {
        for (i = 0; i < 4; i++)
            seeds[i] = system->preset_seeds[i];
        fprintf(stdout,
                "MT: seeds set manually\n");
    } else {
        for (i = 0; i < 4; i++)
            seeds[i] = getseed();
    }

#ifdef MPI
    {
        char out[MAXLINE];
        int j;
        sprintf(out,
                "MT: core[%d] seeds are set to %u %u %u %u\n", rank, seeds[0], seeds[1], seeds[2], seeds[3]);
        for (j = 0; j < size; j++) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (j == rank) {
                printf(
                    "%s", out);
                fflush(stdout);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#else
    fprintf(stdout,
            "MT: core[%d] seeds are set to %u %u %u %u\n", rank, seeds[0], seeds[1], seeds[2], seeds[3]);
#endif

    //seed the rng
    dsfmt_init_by_array(state, seeds, 4);
    free(seeds);
    return;
}

int refill_mt_rands(dsfmt_t* state, double** rand_array) {
    static int rand_array_size = 0;

    if (rand_array_size == 0) {
        rand_array_size = dsfmt_get_min_array_size();
        rand_array_size *= ceil((double)RANDARRAYSIZE / (double)rand_array_size);  //generate size rands at a time
        *rand_array = malloc(rand_array_size * sizeof(double));
        memnullcheck(rand_array, rand_array_size * sizeof(double), __LINE__ - 1, __FILE__);

        //check null
    }

    dsfmt_fill_array_close_open(state, *rand_array, rand_array_size);  //fill array

    return rand_array_size;
}

//seed the mersenne twister rng
void seed_rng(system_t* system, int rank) {
    seed_mt_rand(system, &rand_state, rank);
    rand_array_size = 0;
    rand_array_used = INT_MAX;
    return;
}

//get the next random number from the mersenne twister
double get_rand() {
    if (rand_array_used >= rand_array_size) {  //need more rands
        rand_array_used = 0;
        rand_array_size = refill_mt_rands(&rand_state, &rand_array);
    }
    rand_array_used++;

    return rand_array[rand_array_used - 1];
}

void kill_rng() {
    free(rand_array);
    return;
}
