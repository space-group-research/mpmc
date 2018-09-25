/* 

Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

void mpi_copy_histogram_to_sendbuffer(char *snd, int ***grid, system_t *system) {
    int i, j, k;
    int x_dim = system->grids->histogram->x_dim;
    int y_dim = system->grids->histogram->y_dim;
    int z_dim = system->grids->histogram->z_dim;
    int *sndcast = (int *)snd; /* cast the send buffer address to an (int *) */

    /* histogram */
    for (k = 0; k < z_dim; k++) {
        for (j = 0; j < y_dim; j++) {
            for (i = 0; i < x_dim; i++) {
                sndcast[i + j * x_dim + k * x_dim * y_dim] = system->grids->histogram->grid[i][j][k];
            }
        }
    }
}

void mpi_copy_rcv_histogram_to_data(char *rcv, int ***histogram, system_t *system) {
    int i, j, k;
    int x_dim = system->grids->histogram->x_dim;
    int y_dim = system->grids->histogram->y_dim;
    int z_dim = system->grids->histogram->z_dim;
    int *rcvcast = (int *)rcv;

    /* histogram */
    for (k = 0; k < z_dim; k++) {
        for (j = 0; j < y_dim; j++) {
            for (i = 0; i < x_dim; i++) {
                histogram[i][j][k] = rcvcast[i + j * x_dim + k * x_dim * y_dim];
            }
        }
    }
}
