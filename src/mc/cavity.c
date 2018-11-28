/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

/* check whether a point (x,y,z) lies within an empty cavity */
/* if so, return 1 */
int is_point_empty(system_t *system, double x, double y, double z) {
    int i, j, k;
    int incavity;
    double r;

    incavity = 0;
    for (i = 0; i < system->cavity_grid_size; i++) {
        for (j = 0; j < system->cavity_grid_size; j++) {
            for (k = 0; k < system->cavity_grid_size; k++) {
                if (!system->cavity_grid[i][j][k].occupancy) {
                    r = pow((x - system->cavity_grid[i][j][k].pos[0]), 2);
                    r += pow((y - system->cavity_grid[i][j][k].pos[1]), 2);
                    r += pow((z - system->cavity_grid[i][j][k].pos[2]), 2);
                    r = sqrt(r);

                    if (r < system->cavity_radius)
                        incavity = 1;
                }

            } /* end i */
        }     /* end j */
    }         /* end k */

    return (incavity);
}

/* total volume of accessible cavities via Monte Carlo integration */
void cavity_volume(system_t *system) {
    int throw, hits, num_darts;
    int p, q;
    double pos_vec[3], grid_vec[3];
    double fraction_hits;

    /* good rule of thumb is 1 per 10 A^3 */
    num_darts = system->pbc->volume * DARTSCALE;

    /* throw random darts and count the number of hits */
    for (throw = 0, hits = 0; throw < num_darts; throw ++) {
        /* generate a random grid position */
        for (p = 0; p < 3; p++) grid_vec[p] = -0.5 + get_rand(system);

        /* zero the coordinate vector */
        for (p = 0; p < 3; p++) pos_vec[p] = 0;

        /* linear transform vector into real coordinates */
        for (p = 0; p < 3; p++)
            for (q = 0; q < 3; q++)
                pos_vec[p] += system->pbc->basis[q][p] * grid_vec[q];

        /* check if the random point lies within an empty cavity */
        if (is_point_empty(system, pos_vec[0], pos_vec[1], pos_vec[2])) ++hits;
    }

    /* determine the percentage of free cavity space */
    fraction_hits = ((double)hits) / ((double)num_darts);
    system->cavity_volume = fraction_hits * system->pbc->volume; /* normalize w.r.t. the cell volume */
}

/* probability of finding an empty cavity on the grid */
void cavity_probability(system_t *system) {
    int i, j, k;
    double probability;
    int total_points;

    /* total number of potential cavities */
    total_points = system->cavity_grid_size * system->cavity_grid_size * system->cavity_grid_size;

    /* find the number of open cavities */
    system->cavities_open = 0;
    for (i = 0; i < system->cavity_grid_size; i++)
        for (j = 0; j < system->cavity_grid_size; j++)
            for (k = 0; k < system->cavity_grid_size; k++)
                if (!system->cavity_grid[i][j][k].occupancy) ++system->cavities_open;

    /* the overall probability ratio */
    probability = ((double)system->cavities_open) / ((double)total_points);

    /* update the observable */
    system->nodestats->cavity_bias_probability = probability;
}

/* allocate the grid */
void setup_cavity_grid(system_t *system) {
    int i, j;

    system->cavity_grid = calloc(system->cavity_grid_size, sizeof(cavity_t **));
    memnullcheck(system->cavity_grid, system->cavity_grid_size * sizeof(cavity_t **), __LINE__ - 1, __FILE__);
    for (i = 0; i < system->cavity_grid_size; i++) {
        system->cavity_grid[i] = calloc(system->cavity_grid_size, sizeof(cavity_t *));
        memnullcheck(system->cavity_grid[i], system->cavity_grid_size * sizeof(cavity_t *), __LINE__ - 1, __FILE__);
        for (j = 0; j < system->cavity_grid_size; j++) {
            system->cavity_grid[i][j] = calloc(system->cavity_grid_size, sizeof(cavity_t));
            memnullcheck(system->cavity_grid[i][j], system->cavity_grid_size * sizeof(cavity_t), __LINE__ - 1, __FILE__);
        }
    }
}

/* create a 3D histogram of atoms lying within a sphere centered at each grid point */
void cavity_update_grid(system_t *system) {
    int i, j, k, G;
    int p, q;
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    double grid_component[3];
    double grid_vector[3];
    double r;

    G = system->cavity_grid_size;

    /* clear the grid */
    for (i = 0; i < G; i++) {
        for (j = 0; j < G; j++) {
            for (k = 0; k < G; k++) {
                system->cavity_grid[i][j][k].occupancy = 0;
                for (p = 0; p < 3; p++)
                    system->cavity_grid[i][j][k].pos[p] = 0;
            }
        }
    }

    /* loop over the grid, bin a sphere when needed */
    for (i = 0; i < G; i++) {
        for (j = 0; j < G; j++) {
            for (k = 0; k < G; k++) {
                /* divide up each grid component */
                grid_component[0] = ((double)(i + 1)) / ((double)(G + 1));
                grid_component[1] = ((double)(j + 1)) / ((double)(G + 1));
                grid_component[2] = ((double)(k + 1)) / ((double)(G + 1));

                /* project the grid point onto our actual basis */
                for (p = 0; p < 3; p++)
                    for (q = 0, grid_vector[p] = 0; q < 3; q++)
                        grid_vector[p] += system->pbc->basis[q][p] * grid_component[q];

                /* put into real coordinates */
                for (p = 0; p < 3; p++)
                    for (q = 0; q < 3; q++)
                        grid_vector[p] -= 0.5 * system->pbc->basis[q][p];

                /* if an atomic coordinate lies within a sphere centered on the grid point, then bin it */
                for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
                    for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
                        /* get the displacement from the grid point */
                        for (p = 0, r = 0; p < 3; p++)
                            r += (grid_vector[p] - atom_ptr->wrapped_pos[p]) * (grid_vector[p] - atom_ptr->wrapped_pos[p]);
                        r = sqrt(r);

                        /* inside the sphere? */
                        if (r < system->cavity_radius) ++system->cavity_grid[i][j][k].occupancy;

                    } /* for atom */
                }     /* for molecule */

                /* store the location of this grid point */
                for (p = 0; p < 3; p++)
                    system->cavity_grid[i][j][k].pos[p] = grid_vector[p];

            } /* for k */
        }     /* for j */
    }         /* for i */

    /* update the cavity insertion probability estimate */
    cavity_probability(system);
    /* update the accessible insertion volume */
    cavity_volume(system);
}

#ifdef DEBUG
void test_cavity_grid(system_t *system) {
    int i, j, k, G;

    G = system->cavity_grid_size;

    printf(
        "\n");
    for (i = 0; i < G; i++) {
        printf(
            "DEBUG_CAVITY: ");
        for (j = 0; j < G; j++) {
            for (k = 0; k < G; k++)
                printf(
                    "%d ", system->cavity_grid[i][j][k].occupancy);
        }
        printf(
            "\n");
    }
    printf(
        "\n");
    fflush(stdout);
}
#endif /* DEBUG */
