/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

/* calculate the molecular polarizability tensor from the B matrix */

void thole_polarizability_tensor(system_t *system) {
    int i, j, ii, jj, N;
    int p, q;
    double isotropic;

    N = system->checkpoint->thole_N_atom;

    /* clear the polarizability tensor */
    for (p = 0; p < 3; p++)
        for (q = 0; q < 3; q++)
            system->C_matrix[p][q] = 0;

    /* sum the block terms for the 3x3 molecular tensor */
    for (p = 0; p < 3; p++) {
        for (q = 0; q < 3; q++) {
            for (i = 0; i < N; i++) {
                ii = i * 3;
                for (j = 0; j < N; j++) {
                    jj = j * 3;
                    system->C_matrix[p][q] += system->B_matrix[ii + p][jj + q];
                }
            }
        }
    }

    /* get the isotropic term */
    for (p = 0, isotropic = 0; p < 3; p++)
        isotropic += system->C_matrix[p][p];
    isotropic /= 3.0;

    printf(
        "POLARIZATION: polarizability tensor (A^3):\n");
    fflush(stdout);
    printf(
        "##########################\n");
    for (p = 0; p < 3; p++) {
        for (q = 0; q < 3; q++)
            printf(
                "%.4f ", system->C_matrix[p][q]);
        printf(
            "\n");
    }
    printf(
        "##########################\n");
    printf(
        "isotropic = %.4f\n", isotropic);
    printf(
        "XX/ZZ = %.4f\n", system->C_matrix[0][0] / system->C_matrix[2][2]);
}
