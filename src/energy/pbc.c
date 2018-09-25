/* 

Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

// calculate the minimum cutoff radius from the basis lattice
// aka shortest vector problem
double pbc_cutoff(pbc_t *pbc) {
    int i, j, k, p;
    double curr_mag;
    double short_mag = MAXVALUE;
    double curr_vec[3];

    if (pbc->volume <= 0) return MAXVALUE;

    for (i = -MAX_VECT_COEF; i <= MAX_VECT_COEF; i++) {
        for (j = -MAX_VECT_COEF; j <= MAX_VECT_COEF; j++) {
            for (k = -MAX_VECT_COEF; k <= MAX_VECT_COEF; k++) {
                if (i == 0 && j == 0 && k == 0) continue;
                for (p = 0; p < 3; p++)
                    curr_vec[p] = i * pbc->basis[0][p] + j * pbc->basis[1][p] + k * pbc->basis[2][p];
                curr_mag = sqrt(dddotprod(curr_vec, curr_vec));
                if (curr_mag < short_mag) short_mag = curr_mag;
            }
        }
    }
    return (0.5 * short_mag);
}

/* take the determinant of the basis matrix */
double pbc_volume(pbc_t *pbc) {
    double volume;

    volume = pbc->basis[0][0] * (pbc->basis[1][1] * pbc->basis[2][2] - pbc->basis[1][2] * pbc->basis[2][1]);
    volume += pbc->basis[0][1] * (pbc->basis[1][2] * pbc->basis[2][0] - pbc->basis[1][0] * pbc->basis[2][2]);
    volume += pbc->basis[0][2] * (pbc->basis[1][0] * pbc->basis[2][1] - pbc->basis[1][1] * pbc->basis[2][0]);

    return (volume);
}

/* get the reciprocal space basis */
void pbc_reciprocal(pbc_t *pbc) {
    double inverse_volume;

    inverse_volume = 1.0 / pbc_volume(pbc);

    pbc->reciprocal_basis[0][0] = inverse_volume * (pbc->basis[1][1] * pbc->basis[2][2] - pbc->basis[1][2] * pbc->basis[2][1]);
    pbc->reciprocal_basis[0][1] = inverse_volume * (pbc->basis[0][2] * pbc->basis[2][1] - pbc->basis[0][1] * pbc->basis[2][2]);
    pbc->reciprocal_basis[0][2] = inverse_volume * (pbc->basis[0][1] * pbc->basis[1][2] - pbc->basis[0][2] * pbc->basis[1][1]);

    pbc->reciprocal_basis[1][0] = inverse_volume * (pbc->basis[1][2] * pbc->basis[2][0] - pbc->basis[1][0] * pbc->basis[2][2]);
    pbc->reciprocal_basis[1][1] = inverse_volume * (pbc->basis[0][0] * pbc->basis[2][2] - pbc->basis[0][2] * pbc->basis[2][0]);
    pbc->reciprocal_basis[1][2] = inverse_volume * (pbc->basis[0][2] * pbc->basis[1][0] - pbc->basis[0][0] * pbc->basis[1][2]);

    pbc->reciprocal_basis[2][0] = inverse_volume * (pbc->basis[1][0] * pbc->basis[2][1] - pbc->basis[1][1] * pbc->basis[2][0]);
    pbc->reciprocal_basis[2][1] = inverse_volume * (pbc->basis[0][1] * pbc->basis[2][0] - pbc->basis[0][0] * pbc->basis[2][1]);
    pbc->reciprocal_basis[2][2] = inverse_volume * (pbc->basis[0][0] * pbc->basis[1][1] - pbc->basis[0][1] * pbc->basis[1][0]);
}

void pbc(system_t *system) {
    pbc_t *pbc = system->pbc;

    /* get the unit cell volume and cutoff */
    pbc->volume = pbc_volume(pbc);

    if (pbc->cutoff == 0.) pbc->cutoff = pbc_cutoff(pbc);

    // calculate ewald_alpha and polar_ewald_alpha unless manually set
    if (system->ewald_alpha_set != 1)
        system->ewald_alpha = 3.5 / system->pbc->cutoff;
    if (system->polar_ewald_alpha_set != 1)
        system->polar_ewald_alpha = 3.5 / system->pbc->cutoff;

    /* get the reciprocal space lattice */
    pbc_reciprocal(pbc);

    return;
}
