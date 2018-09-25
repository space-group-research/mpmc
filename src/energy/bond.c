/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

/* energy of molecule in a 1D anharmonic well */
double anharmonic(system_t *system) {
    double k, g, x;
    double energy;
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;

    k = system->rd_anharmonic_k;
    g = system->rd_anharmonic_g;

    energy = 0;
    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            x = atom_ptr->pos[0];
            energy += anharmonic_energy(k, g, x);

            if (system->feynman_hibbs) {
                if (system->feynman_kleinert) {
                    energy = anharmonic_fk(system->temperature, atom_ptr->mass, k, g, x);

                } else {
                    if (system->feynman_hibbs_order == 2)
                        energy += anharmonic_fh_second_order(system->temperature, atom_ptr->mass, k, g, x);
                    else if (system->feynman_hibbs_order == 4)
                        energy += anharmonic_fh_fourth_order(system->temperature, atom_ptr->mass, k, g, x);
                }
            }
        }
    }

    return (energy);
}

/* Feynman-Kleinert iterative method of effective potential */
double anharmonic_fk(double temperature, double mass, double k, double g, double x) {
    int keep_iterating;
    double a_sq;              /* width a^2 (A^2) */
    double omega, omega_sq;   /* spacing Omega^2 (K/A^2) */
    double prev_a_sq;         /* last a_sq */
    double tolerance;         /* iterative tolerance */
    double V_a;               /* V_a^2 (K) */
    double potential;         /* W_1 (K) */
    double conversion_factor; /* hbar^2*(m2A)^2/(k*m) */

    /* convert the mass to kg */
    mass *= AMU2KG;

    conversion_factor = pow(METER2ANGSTROM, 2) * pow(HBAR, 2) / (KB * mass);

    /* initial guess a^2 = beta/12 */
    a_sq = pow(METER2ANGSTROM, 2) * pow(HBAR, 2) / (12.0 * KB * temperature * mass);

    /* solve self-consistently */
    keep_iterating = 1;
    while (keep_iterating) {
        /* save the last a_sq for tolerance */
        prev_a_sq = a_sq;

        omega_sq = conversion_factor * (k + 3.0 * g * a_sq + 3.0 * g * pow(x, 2));
        omega = sqrt(omega_sq);
        a_sq = conversion_factor * (temperature / omega_sq) * ((omega / (2.0 * temperature)) * (1.0 / tanh(omega / (2.0 * temperature))) - 1.0);

        tolerance = fabs(prev_a_sq - a_sq);
        if (tolerance < FEYNMAN_KLEINERT_TOLERANCE)
            keep_iterating = 0;
    }

    V_a = 0.5 * a_sq * k + 0.75 * g * pow(a_sq, 2) + 0.5 * (k + 3.0 * g * a_sq) * pow(x, 2) + 0.25 * g * pow(x, 4);

    potential = temperature * log(sinh(omega / (2.0 * temperature)) / (omega / (2.0 * temperature))) - 0.5 * omega_sq * a_sq / conversion_factor + V_a;

    return (potential);
}

/* up to FH h^2 effective potential term */
double anharmonic_fh_second_order(double temperature, double mass, double k, double g, double x) {
    double first_derivative;
    double second_derivative;
    double potential;

    mass *= AMU2KG;

    first_derivative = k * x + g * pow(x, 3);
    second_derivative = k + 3.0 * g * pow(x, 2);

    potential = pow(METER2ANGSTROM, 2) * pow(HBAR, 2) / (24.0 * KB * temperature * mass) * (second_derivative + 2.0 * first_derivative / x);

    return (potential);
}

/* up to FH h^4 effective potential term */
double anharmonic_fh_fourth_order(double temperature, double mass, double k, double g, double x) {
    double first_derivative, second_derivative;
    double other_derivatives;
    double potential;

    mass *= AMU2KG;

    first_derivative = k * x + g * pow(x, 3);
    second_derivative = k + 3.0 * g * pow(x, 2);
    other_derivatives = 15.0 * k / pow(x, 2) + 45.0 * g;

    potential = pow(METER2ANGSTROM, 2) * pow(HBAR, 2) / (24.0 * KB * temperature * mass) * (second_derivative + 2.0 * first_derivative / x);
    potential += pow(METER2ANGSTROM, 4) * pow(HBAR, 4) / (1152.0 * pow(KB * temperature * mass, 2)) * other_derivatives;

    return (potential);
}

/* anharmonic potential */
double anharmonic_energy(double k, double g, double x) {
    double potential;

    potential = 0.5 * k * pow(x, 2) + 0.25 * g * pow(x, 4);

    return (potential);
}

/* H2 specific function */
double h2_bond_energy(double r) {
    return (morse_energy(H2_SPRING_CONSTANT, H2_De, H2_R0, r));
}

/* morse potential */
/* input args are spring constant k in (eV*A^-2), */
/* dissociation energy (in eV), equil r (in A) and the r value (in A) */
double morse_energy(double k, double de, double r0, double r) {
    double a, potential;

    a = sqrt(k / (2.0 * de));

    potential = (EV2K * de) * pow((1.0 - exp(-a * (r - r0))), 2);

    return (potential);
}
