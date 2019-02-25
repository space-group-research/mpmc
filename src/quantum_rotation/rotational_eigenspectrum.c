/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

/* returns the hamiltonian matrix for the potential integrated over the Ylm basis */
complex_t **rotational_hamiltonian(system_t *system, molecule_t *molecule, int level_max, int l_max) {
    int dim;
    int i, j;
    int li, mi, lj, mj;
    complex_t **hamiltonian;

    /* determine the size of the symmetric matrix */
    dim = (l_max + 1) * (l_max + 1);

    /* allocate our hamiltonian matrix */
    hamiltonian = calloc(dim, sizeof(complex_t *));
    memnullcheck(hamiltonian, dim * sizeof(complex_t *), __LINE__ - 1, __FILE__ - 1);
    for (i = 0; i < dim; i++) {
        hamiltonian[i] = calloc(dim, sizeof(complex_t));
        memnullcheck(hamiltonian[i], dim * sizeof(complex_t), __LINE__ - 1, __FILE__ - 1);
    }

    /* construct the hamiltonian matrix, integrate to get each element */
    for (i = 0, li = 0; i < dim; i++) {
        for (j = i, lj = li; j < dim; j++) {
            /* determine m values */
            mi = i - li * (li + 1);
            mj = j - lj * (lj + 1);

            /* calculate the dirac braket of <li,mi|V|lj,mj> */
            hamiltonian[i][j].real = rotational_integrate(system, molecule, REAL, li, mi, lj, mj);
            hamiltonian[i][j].imaginary = rotational_integrate(system, molecule, IMAGINARY, li, mi, lj, mj);

            /* store the lower half */
            if (i != j) {
                hamiltonian[j][i].real = hamiltonian[i][j].real;
                hamiltonian[j][i].imaginary = -1.0 * hamiltonian[i][j].imaginary;
            } else {
                /* add the kinetic part to the diagonal elements */
                hamiltonian[i][j].real += system->quantum_rotation_B * ((double)(li * (li + 1)));
            }

/* output the matrix */
#ifdef DEBUG
            printf(
                "DEBUG_QROT: rotmat %d %d %.6f %.6f K (%.6f %.6f cm^-1) <%d %d | V | %d %d>\n",
                (i + 1), (j + 1), hamiltonian[i][j].real, hamiltonian[i][j].imaginary,
                hamiltonian[i][j].real * K2WN, hamiltonian[i][j].imaginary * K2WN, li, mi, lj, mj);
            fflush(stdout);
#endif /* DEBUG */

            /* advance the lj when mj hits the top */
            if (mj == lj) lj++;

        } /* for j */

        /* advance the li when mi hits the top */
        if (mi == li) li++;

    } /* for i */
#ifdef DEBUG
    printf(
        "\n");
    fflush(stdout);
#endif /* DEBUG */

    return (hamiltonian);
}

/* check the wavefunction for g or u symmetry */
int determine_rotational_eigensymmetry(molecule_t *molecule, int level, int l_max) {
    int symmetry;
    double theta, phi;
    int i, l, m, index;
    complex_t wavefunction[QUANTUM_ROTATION_SYMMETRY_POINTS];
    double sqmod[QUANTUM_ROTATION_SYMMETRY_POINTS];
    double max_sqmod, max_theta, max_phi;
    complex_t max_wavefunction, inv_wavefunction;

    /* scan a few random points, pick the one with the largest square of the real part */
    for (i = 0, max_sqmod = 0; i < QUANTUM_ROTATION_SYMMETRY_POINTS; i++) {
        /* get random theta and phi */
        theta = get_rand(system) * M_PI;
        phi = get_rand(system) * 2.0 * M_PI;

        /* project the eigenvector into the basis */
        wavefunction[i].real = 0;
        for (l = 0, index = 0; l <= l_max; l++)
            for (m = -l; m <= l; m++, index++)
                wavefunction[i].real += molecule->quantum_rotational_eigenvectors[level][index].real * rotational_basis(REAL, l, m, theta, phi);

        /* get the square modulus */
        sqmod[i] = wavefunction[i].real * wavefunction[i].real;

        /* if we have a new max, save the domain point */
        if (sqmod[i] > max_sqmod) {
            max_sqmod = sqmod[i];
            max_theta = theta;
            max_phi = phi;
        }
    }

    /* check the symmetry at the maximum */
    max_wavefunction.real = 0;
    inv_wavefunction.real = 0;
    for (l = 0, index = 0; l <= l_max; l++) {
        for (m = -l; m <= l; m++, index++) {
            /* get the value of the wavefunction at the maximum */
            max_wavefunction.real += molecule->quantum_rotational_eigenvectors[level][index].real * rotational_basis(REAL, l, m, max_theta, max_phi);

            /* get the value of it's inversion */
            inv_wavefunction.real += molecule->quantum_rotational_eigenvectors[level][index].real * rotational_basis(REAL, l, m, (max_theta + M_PI), (max_phi + M_PI));
        }
    }

    /* did we change sign? */
    if (max_wavefunction.real * inv_wavefunction.real < 0.0)
        symmetry = QUANTUM_ROTATION_ANTISYMMETRIC;
    else
        symmetry = QUANTUM_ROTATION_SYMMETRIC;

    return (symmetry);
}

/* get the rotational energy levels of a single rotor in an external potential */
void quantum_rotational_energies(system_t *system, molecule_t *molecule, int level_max, int l_max) {
    complex_t **hamiltonian;
    int i, j, dim, p;
    /* variables needed for zhpevx() */
    char jobz, range, uplo;
    int vl, vu, il, iu, m, ldz, *iwork, *ifail, info;
    double *hamiltonian_packed, *eigenvalues, abstol, *z, *work, *rwork;

    /* determine the size of the symmetric matrix */
    dim = (l_max + 1) * (l_max + 1);

    /* setup the lapack arguments */
    hamiltonian_packed = calloc((int)(dim * (dim + 1.0) / 2.0), sizeof(complex_t));
    memnullcheck(hamiltonian_packed, (int)(dim * (dim + 1.0) / 2.0) * sizeof(complex_t), __LINE__ - 1, __FILE__ - 1);
    eigenvalues = calloc(dim, sizeof(double));
    memnullcheck(eigenvalues, dim * sizeof(double), __LINE__ - 1, __FILE__ - 1);
    z = calloc(dim * dim, sizeof(complex_t));
    memnullcheck(z, dim * dim * sizeof(complex_t), __LINE__ - 1, __FILE__ - 1);
    work = calloc((2 * dim), sizeof(complex_t));
    memnullcheck(work, 2 * dim * sizeof(complex_t), __LINE__ - 1, __FILE__ - 1);
    rwork = calloc((7 * dim), sizeof(double));
    memnullcheck(rwork, 7 * dim * sizeof(double), __LINE__ - 1, __FILE__ - 1);
    iwork = calloc((5 * dim), sizeof(int));
    memnullcheck(iwork, 5 * dim * sizeof(int), __LINE__ - 1, __FILE__ - 1);
    ifail = calloc(dim, sizeof(int));
    memnullcheck(ifail, dim * sizeof(int), __LINE__ - 1, __FILE__ - 1);

    jobz = 'V';
    range = 'A';
    uplo = 'U';
    vl = 0;
    vu = 0;
    il = 0;
    iu = 0;
    abstol = 1e-300;
    m = 0;
    ldz = dim;
    info = 0;

    /* build the rotational hamiltonian */
    hamiltonian = rotational_hamiltonian(system, molecule, level_max, l_max);

    /* pack the lapack array */
    for (j = 0, p = 0; j < dim; j++) {
        for (i = 0; i <= j; i++, p += 2) {
            hamiltonian_packed[p + 0] = hamiltonian[i][j].real;
            hamiltonian_packed[p + 1] = hamiltonian[i][j].imaginary;
        }
    }

/* diagonalize the hamiltonian */
#ifdef ACML_NOUNDERSCORE
    zhpevx(&jobz, &range, &uplo, &dim, hamiltonian_packed, &vl, &vu, &il, &iu, &abstol, &m, eigenvalues, z, &ldz, work, rwork, iwork, ifail, &info);
#else
    zhpevx_(&jobz, &range, &uplo, &dim, hamiltonian_packed, &vl, &vu, &il, &iu, &abstol, &m, eigenvalues, z, &ldz, work, rwork, iwork, ifail, &info);
#endif /* ACML_NOUNDERSCORE */

    /* store the energy levels  */
    for (i = 0; i < level_max; i++) molecule->quantum_rotational_energies[i] = eigenvalues[i];

    /* store the eigenvectors */
    for (i = 0, p = 0; i < level_max; i++) {
        for (j = 0; j < dim; j++, p += 2) {
            molecule->quantum_rotational_eigenvectors[i][j].real = z[p];
            molecule->quantum_rotational_eigenvectors[i][j].imaginary = z[p + 1];
        }
    }

    /* get the symmetry of each eigenvector */
    for (i = 0; i < level_max; i++)
        molecule->quantum_rotational_eigensymmetry[i] = determine_rotational_eigensymmetry(molecule, i, l_max);

    /* free our arrays */
    for (i = 0; i < dim; i++) free(hamiltonian[i]);
    free(hamiltonian);
    free(hamiltonian_packed);
    free(eigenvalues);
    free(z);
    free(work);
    free(rwork);
    free(iwork);
    free(ifail);
}

/* generate the potential over the quadrature grid */
/* Gauss-Legendre abscissas+weights courtesy of Abramowitz & Stegun, "Handbook of Mathematical Functions", 9th Ed., p.916 */
void quantum_rotational_grid(system_t *system, molecule_t *molecule) {
    int t, p;
    double theta, phi;
    double potential;
    /* N = 16 */
    double roots[QUANTUM_ROTATION_GRID] = {-0.989400934991649932596,
                                           -0.944575023073232576078,
                                           -0.865631202387831743880,
                                           -0.755404408355003033895,
                                           -0.617876244402643748447,
                                           -0.458016777657227386342,
                                           -0.281603550779258913230,
                                           -0.095012509837637440185,
                                           0.095012509837637440185,
                                           0.281603550779258913230,
                                           0.458016777657227386342,
                                           0.617876244402643748447,
                                           0.755404408355003033895,
                                           0.865631202387831743880,
                                           0.944575023073232576078,
                                           0.989400934991649932596};

#ifdef XXX
    /* N = 8 */
    double roots[QUANTUM_ROTATION_GRID] = {-0.960289856497536, -0.796666477413627, -0.525532409916329, -0.183434642495650,
                                           0.183434642495650, 0.525532409916329, 0.796666477413627, 0.960289856497536};
    /* N = 32 */
    double roots[QUANTUM_ROTATION_GRID] = {-0.997263861849481563545, -0.985611511545268335400, -0.964762255587506430774, -0.934906075937739689171,
                                           -0.896321155766052123965, -0.849367613732569970134, -0.794483795967942406963, -0.732182118740289680387,
                                           -0.663044266930215200975, -0.587715757240762329041, -0.506899908932229390024, -0.421351276130635345364,
                                           -0.331868602282127649780, -0.239287362252137074545, -0.144471961582796493485, -0.048307665687738316235,
                                           0.048307665687738316235, 0.144471961582796493485, 0.239287362252137074545, 0.331868602282127649780,
                                           0.421351276130635345364, 0.506899908932229390024, 0.587715757240762329041, 0.663044266930215200975,
                                           0.732182118740289680387, 0.794483795967942406963, 0.849367613732569970134, 0.896321155766052123965,
                                           0.934906075937739689171, 0.964762255587506430774, 0.985611511545268335400, 0.997263861849481563545};
#endif /* XXX */

    /* get the potential on an angular grid */
    for (p = 0; p < QUANTUM_ROTATION_GRID; p++) {
        phi = M_PI * roots[p] + M_PI;
        for (t = 0; t < QUANTUM_ROTATION_GRID; t++) {
            theta = 0.5 * M_PI * roots[t] + 0.5 * M_PI;
            molecule->quantum_rotational_potential_grid[t][p] = sin(theta) * rotational_potential(system, molecule, theta, phi);

        } /* for t */

    } /* for p */
}

/* find the rotational 1-body energies for each molecule in the system */
void quantum_system_rotational_energies(system_t *system) {
    int i, j;
    molecule_t *molecule_ptr;

    /* get the rotational eigenspectrum for each moveable molecule */
    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        if (!(molecule_ptr->frozen || molecule_ptr->adiabatic)) {
            /* generate the necessary grids */
            quantum_rotational_grid(system, molecule_ptr);

            /* get the eigenspectrum */
            quantum_rotational_energies(system, molecule_ptr, system->quantum_rotation_level_max, system->quantum_rotation_l_max);

            /* sum over states into molecular rotational partition functions */
            molecule_ptr->rot_partfunc_g = 0;
            molecule_ptr->rot_partfunc_u = 0;
            for (i = 0; i < system->quantum_rotation_sum; i++) {
                /* symmetry-specific part funcs */
                /* factor of 3 for ungerand comes from degeneracy of nuclear wavefunction */
                if (molecule_ptr->quantum_rotational_eigensymmetry[i] == QUANTUM_ROTATION_SYMMETRIC)
                    molecule_ptr->rot_partfunc_g += exp(-molecule_ptr->quantum_rotational_energies[i] / system->temperature);
                else
                    molecule_ptr->rot_partfunc_u += 3. * exp(-molecule_ptr->quantum_rotational_energies[i] / system->temperature);
            }

            /* assign partition function ratio based upon spin symmetry */
            if (molecule_ptr->nuclear_spin == NUCLEAR_SPIN_PARA) {
                molecule_ptr->rot_partfunc = molecule_ptr->rot_partfunc_g / (molecule_ptr->rot_partfunc_g + molecule_ptr->rot_partfunc_u);
            } else {
                molecule_ptr->rot_partfunc = molecule_ptr->rot_partfunc_u / (molecule_ptr->rot_partfunc_g + molecule_ptr->rot_partfunc_u);
            }
        }
    }

#ifdef DEBUG
    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        if (!(molecule_ptr->frozen || molecule_ptr->adiabatic)) {
            for (i = 0; i < system->quantum_rotation_level_max; i++) {
                printf(
                    "DEBUG_QROT: molecule #%d (%s) rotational level %d = %.6f K (%.6f cm^-1 or %.6f meV or %.6f / B) ",
                    molecule_ptr->id, molecule_ptr->moleculetype, i, molecule_ptr->quantum_rotational_energies[i],
                    molecule_ptr->quantum_rotational_energies[i] * KB / (100.0 * H * C),
                    8.61733238e-2 * (molecule_ptr->quantum_rotational_energies[i] - molecule_ptr->quantum_rotational_energies[0]),
                    molecule_ptr->quantum_rotational_energies[i] / system->quantum_rotation_B);

                if (molecule_ptr->quantum_rotational_eigensymmetry[i] == QUANTUM_ROTATION_SYMMETRIC)
                    printf(
                        "*** symmetric ***");
                else if (molecule_ptr->quantum_rotational_eigensymmetry[i] == QUANTUM_ROTATION_ANTISYMMETRIC)
                    printf(
                        "*** antisymmetric ***");

                printf(
                    "\n");
            }

            for (i = 0; i < system->quantum_rotation_level_max; i++) {
                printf(
                    "DEBUG_QROT: molecule #%d (%s) eigenvec rot level %d\n", molecule_ptr->id, molecule_ptr->moleculetype, i);
                for (j = 0; j < (system->quantum_rotation_l_max + 1) * (system->quantum_rotation_l_max + 1); j++)
                    printf(
                        "\tj=%d %.6f %.6f\n", j, molecule_ptr->quantum_rotational_eigenvectors[i][j].real,
                        molecule_ptr->quantum_rotational_eigenvectors[i][j].imaginary);
                printf(
                    "\n");
            }
        }
    }

#endif /* DEBUG */
}
