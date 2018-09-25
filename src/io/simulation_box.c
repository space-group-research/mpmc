#include <mc.h>

void printboxdim(pbc_t* pbc) {
    char buffer[MAXLINE];

    sprintf(buffer,
            "INPUT: unit cell volume = %.3f A^3 (cutoff = %.3f A)\n", pbc->volume, pbc->cutoff);
    output(buffer);

    sprintf(buffer,
            "INPUT: Basis vectors { a = %8.3lf \tb = %8.3lf \tc = %8.3lf }\n",
            sqrt(dddotprod(pbc->basis[0], pbc->basis[0])),
            sqrt(dddotprod(pbc->basis[1], pbc->basis[1])),
            sqrt(dddotprod(pbc->basis[2], pbc->basis[2])));
    output(buffer);

    sprintf(buffer,
            "INPUT: Basis angles  { α = %8.3lf \tβ = %8.3lf \tγ = %8.3lf }\n",
            180.0 / M_PI * acos(dddotprod(pbc->basis[1], pbc->basis[2]) / sqrt(dddotprod(pbc->basis[1], pbc->basis[1]) * dddotprod(pbc->basis[2], pbc->basis[2]))),
            180.0 / M_PI * acos(dddotprod(pbc->basis[2], pbc->basis[0]) / sqrt(dddotprod(pbc->basis[2], pbc->basis[2]) * dddotprod(pbc->basis[0], pbc->basis[0]))),
            180.0 / M_PI * acos(dddotprod(pbc->basis[0], pbc->basis[1]) / sqrt(dddotprod(pbc->basis[0], pbc->basis[0]) * dddotprod(pbc->basis[1], pbc->basis[1]))));
    output(buffer);

    return;
}

int setup_simulation_box(FILE* finput, system_t* system) {
    // read in input.pqr molecules
    system->molecules = read_molecules(finput, system);
    if (!system->molecules && system->ensemble == ENSEMBLE_REPLAY) {
        output(
            "INPUT: end of trajectory file\n");
        return 1;
    } else if (!system->molecules) {
        error(
            "INPUT: error reading in input molecules\n");
        return (-1);
    } else
        output(
            "INPUT: finished reading in molecules\n");

    //read in pqr box and calculate periodic boundary conditions if neccessary
    if (system->read_pqr_box_on)
        read_pqr_box(finput, system);

    pbc(system);
    if (system->ensemble != ENSEMBLE_SURF && system->ensemble != ENSEMBLE_SURF_FIT) {
        if ((system->pbc->volume <= 0.0) || (system->pbc->cutoff <= 0.0)) {
            error(
                "INPUT: invalid simulation box dimensions.\n");
            return -1;
            ;
        }
    }
    if (system->pbc->volume > 0)
        printboxdim(system->pbc);

    // read in the insertion molecules
    if (system->insert_input) {
        system->insertion_molecules = read_insertion_molecules(system);
        if (!system->insertion_molecules) {
            error(
                "INPUT: error read in insertion molecules\n");
            return -1;
        } else
            output(
                "INPUT: finished reading in insertion molecules\n");
    } else  //else only 1 sorbate type
        system->sorbateCount = 1;

    // now that we've read in the sorbates, we can check that user_fugacities is properly set (if used)
    if (system->user_fugacities) {
        if (system->fugacitiesCount != system->sorbateCount) {
            error(
                "INPUT: number of fugacities set via user_fugacities does not match the number of sorbates.\n");
            return -1;
        }
    }

    return 0;
}
