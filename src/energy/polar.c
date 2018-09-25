/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

/* RRMS of dipoles */
double get_dipole_rrms(system_t *system) {
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    double N, dipole_rrms;

    dipole_rrms = N = 0;

    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            if (isfinite(atom_ptr->dipole_rrms))
                dipole_rrms += atom_ptr->dipole_rrms;
            N++;
        }
    }
    return dipole_rrms / N;
}

/* get the induction energy */
double polar(system_t *system) {
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    int num_iterations;
    double potential;
    char linebuf[MAXLINE];

    /* take measures to let N fluctuate */
    if ((system->ensemble == ENSEMBLE_UVT || system->ensemble == ENSEMBLE_REPLAY) && !system->polar_zodid)
        thole_resize_matrices(system);

    /* get the A matrix */
    if (!system->polar_zodid) {
        thole_amatrix(system);
        if (system->polarizability_tensor) {
            output(
                "POLAR: A matrix:\n");
            print_matrix(3 * ((int)system->checkpoint->thole_N_atom), system->A_matrix);
        }
    }

    /* find the dipoles */

    if (system->polar_ewald_full) {
        //do a full-ewald polarization treatment
        ewald_full(system);

    } else if (system->polar_iterative) {
        //solve the self-consistent problem
        thole_field(system);                       //calc e-field
        num_iterations = thole_iterative(system);  //calc dipoles

        system->nodestats->polarization_iterations = (double)num_iterations;  //statistics
        system->observables->dipole_rrms = get_dipole_rrms(system);

        if (system->iter_success) {
            switch (system->ensemble) {
                case ENSEMBLE_UVT:
                case ENSEMBLE_NVT:
                case ENSEMBLE_NVE:
                case ENSEMBLE_NPT:
                    sprintf(linebuf,
                            "POLAR: polarization iterative solver convergence failure on mc step %d.\n", system->step);
                    error(linebuf);
                    break;
                case ENSEMBLE_REPLAY:
                    sprintf(linebuf,
                            "POLAR: polarization iterative solver convergence failure on configuration %d.\n", system->step);
                    error(linebuf);
                    break;
                case ENSEMBLE_SURF:
                case ENSEMBLE_SURF_FIT:
                case ENSEMBLE_TE:
                default:
                    sprintf(linebuf,
                            "POLAR: polarization iterative solver failed to reach convergence.\n");
                    error(linebuf);
            }
        }

    } else {
        //do matrix inversion
        thole_field(system);            //calc e-field
        thole_bmatrix(system);          //matrix inversion
        thole_bmatrix_dipoles(system);  //get dipoles

        /* output the 3x3 molecular polarizability tensor */
        if (system->polarizability_tensor) {
            output(
                "POLAR: B matrix:\n");
            print_matrix(3 * ((int)system->checkpoint->thole_N_atom), system->B_matrix);
            thole_polarizability_tensor(system);
            die(0);
        }
    }

    /* calculate the polarization energy as 1/2 mu*E */
    potential = 0;
    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            potential += dddotprod(atom_ptr->mu, atom_ptr->ef_static);
            if (system->polar_palmo)
                potential += dddotprod(atom_ptr->mu, atom_ptr->ef_induced_change);
        }
    }
    potential *= -0.5;

#ifdef DEBUG
    fprintf(stderr,
            "mu MOLECULE ATOM * DIPOLES * STATIC * INDUCED * pot/atom -0.5*mu*E_s\n");
    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            fprintf(stderr,
                    "mu %4d %4d * %8.5lf %8.5lf %8.5lf * %8.5lf %8.5lf %8.5lf * %8.5lf %8.5lf %8.5lf * %lf %lf\n",
                    molecule_ptr->id, atom_ptr->id,
                    atom_ptr->mu[0], atom_ptr->mu[1], atom_ptr->mu[2],
                    atom_ptr->ef_static[0], atom_ptr->ef_static[1], atom_ptr->ef_static[2],
                    atom_ptr->ef_induced[0], atom_ptr->ef_induced[1], atom_ptr->ef_induced[2],
                    potential / system->natoms, -0.5 * atom_ptr->mu[0] * atom_ptr->ef_static[0]);
        }
    }
#endif

    return (potential);
}
