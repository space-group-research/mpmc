/* 

Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>
#define OneOverSqrtPi 0.56418958354

//called from energy/polar.c
/* calculate the field with periodic boundaries */
void thole_field(system_t *system) {
    int p;
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;

    /* zero the field vectors */
    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            for (p = 0; p < 3; p++) {
                atom_ptr->ef_static[p] = 0;
                atom_ptr->ef_static_self[p] = 0;
            }
        }
    }

    /* calculate the electrostatic field */
    if (system->polar_ewald)
        ewald_estatic(system);
    else if (system->polar_wolf || system->polar_wolf_full)
        thole_field_wolf(system);
    else
        thole_field_nopbc(system);
}

/* calculate the field without ewald summation/wolf */
void thole_field_nopbc(system_t *system) {
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr;
    int p;
    double r;

    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
                if (pair_ptr->frozen) continue;
                if (molecule_ptr == pair_ptr->molecule) continue;  //don't let molecules polarize themselves

                r = pair_ptr->rimg;

                //inclusive near the cutoff
                if ((r - SMALL_dR < system->pbc->cutoff) && (r != 0.)) {
                    for (p = 0; p < 3; p++) {
                        atom_ptr->ef_static[p] += pair_ptr->atom->charge * pair_ptr->dimg[p] / (r * r * r);
                        pair_ptr->atom->ef_static[p] -= atom_ptr->charge * pair_ptr->dimg[p] / (r * r * r);
                    }

                } /* cutoff */

            } /* pair */
        }     /* atom */
    }         /* molecule */

    return;
}

// calc field using wolf sum (JCP 124 234104 (2006) equation 19
void thole_field_wolf(system_t *system) {
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr;
    int p;         //dimensionality
    double r, rr;  //r and 1/r (reciprocal of r)
    double R = system->pbc->cutoff;
    double rR = 1. / R;
    //used for polar_wolf_alpha (aka polar_wolf_damp)
    double a = system->polar_wolf_alpha;
    double erR;  //complementary error functions
    erR = erfc(a * R);
    double cutoffterm = (erR * rR * rR + 2.0 * a * OneOverSqrtPi * exp(-a * a * R * R) * rR);
    double bigmess = 0;

    //init lookup table if needed
    if (system->polar_wolf_alpha_lookup && !(system->polar_wolf_alpha_table))
        system->polar_wolf_alpha_table = polar_wolf_alpha_lookup_init(system);

    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
                if (molecule_ptr == pair_ptr->molecule) continue;  //don't let molecules polarize themselves
                if (pair_ptr->frozen) continue;                    //don't let the MOF polarize itself

                r = pair_ptr->rimg;

                if ((r - SMALL_dR < system->pbc->cutoff) && (r != 0.)) {
                    rr = 1. / r;

                    //we will need this shit if wolf alpha != 0
                    if ((a != 0) & system->polar_wolf_alpha_lookup)
                        bigmess = polar_wolf_alpha_getval(system, r);
                    else if (a != 0)  //no lookup
                        bigmess = (erfc(a * r) * rr * rr + 2.0 * a * OneOverSqrtPi * exp(-a * a * r * r) * rr);

                    for (p = 0; p < 3; p++) {
                        //see JCP 124 (234104)
                        if (a == 0) {
                            atom_ptr->ef_static[p] += (pair_ptr->atom->charge) * (rr * rr - rR * rR) * pair_ptr->dimg[p] * rr;
                            pair_ptr->atom->ef_static[p] -= (atom_ptr->charge) * (rr * rr - rR * rR) * pair_ptr->dimg[p] * rr;
                        } else {
                            atom_ptr->ef_static[p] += pair_ptr->atom->charge * (bigmess - cutoffterm) * pair_ptr->dimg[p] * rr;
                            pair_ptr->atom->ef_static[p] -= atom_ptr->charge * (bigmess - cutoffterm) * pair_ptr->dimg[p] * rr;
                        }
                    }

                } /* no lookup table */
            }     /* pair */
        }         /* atom */
    }             /* molecule */

    return;
}
