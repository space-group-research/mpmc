#include <mc.h>

//Copyright 2013-2019 Adam Hogan

double de_fx(double r, double b, double sigma, double c6, double c8, double c10, int damp) {
    double energy = 0.0, repulsion = 0.0;

    double r2 = r * r;
    double r4 = r2 * r2;
    double r6 = r4 * r2;
    double r8 = r6 * r2;
    double r10 = r8 * r2;

    if (b != 0.0 && sigma != 0.0)
        repulsion = 596.725194095 * 1.0 / b * exp(-b * (r - sigma));

    if (damp)
        energy = -tt_damping(6, b * r) * c6 / r6 - tt_damping(8, b * r) * c8 / r8 - tt_damping(10, b * r) * c10 / r10 + repulsion;
    else
        energy = -c6 / r6 - c8 / r8 - c10 / r10 + repulsion;

    return energy;
}

double disp_expansion_fh_corr(system_t *system, molecule_t *molecule_ptr, pair_t *pair_ptr, int order) {
    double reduced_mass;
    double dE, d2E, d3E, d4E;  //energy derivatives
    double corr;
    double ir = 1.0 / pair_ptr->rimg;
    double ir2 = ir * ir;
    double ir3 = ir2 * ir;

    if ((order != 2) && (order != 4)) return NAN;  //must be order 2 or 4

    reduced_mass = AMU2KG * molecule_ptr->mass * pair_ptr->molecule->mass /
                   (molecule_ptr->mass + pair_ptr->molecule->mass);

    // I don't really feel like doing the analytical derivatives for this because there will be a million terms in the damping function
    double h = 0.001;                                                                                                                                // this seems to work fine and have atleast an order of magnitude on both sides
    double pm1 = de_fx(pair_ptr->rimg - h, pair_ptr->epsilon, pair_ptr->sigma, pair_ptr->c6, pair_ptr->c8, pair_ptr->c10, system->damp_dispersion);  // back one
    double p0 = de_fx(pair_ptr->rimg, pair_ptr->epsilon, pair_ptr->sigma, pair_ptr->c6, pair_ptr->c8, pair_ptr->c10, system->damp_dispersion);       // center
    double pp1 = de_fx(pair_ptr->rimg + h, pair_ptr->epsilon, pair_ptr->sigma, pair_ptr->c6, pair_ptr->c8, pair_ptr->c10, system->damp_dispersion);  // forward one

    dE = (-0.5 * pm1 + 0.5 * pp1) / h;
    d2E = (pm1 - 2.0 * p0 + pp1) / (h * h);

    //2nd order correction
    corr = M2A2 *
           (HBAR2 / (24.0 * KB * system->temperature * reduced_mass)) *
           (d2E + 2.0 * dE / pair_ptr->rimg);

    if (order >= 4) {
        double pm2 = de_fx(pair_ptr->rimg - 2.0 * h, pair_ptr->epsilon, pair_ptr->sigma, pair_ptr->c6, pair_ptr->c8, pair_ptr->c10, system->damp_dispersion);  // back two
        double pp2 = de_fx(pair_ptr->rimg + 2.0 * h, pair_ptr->epsilon, pair_ptr->sigma, pair_ptr->c6, pair_ptr->c8, pair_ptr->c10, system->damp_dispersion);  // forward two

        d3E = (-0.5 * pm2 + pm1 - pp1 + 0.5 * pp2) / (h * h * h);
        d4E = (pm2 - 4.0 * pm1 + 6.0 * p0 - 4.0 * pp1 + pp2) / (h * h * h * h);

        //4th order corection
        corr += M2A4 *
                (HBAR4 / (1152.0 * KB2 * system->temperature * system->temperature * reduced_mass * reduced_mass)) *
                (15.0 * dE * ir3 + 4.0 * d3E * ir + d4E);
    }

    return corr;
}

double disp_expansion_lrc(const system_t *system, pair_t *pair_ptr, const double cutoff) /* ignoring the exponential repulsion bit because it decays exponentially */
{
    if (!(pair_ptr->frozen) &&                                                      /* disqualify frozen pairs */
        ((pair_ptr->lrc == 0.0) || pair_ptr->last_volume != system->pbc->volume)) { /* LRC only changes if the volume change */

        pair_ptr->last_volume = system->pbc->volume;

        return -4.0 * M_PI * (pair_ptr->c6 / (3.0 * cutoff * cutoff * cutoff) + pair_ptr->c8 / (5.0 * cutoff * cutoff * cutoff * cutoff * cutoff) + pair_ptr->c10 / (7.0 * cutoff * cutoff * cutoff * cutoff * cutoff * cutoff * cutoff)) / system->pbc->volume;
    }

    else
        return pair_ptr->lrc; /* use stored value */
}

double disp_expansion_lrc_self(const system_t *system, atom_t *atom_ptr, const double cutoff) {
    if (!(atom_ptr->frozen) &&                                                           /* disqualify frozen atoms */
        ((atom_ptr->lrc_self == 0.0) || atom_ptr->last_volume != system->pbc->volume)) { /* LRC only changes if the volume change) */

        atom_ptr->last_volume = system->pbc->volume;

        if (system->extrapolate_disp_coeffs) {
            double c10;
            if (atom_ptr->c6 != 0.0 && atom_ptr->c8 != 0.0)
                c10 = 49.0 / 40.0 * atom_ptr->c8 * atom_ptr->c8 / atom_ptr->c6;
            else
                c10 = 0.0;

            return -4.0 * M_PI * (atom_ptr->c6 / (3.0 * cutoff * cutoff * cutoff) + atom_ptr->c8 / (5.0 * cutoff * cutoff * cutoff * cutoff * cutoff) + c10 / (7.0 * cutoff * cutoff * cutoff * cutoff * cutoff * cutoff * cutoff)) / system->pbc->volume;
        } else
            return -4.0 * M_PI * (atom_ptr->c6 / (3.0 * cutoff * cutoff * cutoff) + atom_ptr->c8 / (5.0 * cutoff * cutoff * cutoff * cutoff * cutoff) + atom_ptr->c10 / (7.0 * cutoff * cutoff * cutoff * cutoff * cutoff * cutoff * cutoff)) / system->pbc->volume;
    }

    return atom_ptr->lrc_self; /* use stored value */
}

double disp_expansion(system_t *system) {
    double potential = 0.0;

    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr;

    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
                if (pair_ptr->recalculate_energy) {
                    /* pair LRC */
                    if (system->rd_lrc)
                        pair_ptr->lrc = disp_expansion_lrc(system, pair_ptr, system->pbc->cutoff);

                    /* make sure we're not excluded or beyond the cutoff */
                    if (!(pair_ptr->rd_excluded || pair_ptr->frozen)) {
                        const double r = pair_ptr->rimg;
                        const double r2 = r * r;
                        const double r4 = r2 * r2;
                        const double r6 = r4 * r2;
                        const double r8 = r6 * r2;
                        const double r10 = r8 * r2;

                        double c6 = pair_ptr->c6;
                        const double c8 = pair_ptr->c8;
                        const double c10 = pair_ptr->c10;

                        if (system->disp_expansion_mbvdw == 1)
                            c6 = 0.0;

                        double repulsion = 0.0;

                        // F0 = 0.001 Eh/bohr aka a.u.
                        // .001/3.166811429E-6*1.8897161646321 = 596.725194095
                        if (pair_ptr->epsilon != 0.0 && pair_ptr->sigma != 0.0)
                            repulsion = 596.725194095 * 1.0 / pair_ptr->epsilon * exp(-pair_ptr->epsilon * (r - pair_ptr->sigma));

                        if (system->damp_dispersion)
                            pair_ptr->rd_energy = -tt_damping(6, pair_ptr->epsilon * r) * c6 / r6 - tt_damping(8, pair_ptr->epsilon * r) * c8 / r8 - tt_damping(10, pair_ptr->epsilon * r) * c10 / r10 + repulsion;
                        else
                            pair_ptr->rd_energy = -c6 / r6 - c8 / r8 - c10 / r10 + repulsion;

                        if (system->feynman_hibbs)
                            pair_ptr->rd_energy += disp_expansion_fh_corr(system, molecule_ptr, pair_ptr, system->feynman_hibbs_order);

                        if (system->cavity_autoreject_repulsion != 0.0)
                            if (repulsion > system->cavity_autoreject_repulsion)
                                pair_ptr->rd_energy = MAXVALUE;
                    }
                }
                potential += pair_ptr->rd_energy + pair_ptr->lrc;
            }
        }
    }

    if (system->disp_expansion_mbvdw == 1) {
        thole_amatrix(system);
        potential += vdw(system);
    }

    /* calculate self LRC interaction */
    if (system->rd_lrc) {
        for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
            for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
                atom_ptr->lrc_self = disp_expansion_lrc_self(system, atom_ptr, system->pbc->cutoff);
                potential += atom_ptr->lrc_self;
            }
        }
    }

    return potential;
}

double disp_expansion_nopbc(system_t *system) {
    double potential = 0.0;

    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr;

    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
                if (pair_ptr->recalculate_energy) {
                    /* make sure we're not excluded or beyond the cutoff */
                    if (!(pair_ptr->rd_excluded || pair_ptr->frozen)) {
                        const double r = pair_ptr->rimg;
                        const double r2 = r * r;
                        const double r4 = r2 * r2;
                        const double r6 = r4 * r2;
                        const double r8 = r6 * r2;
                        const double r10 = r8 * r2;

                        double c6 = pair_ptr->c6;
                        const double c8 = pair_ptr->c8;
                        const double c10 = pair_ptr->c10;

                        if (system->disp_expansion_mbvdw == 1)
                            c6 = 0.0;

                        double repulsion = 0.0;

                        if (pair_ptr->epsilon != 0.0 && pair_ptr->sigma != 0.0)
                            repulsion = 596.725194095 * 1.0 / pair_ptr->epsilon * exp(-pair_ptr->epsilon * (r - pair_ptr->sigma));

                        if (system->damp_dispersion)
                            pair_ptr->rd_energy = -tt_damping(6, pair_ptr->epsilon * r) * c6 / r6 - tt_damping(8, pair_ptr->epsilon * r) * c8 / r8 - tt_damping(10, pair_ptr->epsilon * r) * c10 / r10 + repulsion;
                        else
                            pair_ptr->rd_energy = -c6 / r6 - c8 / r8 - c10 / r10 + repulsion;
                    }
                }
                potential += pair_ptr->rd_energy;
            }
        }
    }

    if (system->disp_expansion_mbvdw == 1) {
        thole_amatrix(system);
        potential += vdw(system);
    }

    return potential;
}

double factorial(int n) {
    int i;
    double fac = 1.0;
    for (i = 2; i <= n; i++)
        fac *= i;
    return fac;
}

double tt_damping(int n, double br) {
    double sum = 1.0, running_br = br;
    int i;
    for (i = 1; i <= n; i++) {
        sum += running_br / factorial(i);
        running_br *= br;
    }

    const double result = 1.0 - exp(-br) * sum;

    if (result > 0.000000001)
        return result;
    else
        return 0.0; /* This is so close to zero lets just call it zero to avoid rounding error and the simulation blowing up */
}
