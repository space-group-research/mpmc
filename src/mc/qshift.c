#include <stdio.h>
#include <stdlib.h>
#include <mc.h>

//perturb the position of H2Q sites, and recalculate the
// charges on H2Q and H2G sites
void qshift_do(
    system_t *system,
    qshiftData_t *q,
    double scale_r,
    double Quadrupole) {
    double xH2Qold = 0, xH2Qnew;
    int FLAG = 0;
    atom_t *atom_ptr;

    for (atom_ptr = system->molecules->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
        if (!strncasecmp(atom_ptr->atomtype,
                         "H2Q", 3)) {
            xH2Qold = atom_ptr->pos[0];
            FLAG = 1;
        }
    }

    if (FLAG == 0) {
        error(
            "SURFACE: QSHIFT ERROR: NOT ATOM WITH LABEL H2Q. TURN OFF QSHIFT!\n");
        die(-1);
    }
    if (xH2Qold == 0) {
        error(
            "SURFACE: QSHIFT ERROR: H2Q HAS POSITION->X == 0.\n");
        error(
            "ALIGN YOUR LINEAR MOLECULE ALONG THE X-AXIS\n");
        die(-1);
    }

    // use surf_scale_r and get_rand to adjust position of H2Q
    q->drH2Q = scale_r * (0.5 - get_rand(system));
    xH2Qnew = xH2Qold - q->drH2Q;

    // set qH2Q to conserve quadrupole
    q->qH2Q = 0.5 * Quadrupole / (xH2Qnew * xH2Qnew);

    // set qH2G to conserve charge
    q->qH2G = -2.0 * q->qH2Q;

    return;
}

//take system variables, return quadrupole
double calcquadrupole(system_t *system) {
    double x = 0, q = 0;
    int FLAG = 0;

    atom_t *atom_ptr;

    for (atom_ptr = system->molecules->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
        if (!strncasecmp(atom_ptr->atomtype,
                         "H2Q", 3)) {
            x = atom_ptr->pos[0];
            q = atom_ptr->charge;
            FLAG = 1;
        }
    }

    if (FLAG == 0) {
        error(
            "SURFACE: QSHIFT ERROR: NOT ATOM WITH LABEL H2Q. TURN OFF QSHIFT!\n");
        die(-1);
    }
    if (q == 0) {
        error(
            "SURFACE: QSHIFT ERROR: H2Q HAS CHARGE == 0."
            " WHAT THE HECK?\n");
        die(-1);
    }
    if (x == 0) {
        error(
            "SURFACE: QSHIFT ERROR: H2Q HAS POSITION->X == 0.");
        error(
            " ALIGN YOUR LINEAR MOLECULE ALONG THE X-AXIS\n");
        die(-1);
    }

    return 2.0 * x * x * q;
}
