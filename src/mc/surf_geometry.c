/*

Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>
#include <quaternion.h>

/* set the distance and angles for a particular dimer orientation */
int surface_dimer_geometry(system_t *system, double r, double alpha_origin, double beta_origin, double gamma_origin, double alpha_move, double beta_move, double gamma_move, int reverse_flag) {
    int i;
    molecule_t *molecule_origin, *molecule_move;
    atom_t *atom_ptr;

    /* make sure that there are only two molecules */
    molecule_origin = system->molecules;
    molecule_move = molecule_origin->next;
    if (!molecule_move) {
        error(
            "SURFACE: the input PQR has only a single molecule\n");
        return (-1);
    }
    if (molecule_move->next) {
        error(
            "SURFACE: the input PQR must contain exactly two molecules\n");
        return (-1);
    }

    /* relocate both molecules to the origin */
    for (atom_ptr = molecule_origin->atoms; atom_ptr; atom_ptr = atom_ptr->next)
        for (i = 0; i < 3; i++)
            atom_ptr->pos[i] -= molecule_origin->com[i];
    for (i = 0; i < 3; i++)
        molecule_origin->com[i] = 0;

    for (atom_ptr = molecule_move->atoms; atom_ptr; atom_ptr = atom_ptr->next)
        for (i = 0; i < 3; i++)
            atom_ptr->pos[i] -= molecule_move->com[i];
    for (i = 0; i < 3; i++)
        molecule_move->com[i] = 0;

    /* relocate the moveable molecule r distance away along the x axis */
    molecule_move->com[0] = r;
    molecule_move->com[1] = 0;
    molecule_move->com[2] = 0;
    for (atom_ptr = molecule_move->atoms; atom_ptr; atom_ptr = atom_ptr->next)
        for (i = 0; i < 3; i++)
            atom_ptr->pos[i] += molecule_move->com[i];

    /* apply the rotation to the second molecule */
    if (system->surf_global_axis_on) {
        molecule_rotate_quaternion(molecule_origin, alpha_origin, beta_origin, gamma_origin, reverse_flag);
        molecule_rotate_quaternion(molecule_move, alpha_move, beta_move, gamma_move, reverse_flag);
    } else {
        molecule_rotate_euler(molecule_origin, alpha_origin, beta_origin, gamma_origin, reverse_flag);
        molecule_rotate_euler(molecule_move, alpha_move, beta_move, gamma_move, reverse_flag);
    }

    return (0);
}

void reset_molecule_position(molecule_t *mptr) {
    atom_t *aptr;
    int i;

    for (aptr = mptr->atoms; aptr; aptr = aptr->next)
        for (i = 0; i < 3; i++)
            aptr->pos[i] -= mptr->com[i];
    for (i = 0; i < 3; i++)
        mptr->com[i] = 0;

    return;
}

/* same idea as surface_dimer_geometry but here we need to keep the molecule at the origin fixed,
 * while moving the second molecule (x,y,z) and rotating it (a,b,c) in order to 
 * compute the virial integral with correct haar measure sin(b) dxdydzdadbdc */

/* set the distance and angles for a particular dimer orientation */
int surface_dimer_geometry_virial(system_t *system, double x, double y, double z, double a, double b, double c, int reverse_flag) {
    int i;
    molecule_t *molecule_origin, *molecule_move;
    atom_t *atom_ptr;

    // make sure that there are only two molecules
    molecule_move = system->molecules->next;
    if (!molecule_move) {
        error(
            "SURFACE: the input PQR has only a single molecule\n");
        return (-1);
    }
    if (molecule_move->next) {
        error(
            "SURFACE: the input PQR must contain exactly two molecules\n");
        return (-1);
    }

    // relocate movable molecule
    for (atom_ptr = molecule_move->atoms; atom_ptr; atom_ptr = atom_ptr->next)
        for (i = 0; i < 3; i++)
            atom_ptr->pos[i] -= molecule_move->com[i];
    for (i = 0; i < 3; i++)
        molecule_move->com[i] = 0;

    // apply the rotation to molecule #2
    if (system->surf_global_axis_on)
        molecule_rotate_quaternion(molecule_move, a, b, c, reverse_flag);
    else
        molecule_rotate_euler(molecule_move, a, b, c, reverse_flag);

    // relocate the moveable molecule
    molecule_move->com[0] = x;
    molecule_move->com[1] = y;
    molecule_move->com[2] = z;
    for (atom_ptr = molecule_move->atoms; atom_ptr; atom_ptr = atom_ptr->next)
        for (i = 0; i < 3; i++)
            atom_ptr->pos[i] += molecule_move->com[i];

    return (0);
}
