#include <mc.h>

#define COLOR_H \
    "0.2 0.2 0.2"
#define COLOR_C \
    "0.1 0.5 0.1"
#define COLOR_N \
    "0.2 0.2 1.0"
#define COLOR_O \
    "1.0 0.0 0.0"
#define COLOR_XXX \
    "0.1 0.1 0.1"

void print_frozen_coords(FILE *fp_frozen, system_t *system) {
    molecule_t *mol;
    atom_t *atom;

    for (mol = system->molecules; mol; mol = mol->next) {
        if (mol->frozen) {
            for (atom = mol->atoms; atom; atom = atom->next) {
                fprintf(fp_frozen,
                        "%f %f %f\n", atom->pos[0], atom->pos[1], atom->pos[2]);
            }
        }
    }
}

int count_frozen(system_t *system) {
    int count = 0;
    molecule_t *mol;
    atom_t *atom;

    for (mol = system->molecules; mol; mol = mol->next) {
        if (mol->frozen) {
            for (atom = mol->atoms; atom; atom = atom->next) count++;
        }
    }
    return count;
}

int bondlength_check(atom_t *atom1, atom_t *atom2, system_t *system) {
    double distance;
    double gm_mass; /* geometric mean of mass */
    int is_bonded;
    double slope = 0.0234; /* tune to meet bond length expectations */
    double yint = 0.603;   /* tune to meet bond length expectations */

    gm_mass = sqrt(atom1->mass * atom2->mass);
    distance = sqrt(pow(atom1->pos[0] - atom2->pos[0], 2) + pow(atom1->pos[1] - atom2->pos[1], 2) + pow(atom1->pos[2] - atom2->pos[2], 2));

    if (distance < (gm_mass * slope + yint) * system->max_bondlength)
        is_bonded = 1;
    else
        is_bonded = 0;
    return is_bonded;
}

int calculate_bonds(system_t *system) {
    int inner_index = 0, outer_index = 0;
    int bonds = 0;
    // double bondlength; (unused variable)
    molecule_t *mol;
    atom_t *atom;
    atom_t *atom2;

    for (mol = system->molecules; mol; mol = mol->next) {
        if (mol->frozen) {
            for (atom = mol->atoms; atom; atom = atom->next, inner_index++) {
                if (atom->next) {
                    for (atom2 = atom->next, outer_index = inner_index + 1; atom2; atom2 = atom2->next, outer_index++) {
                        if (bondlength_check(atom, atom2, system)) {
                            bonds++;
                        }
                    }
                }
            }
        }
    }
    return bonds;
}

void print_frozen_bonds(FILE *fp_frozen, system_t *system) {
    int inner_index = 0, outer_index = 0;
    int bonds = 0;
    // double bondlength;  (unused variable)
    molecule_t *mol;
    atom_t *atom;
    atom_t *atom2;

    for (mol = system->molecules; mol; mol = mol->next) {
        if (mol->frozen) {
            for (atom = mol->atoms; atom; atom = atom->next, inner_index++) {
                if (atom->next) {
                    for (atom2 = atom->next, outer_index = inner_index + 1; atom2; atom2 = atom2->next, outer_index++) {
                        if (bondlength_check(atom, atom2, system)) {
                            bonds++;
                            fprintf(fp_frozen,
                                    "%d %d\n", inner_index, outer_index);
                        }
                    }
                }
            }
        }
    }
}

void print_frozen_masses(FILE *fp_frozen, system_t *system) {
    molecule_t *mol;
    atom_t *atom;

    for (mol = system->molecules; mol; mol = mol->next) {
        if (mol->frozen) {
            for (atom = mol->atoms; atom; atom = atom->next) {
                fprintf(fp_frozen,
                        "%f\n", atom->mass);
            }
        }
    }
}

void print_frozen_colors(FILE *fp_frozen, system_t *system) {
    double mass;
    molecule_t *mol;
    atom_t *atom;

    for (mol = system->molecules; mol; mol = mol->next) {
        if (mol->frozen) {
            for (atom = mol->atoms; atom; atom = atom->next) {
                mass = atom->mass;
                if (mass < 1.1)
                    fprintf(fp_frozen,
                            "%s\n", COLOR_H);
                else if (mass < 12.2)
                    fprintf(fp_frozen,
                            "%s\n", COLOR_C);
                else if (mass < 14.1)
                    fprintf(fp_frozen,
                            "%s\n", COLOR_N);
                else if (mass < 16.1)
                    fprintf(fp_frozen,
                            "%s\n", COLOR_O);
                else
                    fprintf(fp_frozen,
                            "%s\n", COLOR_XXX);
            }
        }
    }
}

void write_frozen(FILE *fp_frozen, system_t *system) {
    int numatoms;
    int numbonds;

    numatoms = count_frozen(system);
    numbonds = calculate_bonds(system);

    rewind(fp_frozen);
    fprintf(fp_frozen,
            "# OpenDX format coordinate file for frozen atoms\n");
    fprintf(fp_frozen,
            "object 1 class array type float rank 1 shape 3 items %d data follows\n", numatoms);
    print_frozen_coords(fp_frozen, system);
    fprintf(fp_frozen,
            "object 2 class array type int rank 1 shape 2 items %d data follows\n", numbonds);
    print_frozen_bonds(fp_frozen, system);
    fprintf(fp_frozen,
            "attribute \"element type\" string \"lines\"\n");
    fprintf(fp_frozen,
            "attribute \"ref\" string \"positions\"\n");
    fprintf(fp_frozen,
            "object 3 class array type float rank 0 items %d data follows\n", numatoms);
    print_frozen_masses(fp_frozen, system);
    fprintf(fp_frozen,
            "attribute \"dep\" string \"positions\"\n");
    fprintf(fp_frozen,
            "object 4 class array type float rank 1 shape 3 items %d data follows\n", numatoms);
    print_frozen_colors(fp_frozen, system);
    fprintf(fp_frozen,
            "object \"irregular positions irregular connections\" class field\n");
    fprintf(fp_frozen,
            "component \"positions\" value 1\n");
    fprintf(fp_frozen,
            "component \"connections\" value 2\n");
    fprintf(fp_frozen,
            "component \"data\" value 3\n");
    fprintf(fp_frozen,
            "component \"colors\" value 4\n");
    fprintf(fp_frozen,
            "end\n");
}
