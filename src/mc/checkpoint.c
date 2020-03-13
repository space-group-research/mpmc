#include <mc.h>

/* this function (a) determines what move will be made next time make_move() is called and (b) backups up the state to be restore upon rejection */
void checkpoint(system_t *system) {
    int j;
    int i_exchange, i_adiabatic;
    int num_molecules_exchange, num_molecules_adiabatic, altered;
    double num_molecules_adiabatic_double;
    molecule_t **ptr_array_exchange, **ptr_array_adiabatic;
    molecule_t *molecule_ptr, *prev_molecule_ptr;

    int alt;
    // sorbateInfo_t * sorbate_ptr;  (unused variable)

    /* save the current observables */
    memcpy(system->checkpoint->observables, system->observables, sizeof(observables_t));

    /* count exchangeable and adiabatic molecules */
    num_molecules_exchange = 0;
    num_molecules_adiabatic = 0;
    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        if (!(molecule_ptr->frozen || molecule_ptr->adiabatic || molecule_ptr->target))
            ++num_molecules_exchange;
        if (molecule_ptr->adiabatic)
            ++num_molecules_adiabatic;
    }

    /* go thru again, make an array of all eligible molecules */
    ptr_array_exchange = calloc(num_molecules_exchange, sizeof(molecule_t *));
    memnullcheck(ptr_array_exchange, num_molecules_exchange * sizeof(molecule_t *), __LINE__ - 1, __FILE__);
    ptr_array_adiabatic = calloc(num_molecules_adiabatic, sizeof(molecule_t *));
    memnullcheck(ptr_array_adiabatic, num_molecules_adiabatic * sizeof(molecule_t *), __LINE__ - 1, __FILE__);
    for (molecule_ptr = system->molecules, i_exchange = 0, i_adiabatic = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        if (!(molecule_ptr->frozen || molecule_ptr->adiabatic || molecule_ptr->target)) {
            ptr_array_exchange[i_exchange] = molecule_ptr;
            ++i_exchange;
        }

        if (molecule_ptr->adiabatic) {
            ptr_array_adiabatic[i_adiabatic] = molecule_ptr;
            ++i_adiabatic;
        }
    }

    /* determine what kind of move to do */
    switch (system->ensemble) {
        case ENSEMBLE_UVT:

            if (get_rand(system) < system->insert_probability) { /* do a particle insertion/deletion move */

                if (get_rand(system) < 0.5) { /* INSERT */
                    system->checkpoint->movetype = MOVETYPE_INSERT;
                } else { /* REMOVE */
                    system->checkpoint->movetype = MOVETYPE_REMOVE;
                }

            } else if (system->quantum_rotation) {
                if (get_rand(system) < system->spinflip_probability) /* SPINFLIP */
                    system->checkpoint->movetype = MOVETYPE_SPINFLIP;
                else {
                    if (num_molecules_adiabatic && (get_rand(system) < 0.5)) /* DISPLACE */
                        system->checkpoint->movetype = MOVETYPE_ADIABATIC;   /* for the adiabatic mole fraction */
                    else
                        system->checkpoint->movetype = MOVETYPE_DISPLACE;
                }

            } else {                                                     /* DISPLACE */
                if (num_molecules_adiabatic && (get_rand(system) < 0.5)) /* DISPLACE */
                    system->checkpoint->movetype = MOVETYPE_ADIABATIC;   /* for the adiabatic mole fraction */
                else
                    system->checkpoint->movetype = MOVETYPE_DISPLACE;
            }

            break;
        case ENSEMBLE_NVT:
        case ENSEMBLE_NVE:

            if (system->quantum_rotation && (get_rand(system) < system->spinflip_probability))
                system->checkpoint->movetype = MOVETYPE_SPINFLIP;
            else
                system->checkpoint->movetype = MOVETYPE_DISPLACE;

            break;
        case ENSEMBLE_NPT:

            if (system->volume_probability == 0.0) {  //if volume probability isn't set, then do volume moves with prob = 1/nmolecules
                if (get_rand(system) < 1.0 / system->observables->N)
                    system->checkpoint->movetype = MOVETYPE_VOLUME;
                else
                    system->checkpoint->movetype = MOVETYPE_DISPLACE;
            } else {  //if volume probability IS set
                if (get_rand(system) < system->volume_probability)
                    system->checkpoint->movetype = MOVETYPE_VOLUME;
                else
                    system->checkpoint->movetype = MOVETYPE_DISPLACE;
            }

            break;
        default:
            error(
                "CHECKPOINT: invalid ensemble\n");
            die(-1);
    }

    /* if we have any adiabatic molecules, then we have to do some special stuff */
    /* randomly pick a (moveable) atom */
    if (system->checkpoint->movetype == MOVETYPE_ADIABATIC) {
        --num_molecules_adiabatic;
        num_molecules_adiabatic_double = (double)num_molecules_adiabatic;
        altered = num_molecules_adiabatic - (int)rint(num_molecules_adiabatic_double * get_rand(system));
        system->checkpoint->molecule_altered = ptr_array_adiabatic[altered];

    } else {
        if (system->num_insertion_molecules && system->checkpoint->movetype == MOVETYPE_INSERT) {
            // If the move is an insertion and a molecule insertion list
            // has been set up, then select a molecule from said list:
            alt = (int)floor(get_rand(system) * (double)system->num_insertion_molecules);
            system->checkpoint->molecule_altered = system->insertion_molecules_array[alt];
            //needed to calculate boltzmann factor
            system->sorbateInsert = alt;

        }  // multi-sorbate + insert
        else {
            // Otherwise, perform the original MPMC treatment:
            --num_molecules_exchange;
            altered = (int)floor(get_rand(system) * system->observables->N);
            system->checkpoint->molecule_altered = ptr_array_exchange[altered];

            // if multi sorbate, we need to record the type of sorbate removed
            if (system->num_insertion_molecules && system->checkpoint->movetype == MOVETYPE_REMOVE) {
                alt = 0;
                for (j = 0; j < system->sorbateCount; j++) {
                    if (!strcasecmp(system->sorbateInfo[j].id, ptr_array_exchange[altered]->moleculetype)) {
                        system->sorbateInsert = alt;
                        break;
                    }
                    alt++;
                }
            }  // multi-sorbate + remove

        }  //end else multi-sorbate + insert
    }      //end else adiabatic

    /* free our temporary eligibility lists */
    free(ptr_array_exchange);
    if (ptr_array_adiabatic) free(ptr_array_adiabatic);

    /* never completely empty the list */
    if (!num_molecules_exchange && system->checkpoint->movetype == MOVETYPE_REMOVE) {
        if (system->quantum_rotation && (get_rand(system) < system->spinflip_probability))
            system->checkpoint->movetype = MOVETYPE_SPINFLIP;
        else
            system->checkpoint->movetype = MOVETYPE_DISPLACE;
    }

    // Determine the head and tail of the selected molecule
    if (system->num_insertion_molecules && system->checkpoint->movetype == MOVETYPE_INSERT) {
        // When using an insertion list, we will always insert at the end of the list.

        // Advance to the end of the list
        molecule_ptr = system->molecules;
        while (molecule_ptr->next)
            molecule_ptr = molecule_ptr->next;

        // Set head to the last molecule in the list, tail to the list terminator.
        system->checkpoint->head = molecule_ptr;
        system->checkpoint->tail = NULL;

    } else {
        // If this is not a insertion using an insertion list, perform the original MPMC treatment:
        for (molecule_ptr = system->molecules, prev_molecule_ptr = NULL; molecule_ptr; molecule_ptr = molecule_ptr->next) {
            if (molecule_ptr == system->checkpoint->molecule_altered) {
                system->checkpoint->head = prev_molecule_ptr;
                system->checkpoint->tail = molecule_ptr->next;
            }
            prev_molecule_ptr = molecule_ptr;
        }
    }

    /* if we have a molecule already backed up (from a previous accept), go ahead and free it */
    if (system->checkpoint->molecule_backup) free_molecule(system, system->checkpoint->molecule_backup);
    /* backup the state that will be altered */
    system->checkpoint->molecule_backup = copy_molecule(system, system->checkpoint->molecule_altered);

    return;
}
