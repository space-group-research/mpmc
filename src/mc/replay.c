#include <mc.h>

//return less than 0 if error
int replay_trajectory(system_t* system) {
    //need to read in molecules and shit and loop through the trajectory file.

    int errchk;
    double initial_energy, final_energy;
    char linebuf[MAXLINE];
    double p, psum, psum2;
    psum2 = psum = 0;
    FILE* finput = fopen(system->traj_input,
                         "r");
    filecheck(finput, system->traj_input, READ);
    rewind(finput);
    if (system->sorbateCount > 1) {
        system->sorbateGlobal = calloc(system->sorbateCount, sizeof(sorbateAverages_t));
        memnullcheck(system->sorbateGlobal, sizeof(sorbateAverages_t), __LINE__ - 1, __FILE__);
    }

    //now we need to loop. we will follow this procedure
    // 1. setup simulation box
    // 2. calculate shit
    // 3. clean up -> free molecules and atom

    system->step = 0;
    if (!rank)
        if (open_files(system) < 0) {
            error(
                "REPLAY: could not open files\n");
            return (-1);
        }

    while (1) {
        //remove old structures from previous loop or from input.c prior to calling this function
        free_all_pairs(system);
        free_all_molecules(system->molecules);
        if (system->insertion_molecules_array) {
            free(system->insertion_molecules_array);
            system->insertion_molecules_array = NULL;
        }
        if (system->insertion_molecules) {
            free_all_molecules(system->insertion_molecules);
            system->insertion_molecules = NULL;
        }

        errchk = setup_simulation_box(finput, system);
        if (errchk == 1)  //if 1, we are out of molecules
            break;
        if (errchk != 1 && errchk != 0) {  //some error
            output(
                "REPLAY: simulation box not properly set up\n");
            die(-1);
        } else
            system->step++;

        //set up pairs
        setup_pairs(system);
        if (system->spectre) spectre_wrapall(system);

        // set volume observable
        system->observables->volume = system->pbc->volume;

        // calculate the energy of the system
        initial_energy = energy(system);

#ifdef QM_ROTATION
        // solve for the rotational energy levels
        if (system->quantum_rotation) quantum_system_rotational_energies(system);
#endif  // QM_ROTATION

        calc_system_mass(system);
        update_nodestats(system->nodestats, system->avg_nodestats);
        update_root_nodestats(system, system->avg_nodestats, system->avg_observables);
        update_root_averages(system, system->observables, system->avg_observables);
        if (system->sorbateCount > 1) {
            update_sorbate_info(system);
            update_root_sorb_averages(system, system->sorbateInfo);
        }

        if (system->file_pointers.fp_energy)
            write_observables(system->file_pointers.fp_energy, system, system->observables, system->temperature);
        if (system->file_pointers.fp_energy_csv)
            write_observables_csv(system->file_pointers.fp_energy_csv, system, system->observables, system->temperature);

        // calculate pressure if requested
        if (system->calc_pressure) {
            volume_change(system);
            final_energy = energy(system);
            p = exp(-(final_energy - initial_energy) / system->temperature) *
                pow(system->observables->volume / (system->observables->volume - system->calc_pressure_dv), system->observables->N);
            psum += p;
            psum2 += p * p;
        }

        write_averages(system);

    }  //get next trajectory snapshot

    //calculate and output pressure
    if (system->calc_pressure) {
        psum /= (double)(system->step);  //average
        psum2 /= (double)(system->step);
        psum2 = system->temperature / system->calc_pressure_dv * KB * pow(METER2ANGSTROM, 3) / ATM2PASCALS *
                sqrt(psum2 - psum * psum) / psum;
        psum = system->temperature / system->calc_pressure_dv * KB * pow(METER2ANGSTROM, 3) / ATM2PASCALS *
               log(psum);
        sprintf(linebuf,
                "REPLAY: Pressure = %.5lf +- %.5lf atm\n", psum, psum2);
        output(linebuf);
    }

    // close any open files
    fclose(finput);
    if (!rank)
        close_files(system);

    free(system->sorbateGlobal);

    return (0);
}
