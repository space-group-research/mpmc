#include <mc.h>

int calculate_te(system_t *system) {
    energy(system);

#ifdef PRINT_EFIELD
    printf(
        "OUT: MUARRAY %lf %lf\n", system->molecules->atoms->ef_static[0], system->molecules->atoms->mu[0]);
    printf(
        "OUT: MUARRAY %lf %lf\n", system->molecules->next->atoms->ef_static[0], system->molecules->next->atoms->mu[0]);
#endif

#ifdef QM_ROTATION
    // solve for the rotational energy levels
    if (system->quantum_rotation) quantum_system_rotational_energies(system);
#endif  // QM_ROTATION

    // if root, open necessary output files
    if (!rank) {
        print_observables(system);
        if (open_files(system) < 0) {
            error(
                "CALC_TE: could not open files\n");
            return (-1);
        }
    }

    if (system->file_pointers.fp_energy)
        write_observables(system->file_pointers.fp_energy, system, system->observables, system->temperature);
    if (system->file_pointers.fp_energy_csv)
        write_observables_csv(system->file_pointers.fp_energy_csv, system, system->observables, system->temperature);

    // close any open files
    if (!rank)
        close_files(system);

    return (0);
}
