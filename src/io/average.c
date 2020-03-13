/*

Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

/* determine the acceptance rate */
void track_ar(nodestats_t *ns) {
    if (ns->accept + ns->reject)
        ns->acceptance_rate = ((double)ns->accept) / ((double)(ns->accept + ns->reject));
    else
        ns->acceptance_rate = 0;

    if (ns->accept_insert + ns->reject_insert)
        ns->acceptance_rate_insert = ((double)ns->accept_insert) / ((double)(ns->accept_insert + ns->reject_insert));
    else
        ns->acceptance_rate_insert = 0;

    if (ns->accept_remove + ns->reject_remove)
        ns->acceptance_rate_remove = ((double)ns->accept_remove) / ((double)(ns->accept_remove + ns->reject_remove));
    else
        ns->acceptance_rate_remove = 0;

    if (ns->accept_displace + ns->reject_displace)
        ns->acceptance_rate_displace = ((double)ns->accept_displace) / ((double)(ns->accept_displace + ns->reject_displace));
    else
        ns->acceptance_rate_displace = 0;

    if (ns->accept_adiabatic + ns->reject_adiabatic)
        ns->acceptance_rate_adiabatic = ((double)ns->accept_adiabatic) / ((double)(ns->accept_adiabatic + ns->reject_adiabatic));
    else
        ns->acceptance_rate_adiabatic = 0;

    if (ns->accept_spinflip + ns->reject_spinflip)
        ns->acceptance_rate_spinflip = ((double)ns->accept_spinflip) / ((double)(ns->accept_spinflip + ns->reject_spinflip));
    else
        ns->acceptance_rate_spinflip = 0;

    if (ns->accept_volume + ns->reject_volume)
        ns->acceptance_rate_volume = ((double)ns->accept_volume) / ((double)(ns->accept_volume + ns->reject_volume));
    else
        ns->acceptance_rate_volume = 0;

    if (ns->accept_ptemp + ns->reject_ptemp)
        ns->acceptance_rate_ptemp = ((double)ns->accept_ptemp) / ((double)(ns->accept_ptemp + ns->reject_ptemp));
    else
        ns->acceptance_rate_ptemp = 0;
}

/* update node statistics related to the processing */
void update_nodestats(nodestats_t *nodestats, avg_nodestats_t *avg_nodestats) {
    static int counter = 0;
    double factor;

    counter++;
    factor = ((double)(counter - 1)) / ((double)(counter));

    avg_nodestats->boltzmann_factor = factor * avg_nodestats->boltzmann_factor + nodestats->boltzmann_factor / ((double)counter);
    avg_nodestats->boltzmann_factor_sq = factor * avg_nodestats->boltzmann_factor_sq + nodestats->boltzmann_factor * nodestats->boltzmann_factor / ((double)counter);

    /* these really aren't averages, but accumulative values */
    avg_nodestats->acceptance_rate = nodestats->acceptance_rate;
    avg_nodestats->acceptance_rate_insert = nodestats->acceptance_rate_insert;
    avg_nodestats->acceptance_rate_remove = nodestats->acceptance_rate_remove;
    avg_nodestats->acceptance_rate_displace = nodestats->acceptance_rate_displace;
    avg_nodestats->acceptance_rate_adiabatic = nodestats->acceptance_rate_adiabatic;
    avg_nodestats->acceptance_rate_spinflip = nodestats->acceptance_rate_spinflip;
    avg_nodestats->acceptance_rate_volume = nodestats->acceptance_rate_volume;
    avg_nodestats->acceptance_rate_ptemp = nodestats->acceptance_rate_ptemp;

    avg_nodestats->cavity_bias_probability = factor * avg_nodestats->cavity_bias_probability + nodestats->cavity_bias_probability / ((double)counter);
    avg_nodestats->cavity_bias_probability_sq = factor * avg_nodestats->cavity_bias_probability_sq + nodestats->cavity_bias_probability * nodestats->cavity_bias_probability / ((double)counter);

    avg_nodestats->polarization_iterations = factor * avg_nodestats->polarization_iterations + nodestats->polarization_iterations / ((double)counter);
    avg_nodestats->polarization_iterations_sq = factor * avg_nodestats->polarization_iterations_sq + nodestats->polarization_iterations * nodestats->polarization_iterations / ((double)counter);
}

void count_sorbates(system_t *system) {
    int i;
    molecule_t *molecule_ptr;

    // Zero every count (N) for each sorbate in the averages list
    for (i = 0; i < system->sorbateCount; i++)
        system->sorbateInfo[i].currN = 0;

    // Count each sorbate in the system and record the total in
    // the corresponding entry in the sorbate averages list.
    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (i = 0; i < system->sorbateCount; i++) {
            if (!strcasecmp(system->sorbateInfo[i].id, molecule_ptr->moleculetype)) {
                system->sorbateInfo[i].currN++;
                break;
            }
        }
    }

    return;
}

// update the stats for the individual sorbates
void update_sorbate_info(system_t *system) {
    // molecule_t *molecule_ptr;  (unused variable)

    double sorbed_mass, pressure;
    int i;

    //update the number of particles of each sorbate
    count_sorbates(system);

    for (i = 0; i < system->sorbateCount; i++) {
        if (system->h2_fugacity || system->co2_fugacity || system->ch4_fugacity || system->n2_fugacity)
            pressure = system->fugacities[0];
        else if (system->user_fugacities)
            pressure = system->fugacities[i];
        else
            pressure = system->pressure;

        sorbed_mass = system->sorbateInfo[i].currN * system->sorbateInfo[i].mass;

        system->sorbateInfo[i].percent_wt = 100.0 * sorbed_mass / (system->observables->total_mass);
        system->sorbateInfo[i].percent_wt_me = 100.0 * sorbed_mass / (system->observables->frozen_mass);
        system->sorbateInfo[i].excess_ratio = 1000.0 * system->sorbateInfo[i].mass * (system->sorbateInfo[i].currN - system->sorbateInfo[i].mass * system->free_volume * pressure * ATM2REDUCED / system->temperature) /
                                              system->observables->frozen_mass;
        system->sorbateInfo[i].density = sorbed_mass / (system->observables->volume * NA * A32CM3);
        system->sorbateInfo[i].pore_density = sorbed_mass / (system->free_volume * NA * A32CM3);
    }

    return;
}

// calculate and set observable variabels for the frozen and total mass of the system
void calc_system_mass(system_t *system) {
    molecule_t *molecule_ptr;
    // double frozen_mass;  (unused variable)

    system->observables->total_mass = 0;
    system->observables->frozen_mass = 0;

    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        system->observables->total_mass += molecule_ptr->mass;
        if (molecule_ptr->frozen || molecule_ptr->adiabatic)
            system->observables->frozen_mass += molecule_ptr->mass;
    }

    return;
}

void update_root_sorb_averages(system_t *system, sorbateInfo_t *sinfo) {
    //sorbateGlobal is an array. sorbateStats is a linked list.

    sorbateAverages_t *sorbateGlobal = system->sorbateGlobal;

    static int counter = 0;
    double m, factor, sdom;
    double numerator, denominator, relative_err;
    int i, j;

    ++counter;
    m = (double)counter;
    sdom = 2.0 / sqrt(m - 1.0);
    factor = (m - 1.0) / m;

    // for each sorbate
    for (i = 0; i < system->sorbateCount; i++) {
        sorbateGlobal[i].avgN = factor * sorbateGlobal[i].avgN + sinfo[i].currN / m;
        sorbateGlobal[i].avgN_sq = factor * sorbateGlobal[i].avgN_sq + sinfo[i].currN * sinfo[i].currN / m;
        sorbateGlobal[i].avgN_err = sdom * sqrt(sorbateGlobal[i].avgN_sq - sorbateGlobal[i].avgN * sorbateGlobal[i].avgN);

        sorbateGlobal[i].percent_wt = factor * sorbateGlobal[i].percent_wt + sinfo[i].percent_wt / m;
        sorbateGlobal[i].percent_wt_sq = factor * sorbateGlobal[i].percent_wt_sq + sinfo[i].percent_wt * sinfo[i].percent_wt / m;
        sorbateGlobal[i].percent_wt_err = sdom * sqrt(sorbateGlobal[i].percent_wt_sq - sorbateGlobal[i].percent_wt * sorbateGlobal[i].percent_wt);

        sorbateGlobal[i].percent_wt_me = factor * sorbateGlobal[i].percent_wt_me + sinfo[i].percent_wt_me / m;
        sorbateGlobal[i].percent_wt_me_sq = factor * sorbateGlobal[i].percent_wt_me_sq + sinfo[i].percent_wt_me * sinfo[i].percent_wt_me / m;
        sorbateGlobal[i].percent_wt_me_err = sdom * sqrt(sorbateGlobal[i].percent_wt_me_sq - sorbateGlobal[i].percent_wt_me * sorbateGlobal[i].percent_wt_me);

        sorbateGlobal[i].excess_ratio = factor * sorbateGlobal[i].excess_ratio + sinfo[i].excess_ratio / m;
        sorbateGlobal[i].excess_ratio_sq = factor * sorbateGlobal[i].excess_ratio_sq + sinfo[i].excess_ratio * sinfo[i].excess_ratio / m;
        sorbateGlobal[i].excess_ratio_err = sdom * sqrt(sorbateGlobal[i].excess_ratio_sq - sorbateGlobal[i].excess_ratio * sorbateGlobal[i].excess_ratio);

        sorbateGlobal[i].pore_density = factor * sorbateGlobal[i].pore_density + sinfo[i].pore_density / m;
        sorbateGlobal[i].pore_density_sq = factor * sorbateGlobal[i].pore_density_sq + sinfo[i].pore_density * sinfo[i].pore_density / m;
        sorbateGlobal[i].pore_density_err = sdom * sqrt(sorbateGlobal[i].pore_density_sq - sorbateGlobal[i].pore_density * sorbateGlobal[i].pore_density);

        sorbateGlobal[i].density = factor * sorbateGlobal[i].density + sinfo[i].density / m;
        sorbateGlobal[i].density_sq = factor * sorbateGlobal[i].density_sq + sinfo[i].density * sinfo[i].density / m;
        sorbateGlobal[i].density_err = sdom * sqrt(sorbateGlobal[i].density_sq - sorbateGlobal[i].density * sorbateGlobal[i].density);
    }

    // calculate selectivity
    for (i = 0; i < system->sorbateCount; i++) {
        numerator = sorbateGlobal[i].avgN;
        relative_err = sorbateGlobal[i].avgN_err * sorbateGlobal[i].avgN_err / (sorbateGlobal[i].avgN * sorbateGlobal[i].avgN);
        denominator = 0;
        for (j = 0; j < system->sorbateCount; j++) {
            if (j == i) continue;
            denominator += sorbateGlobal[j].avgN;
            relative_err += sorbateGlobal[j].avgN_err * sorbateGlobal[j].avgN_err / (sorbateGlobal[j].avgN * sorbateGlobal[j].avgN);
        }
        sorbateGlobal[i].selectivity = numerator / denominator;
        sorbateGlobal[i].selectivity_err =
            sorbateGlobal[i].selectivity * sqrt(relative_err);
    }

    return;
}

void update_root_averages(system_t *system, observables_t *observables, avg_observables_t *avg_observables) {
    double particle_mass = 0, curr_density;
    double frozen_mass = system->observables->frozen_mass;

    molecule_t *molecule_ptr;
    static int counter = 0;
    double m, factor, gammaratio, sdom;

    ++counter;
    m = (double)counter;
    sdom = 2.0 / sqrt(m - 1.0);
    factor = (m - 1.0) / m;

    /* the physical observables */
    avg_observables->energy = factor * avg_observables->energy + observables->energy / m;
    avg_observables->energy_sq = factor * avg_observables->energy_sq + (observables->energy * observables->energy) / m;
    avg_observables->energy_error = sdom * sqrt(avg_observables->energy_sq - avg_observables->energy * avg_observables->energy);

    avg_observables->energy_sq_sq = factor * avg_observables->energy_sq_sq + pow(observables->energy, 4) / m;
    avg_observables->energy_sq_error = sdom * sqrt(avg_observables->energy_sq_sq - pow(avg_observables->energy, 4));

    avg_observables->coulombic_energy = factor * avg_observables->coulombic_energy + observables->coulombic_energy / m;
    avg_observables->coulombic_energy_sq = factor * avg_observables->coulombic_energy_sq + (observables->coulombic_energy * observables->coulombic_energy) / m;
    avg_observables->coulombic_energy_error = sdom * sqrt(avg_observables->coulombic_energy_sq - avg_observables->coulombic_energy * avg_observables->coulombic_energy);

    avg_observables->rd_energy = factor * avg_observables->rd_energy + observables->rd_energy / m;
    avg_observables->rd_energy_sq = factor * avg_observables->rd_energy_sq + (observables->rd_energy * observables->rd_energy) / m;
    avg_observables->rd_energy_error = sdom * sqrt(avg_observables->rd_energy_sq - avg_observables->rd_energy * avg_observables->rd_energy);

    avg_observables->polarization_energy = factor * avg_observables->polarization_energy + observables->polarization_energy / m;
    avg_observables->polarization_energy_sq = factor * avg_observables->polarization_energy_sq + (observables->polarization_energy * observables->polarization_energy) / m;
    avg_observables->polarization_energy_error = sdom * sqrt(avg_observables->polarization_energy_sq - avg_observables->polarization_energy * avg_observables->polarization_energy);

    avg_observables->vdw_energy = factor * avg_observables->vdw_energy + observables->vdw_energy / m;
    avg_observables->vdw_energy_sq = factor * avg_observables->vdw_energy_sq + (observables->vdw_energy * observables->vdw_energy) / m;
    avg_observables->vdw_energy_error = sdom * sqrt(avg_observables->vdw_energy_sq - avg_observables->vdw_energy * avg_observables->vdw_energy);

    avg_observables->three_body_energy = factor * avg_observables->three_body_energy + observables->three_body_energy / m;
    avg_observables->three_body_energy_sq = factor * avg_observables->three_body_energy_sq + (observables->three_body_energy * observables->three_body_energy) / m;
    avg_observables->three_body_energy_error = sdom * sqrt(avg_observables->three_body_energy_sq - avg_observables->three_body_energy * avg_observables->three_body_energy);

    avg_observables->dipole_rrms = factor * avg_observables->dipole_rrms + observables->dipole_rrms / m;
    avg_observables->dipole_rrms_sq = factor * avg_observables->dipole_rrms_sq + (observables->dipole_rrms * observables->dipole_rrms) / m;
    avg_observables->dipole_rrms_error = sdom * sqrt(avg_observables->dipole_rrms_sq - avg_observables->dipole_rrms * avg_observables->dipole_rrms);

    avg_observables->kinetic_energy = factor * avg_observables->kinetic_energy + observables->kinetic_energy / m;
    avg_observables->kinetic_energy_sq = factor * avg_observables->kinetic_energy_sq + (observables->kinetic_energy * observables->kinetic_energy) / m;
    avg_observables->kinetic_energy_error = sdom * sqrt(avg_observables->kinetic_energy_sq - avg_observables->kinetic_energy * avg_observables->kinetic_energy);

    avg_observables->temperature = factor * avg_observables->temperature + observables->temperature / m;
    avg_observables->temperature_sq = factor * avg_observables->temperature_sq + (observables->temperature * observables->temperature) / m;
    avg_observables->temperature_error = sdom * sqrt(avg_observables->temperature_sq - avg_observables->temperature * avg_observables->temperature);

    avg_observables->volume = factor * avg_observables->volume + observables->volume / m;
    avg_observables->volume_sq = factor * avg_observables->volume_sq + (observables->volume * observables->volume) / m;
    avg_observables->volume_error = sdom * sqrt(avg_observables->volume_sq - avg_observables->volume * avg_observables->volume);

    avg_observables->N = factor * avg_observables->N + observables->N / m;
    avg_observables->N_sq = factor * avg_observables->N_sq + (observables->N * observables->N) / m;
    avg_observables->N_error = sdom * sqrt(avg_observables->N_sq - avg_observables->N * avg_observables->N);

    avg_observables->spin_ratio = factor * avg_observables->spin_ratio + observables->spin_ratio / m;
    avg_observables->spin_ratio_sq = factor * avg_observables->spin_ratio_sq + (observables->spin_ratio * observables->spin_ratio) / m;
    avg_observables->spin_ratio_error = sdom * sqrt(avg_observables->spin_ratio_sq - avg_observables->spin_ratio * avg_observables->spin_ratio);

    avg_observables->NU = factor * avg_observables->NU + observables->NU / m;

    /* particle mass will be used in calculations for single sorbate systems */
    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
        if (!molecule_ptr->frozen && !molecule_ptr->adiabatic)
            particle_mass = molecule_ptr->mass;

    /* density in g/cm^3 */  //had to modify since density isn't neccessarily constant since adding NPT
    curr_density = observables->N * particle_mass / (system->pbc->volume * NA * A32CM3);

    avg_observables->density = factor * avg_observables->density + curr_density / m;
    avg_observables->density_sq = factor * avg_observables->density_sq + (curr_density * curr_density) / m;
    avg_observables->density_error = sdom * sqrt(avg_observables->density_sq - (avg_observables->density) * (avg_observables->density));

    /* needed for calculating sstdev (stdev of stdev) */
    /*
	gammaratio = tgamma(0.5*(double)counter) / tgamma(0.5*((double)counter-1));
	gammaratio = sqrt(1.0/counter * ((double)counter - 1.0 - 2.0*gammaratio*gammaratio) );
	*/
    //stirling approx to avoid numerical overflow
    gammaratio = pow((m - 2.0) / (m - 1.0), 0.5 * m - 1.0) * sqrt(0.5 * (m - 2.0)) * exp(0.5);
    gammaratio = sqrt(1.0 / counter * (m - 1.0 - 2.0 * gammaratio * gammaratio));

    /* heat capacity in kJ/mol K */
    avg_observables->heat_capacity = (KB * NA / 1000.0) * (avg_observables->energy_sq - avg_observables->energy * avg_observables->energy) / (system->temperature * system->temperature);
    /* error in heat capacity is the standard deviation of the variance, = 2*sstdev */
    avg_observables->heat_capacity_error = sdom * 2.0 * gammaratio * avg_observables->heat_capacity;

    /* compressibility */
    if (system->ensemble != ENSEMBLE_NPT)
        avg_observables->compressibility = ATM2PASCALS * (system->pbc->volume / pow(METER2ANGSTROM, 3)) * (avg_observables->N_sq - avg_observables->N * avg_observables->N) / (KB * system->temperature * avg_observables->N * avg_observables->N);
    else
        avg_observables->compressibility = ATM2PASCALS * pow(METER2ANGSTROM, -3) *
                                           (avg_observables->volume_sq - avg_observables->volume * avg_observables->volume) /
                                           (KB * system->temperature * avg_observables->volume);
    avg_observables->compressibility_error = sdom * 2.0 * gammaratio * avg_observables->compressibility;

    /* we have a solid phase */
    if (frozen_mass > 0.0) {
        /* percent weight */
        avg_observables->percent_wt = 100.0 * avg_observables->N * particle_mass / (frozen_mass + avg_observables->N * particle_mass);
        avg_observables->percent_wt_error = sdom * 100.0 * avg_observables->N_error * particle_mass / (frozen_mass + avg_observables->N_error * particle_mass);

        /* percent weight like ME*/
        avg_observables->percent_wt_me = 100.0 * avg_observables->N * particle_mass / frozen_mass;
        avg_observables->percent_wt_me_error = sdom * 100.0 * avg_observables->N_error * particle_mass / frozen_mass;

        /* excess weight mg/g */
        if (system->free_volume > 0.0) {
            if (system->fugacities)
                avg_observables->excess_ratio = 1000.0 * (avg_observables->N * particle_mass - (particle_mass * system->free_volume * system->fugacities[0] * ATM2REDUCED) / system->temperature) / frozen_mass;
            else
                avg_observables->excess_ratio = 1000.0 * (avg_observables->N * particle_mass - (particle_mass * system->free_volume * system->pressure * ATM2REDUCED) / system->temperature) / frozen_mass;
            avg_observables->excess_ratio_error = sdom * 1000.0 * avg_observables->N_error * particle_mass / frozen_mass;

            /* pore density */  //only valid for constant V, pure systems
            avg_observables->pore_density = curr_density * system->pbc->volume / system->free_volume;
            avg_observables->pore_density_error = sdom * avg_observables->N_error * particle_mass / (system->free_volume * NA * A32CM3);
        }

        /* calculate the isosteric heat */
        avg_observables->qst = -(avg_observables->NU - avg_observables->N * avg_observables->energy);
        avg_observables->qst /= (avg_observables->N_sq - avg_observables->N * avg_observables->N);
        avg_observables->qst += system->temperature;
        avg_observables->qst *= KB * NA / 1000.0; /* convert to kJ/mol */

        /* and the NVT-style Qst */
        avg_observables->qst_nvt = -avg_observables->energy * KB * NA / 1000.0 / avg_observables->N;
    }

    return;
}

void clear_avg_nodestats(system_t *system) {
    avg_observables_t *avg_observables = system->avg_observables;
    avg_nodestats_t *avg_nodestats = system->avg_nodestats;

    avg_nodestats->counter = 0;

    avg_observables->boltzmann_factor = 0;
    avg_observables->boltzmann_factor_sq = 0;

    avg_observables->acceptance_rate = 0;
    avg_observables->acceptance_rate_insert = 0;
    avg_observables->acceptance_rate_remove = 0;
    avg_observables->acceptance_rate_displace = 0;
    avg_observables->acceptance_rate_adiabatic = 0;
    avg_observables->acceptance_rate_spinflip = 0;
    avg_observables->acceptance_rate_volume = 0;
    avg_observables->acceptance_rate_ptemp = 0;

    avg_observables->cavity_bias_probability = 0;
    avg_observables->cavity_bias_probability_sq = 0;

    avg_observables->polarization_iterations = 0;
    avg_observables->polarization_iterations_sq = 0;

    return;
}

void update_root_nodestats(system_t *system, avg_nodestats_t *avg_nodestats, avg_observables_t *avg_observables) {
    double m, factor, sdom;

    m = (double)(++system->avg_nodestats->counter);
    sdom = 2.0 / sqrt(floor((double)((system->step + 1.0) * size) / (double)(system->corrtime)) - 1.0);
    factor = (m - 1.0) / m;

    avg_observables->boltzmann_factor = factor * avg_observables->boltzmann_factor + avg_nodestats->boltzmann_factor / m;
    avg_observables->boltzmann_factor_sq = factor * avg_observables->boltzmann_factor_sq + avg_nodestats->boltzmann_factor_sq / m;
    avg_observables->boltzmann_factor_error = sdom * sqrt(avg_observables->boltzmann_factor_sq - avg_observables->boltzmann_factor * avg_observables->boltzmann_factor);

    avg_observables->acceptance_rate = factor * avg_observables->acceptance_rate + avg_nodestats->acceptance_rate / m;
    avg_observables->acceptance_rate_insert = factor * avg_observables->acceptance_rate_insert + avg_nodestats->acceptance_rate_insert / m;
    avg_observables->acceptance_rate_remove = factor * avg_observables->acceptance_rate_remove + avg_nodestats->acceptance_rate_remove / m;
    avg_observables->acceptance_rate_displace = factor * avg_observables->acceptance_rate_displace + avg_nodestats->acceptance_rate_displace / m;
    avg_observables->acceptance_rate_adiabatic = factor * avg_observables->acceptance_rate_adiabatic + avg_nodestats->acceptance_rate_adiabatic / m;
    avg_observables->acceptance_rate_spinflip = factor * avg_observables->acceptance_rate_spinflip + avg_nodestats->acceptance_rate_spinflip / m;
    avg_observables->acceptance_rate_volume = factor * avg_observables->acceptance_rate_volume + avg_nodestats->acceptance_rate_volume / m;
    avg_observables->acceptance_rate_ptemp = factor * avg_observables->acceptance_rate_ptemp + avg_nodestats->acceptance_rate_ptemp / m;

    avg_observables->cavity_bias_probability = factor * avg_observables->cavity_bias_probability + avg_nodestats->cavity_bias_probability / m;
    avg_observables->cavity_bias_probability_sq = factor * avg_observables->cavity_bias_probability_sq + avg_nodestats->cavity_bias_probability_sq / m;
    avg_observables->cavity_bias_probability_error = sdom * sqrt(avg_observables->cavity_bias_probability_sq - avg_observables->cavity_bias_probability * avg_observables->cavity_bias_probability);

    avg_observables->polarization_iterations = factor * avg_observables->polarization_iterations + avg_nodestats->polarization_iterations / m;
    avg_observables->polarization_iterations_sq = factor * avg_observables->polarization_iterations_sq + avg_nodestats->polarization_iterations_sq / m;
    avg_observables->polarization_iterations_error = sdom * sqrt(avg_observables->polarization_iterations_sq - avg_observables->polarization_iterations * avg_observables->polarization_iterations);

    return;
}
