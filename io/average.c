/*

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

void clear_nodestats(nodestats_t *stats) {

	memset(stats, 0, sizeof(nodestats_t));

}

void clear_node_averages(avg_nodestats_t *avgstats) {

	memset(avgstats, 0, sizeof(avg_nodestats_t));

}

void clear_observables(observables_t *obs) {

	memset(obs, 0, sizeof(observables_t));

}

void clear_root_averages(avg_observables_t *avgobs) {

	memset(avgobs, 0, sizeof(avg_observables_t));

}

/* determine the acceptance rate */
void track_ar(nodestats_t *ns) {


	if(ns->accept + ns->reject)
		ns->acceptance_rate = ((double)ns->accept) / ((double)(ns->accept + ns->reject));
	else
		ns->acceptance_rate = 0;

	if(ns->accept_insert + ns->reject_insert)
		ns->acceptance_rate_insert = ((double)ns->accept_insert) / ((double)(ns->accept_insert + ns->reject_insert));
	else
		ns->acceptance_rate_insert = 0;

	if(ns->accept_remove + ns->reject_remove)
		ns->acceptance_rate_remove = ((double)ns->accept_remove) / ((double)(ns->accept_remove + ns->reject_remove));
	else
		ns->acceptance_rate_remove = 0;

	if(ns->accept_displace + ns->reject_displace)
		ns->acceptance_rate_displace = ((double)ns->accept_displace) / ((double)(ns->accept_displace + ns->reject_displace));
	else
		ns->acceptance_rate_displace = 0;

	if(ns->accept_adiabatic + ns->reject_adiabatic)
		ns->acceptance_rate_adiabatic = ((double)ns->accept_adiabatic) / ((double)(ns->accept_adiabatic + ns->reject_adiabatic));
	else
		ns->acceptance_rate_adiabatic = 0;

	if(ns->accept_spinflip + ns->reject_spinflip)
		ns->acceptance_rate_spinflip = ((double)ns->accept_spinflip) / ((double)(ns->accept_spinflip + ns->reject_spinflip));
	else
		ns->acceptance_rate_spinflip = 0;

}

/* update node statistics related to the processing */
void update_nodestats(nodestats_t *nodestats, avg_nodestats_t *avg_nodestats) {

	static int counter = 0;
	double factor;

	++counter;
	factor = ((double)(counter - 1))/((double)(counter));

	avg_nodestats->boltzmann_factor = factor*avg_nodestats->boltzmann_factor + nodestats->boltzmann_factor / ((double)counter);
	avg_nodestats->boltzmann_factor_sq = factor*avg_nodestats->boltzmann_factor_sq + nodestats->boltzmann_factor*nodestats->boltzmann_factor / ((double)counter);

	/* these really aren't averages, but accumulative values */
	avg_nodestats->acceptance_rate = nodestats->acceptance_rate;
	avg_nodestats->acceptance_rate_insert = nodestats->acceptance_rate_insert;
	avg_nodestats->acceptance_rate_remove = nodestats->acceptance_rate_remove;
	avg_nodestats->acceptance_rate_displace = nodestats->acceptance_rate_displace;
	avg_nodestats->acceptance_rate_adiabatic = nodestats->acceptance_rate_adiabatic;
	avg_nodestats->acceptance_rate_spinflip = nodestats->acceptance_rate_spinflip;

	avg_nodestats->cavity_bias_probability = factor*avg_nodestats->cavity_bias_probability + nodestats->cavity_bias_probability / ((double)counter);
	avg_nodestats->cavity_bias_probability_sq = factor*avg_nodestats->cavity_bias_probability_sq + nodestats->cavity_bias_probability*nodestats->cavity_bias_probability / ((double)counter);

	avg_nodestats->polarization_iterations = factor*avg_nodestats->polarization_iterations + nodestats->polarization_iterations / ((double)counter);
	avg_nodestats->polarization_iterations_sq = factor*avg_nodestats->polarization_iterations_sq + nodestats->polarization_iterations*nodestats->polarization_iterations / ((double)counter);


}



void update_root_averages(system_t *system, observables_t *observables, avg_nodestats_t *avg_nodestats, avg_observables_t *avg_observables) {

	double particle_mass, frozen_mass;

	molecule_t *molecule_ptr;
	static int counter = 0;
	double factor;

	++counter;
	factor = ((double)(counter - 1))/((double)(counter));


	/* the physical observables */
	avg_observables->energy = factor*avg_observables->energy + observables->energy / ((double)counter);
	avg_observables->energy_sq = factor*avg_observables->energy_sq + (observables->energy*observables->energy) / ((double)counter);
	avg_observables->energy_error = 0.5*sqrt(avg_observables->energy_sq  - avg_observables->energy*avg_observables->energy);

	avg_observables->energy_sq_sq = factor*avg_observables->energy_sq_sq + pow(observables->energy, 4.0) / ((double)counter);
	avg_observables->energy_sq_error = 0.5*sqrt(avg_observables->energy_sq_sq  - pow(avg_observables->energy, 4.0));

	avg_observables->coulombic_energy = factor*avg_observables->coulombic_energy + observables->coulombic_energy / ((double)counter);
	avg_observables->coulombic_energy_sq = factor*avg_observables->coulombic_energy_sq + (observables->coulombic_energy*observables->coulombic_energy) / ((double)counter);
	avg_observables->coulombic_energy_error = 0.5*sqrt(avg_observables->coulombic_energy_sq  - avg_observables->coulombic_energy*avg_observables->coulombic_energy);

	avg_observables->rd_energy = factor*avg_observables->rd_energy + observables->rd_energy / ((double)counter);
	avg_observables->rd_energy_sq = factor*avg_observables->rd_energy_sq + (observables->rd_energy*observables->rd_energy) / ((double)counter);
	avg_observables->rd_energy_error = 0.5*sqrt(avg_observables->rd_energy_sq  - avg_observables->rd_energy*avg_observables->rd_energy);

	avg_observables->polarization_energy = factor*avg_observables->polarization_energy + observables->polarization_energy / ((double)counter);
	avg_observables->polarization_energy_sq = factor*avg_observables->polarization_energy_sq + (observables->polarization_energy*observables->polarization_energy) / ((double)counter);
	avg_observables->polarization_energy_error = 0.5*sqrt(avg_observables->polarization_energy_sq  - avg_observables->polarization_energy*avg_observables->polarization_energy);

	avg_observables->vdw_energy = factor*avg_observables->vdw_energy + observables->vdw_energy / ((double)counter);
	avg_observables->vdw_energy_sq = factor*avg_observables->vdw_energy_sq + (observables->vdw_energy*observables->vdw_energy) / ((double)counter);
	avg_observables->vdw_energy_error = 0.5*sqrt(avg_observables->vdw_energy_sq  - avg_observables->vdw_energy*avg_observables->vdw_energy);



	avg_observables->dipole_rrms = factor*avg_observables->dipole_rrms + observables->dipole_rrms / ((double)counter);
	avg_observables->dipole_rrms_sq = factor*avg_observables->dipole_rrms_sq + (observables->dipole_rrms*observables->dipole_rrms) / ((double)counter);
	avg_observables->dipole_rrms_error = 0.5*sqrt(avg_observables->dipole_rrms_sq  - avg_observables->dipole_rrms*avg_observables->dipole_rrms);

	avg_observables->kinetic_energy = factor*avg_observables->kinetic_energy + observables->kinetic_energy / ((double)counter);
	avg_observables->kinetic_energy_sq = factor*avg_observables->kinetic_energy_sq + (observables->kinetic_energy*observables->kinetic_energy) / ((double)counter);
	avg_observables->kinetic_energy_error = 0.5*sqrt(avg_observables->kinetic_energy_sq  - avg_observables->kinetic_energy*avg_observables->kinetic_energy);

	avg_observables->temperature = factor*avg_observables->temperature + observables->temperature / ((double)counter);
	avg_observables->temperature_sq = factor*avg_observables->temperature_sq + (observables->temperature*observables->temperature) / ((double)counter);
	avg_observables->temperature_error = 0.5*sqrt(avg_observables->temperature_sq  - avg_observables->temperature*avg_observables->temperature);

	avg_observables->N = factor*avg_observables->N + observables->N / ((double)counter);
	avg_observables->N_sq = factor*avg_observables->N_sq + (observables->N*observables->N) / ((double)counter);
	avg_observables->N_error = 0.5*sqrt(avg_observables->N_sq  - avg_observables->N*avg_observables->N);

	avg_observables->spin_ratio = factor*avg_observables->spin_ratio + observables->spin_ratio / ((double)counter);
	avg_observables->spin_ratio_sq = factor*avg_observables->spin_ratio_sq + (observables->spin_ratio*observables->spin_ratio) / ((double)counter);
	avg_observables->spin_ratio_error = 0.5*sqrt(avg_observables->spin_ratio_sq  - avg_observables->spin_ratio*avg_observables->spin_ratio);

	avg_observables->NU = factor*avg_observables->NU + observables->NU / ((double)counter);


	/* avg in nodestats */
	avg_observables->boltzmann_factor = factor*avg_observables->boltzmann_factor + avg_nodestats->boltzmann_factor / ((double)counter);
	avg_observables->boltzmann_factor_sq = factor*avg_observables->boltzmann_factor_sq + avg_nodestats->boltzmann_factor_sq / ((double)counter);
	avg_observables->boltzmann_factor_error = 0.5*sqrt(avg_observables->boltzmann_factor_sq - avg_observables->boltzmann_factor*avg_observables->boltzmann_factor);


	avg_observables->acceptance_rate = factor*avg_observables->acceptance_rate + avg_nodestats->acceptance_rate / ((double)counter);
	avg_observables->acceptance_rate_insert = factor*avg_observables->acceptance_rate_insert + avg_nodestats->acceptance_rate_insert / ((double)counter);
	avg_observables->acceptance_rate_remove = factor*avg_observables->acceptance_rate_remove + avg_nodestats->acceptance_rate_remove / ((double)counter);
	avg_observables->acceptance_rate_displace = factor*avg_observables->acceptance_rate_displace + avg_nodestats->acceptance_rate_displace / ((double)counter);
	avg_observables->acceptance_rate_adiabatic = factor*avg_observables->acceptance_rate_adiabatic + avg_nodestats->acceptance_rate_adiabatic / ((double)counter);
	avg_observables->acceptance_rate_spinflip = factor*avg_observables->acceptance_rate_spinflip + avg_nodestats->acceptance_rate_spinflip / ((double)counter);

	avg_observables->cavity_bias_probability = factor*avg_observables->cavity_bias_probability + avg_nodestats->cavity_bias_probability / ((double)counter);
	avg_observables->cavity_bias_probability_sq = factor*avg_observables->cavity_bias_probability_sq + avg_nodestats->cavity_bias_probability_sq / ((double)counter);
	avg_observables->cavity_bias_probability_error = 0.5*sqrt(avg_observables->cavity_bias_probability_sq - avg_observables->cavity_bias_probability*avg_observables->cavity_bias_probability);

	avg_observables->polarization_iterations = factor*avg_observables->polarization_iterations + avg_nodestats->polarization_iterations / ((double)counter);
	avg_observables->polarization_iterations_sq = factor*avg_observables->polarization_iterations_sq + avg_nodestats->polarization_iterations_sq / ((double)counter);
	avg_observables->polarization_iterations_error = 0.5*sqrt(avg_observables->polarization_iterations_sq - avg_observables->polarization_iterations*avg_observables->polarization_iterations);


	/* get the mass of the two phases */
	for(molecule_ptr = system->molecules, particle_mass = 0, frozen_mass = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		if(molecule_ptr->frozen || molecule_ptr->adiabatic)
			frozen_mass += molecule_ptr->mass;
		else
			particle_mass = molecule_ptr->mass;
	}

	/* density in g/cm^3 */
	avg_observables->density = avg_observables->N*particle_mass/(system->pbc->volume*NA*A32CM3);
	avg_observables->density_error = avg_observables->N_error*particle_mass/(system->pbc->volume*NA*A32CM3);

	/* heat capacity in kJ/mol K */
	avg_observables->heat_capacity = (KB*NA/1000.0)*(avg_observables->energy_sq - avg_observables->energy*avg_observables->energy)/(system->temperature*system->temperature);
	avg_observables->heat_capacity_sq = factor*avg_observables->heat_capacity_sq + (avg_observables->heat_capacity*avg_observables->heat_capacity) / ((double)counter);
	avg_observables->heat_capacity_error = 0.5*sqrt(avg_observables->heat_capacity_sq - avg_observables->heat_capacity*avg_observables->heat_capacity);

	/* compressibility */
	avg_observables->compressibility = ATM2PASCALS*(system->pbc->volume/pow(METER2ANGSTROM, 3.0))*(avg_observables->N_sq - avg_observables->N*avg_observables->N)/(KB*system->temperature*avg_observables->N*avg_observables->N);
	avg_observables->compressibility_sq = factor*avg_observables->compressibility_sq + (avg_observables->compressibility*avg_observables->compressibility) / ((double)counter);
	avg_observables->compressibility_error = 0.5*sqrt(avg_observables->compressibility_sq - avg_observables->compressibility*avg_observables->compressibility);

	/* we have a solid phase */
	if(frozen_mass > 0.0) {

		/* percent weight */
		avg_observables->percent_wt = 100.0*avg_observables->N*particle_mass/(frozen_mass + avg_observables->N*particle_mass);
		avg_observables->percent_wt_error = 100.0*avg_observables->N_error*particle_mass/(frozen_mass + avg_observables->N_error*particle_mass);

		/* percent weight like ME*/
		avg_observables->percent_wt_me = 100.0*avg_observables->N*particle_mass/frozen_mass;
		avg_observables->percent_wt_me_error = 100.0*avg_observables->N_error*particle_mass/frozen_mass;

		/* excess weight mg/g */
		if(system->free_volume > 0.0) {

			if(system->fugacity > 0.0)
				avg_observables->excess_ratio = 1000.0*(avg_observables->N*particle_mass - (particle_mass*system->free_volume*system->fugacity*ATM2REDUCED)/system->temperature)/frozen_mass;
			else
				avg_observables->excess_ratio = 1000.0*(avg_observables->N*particle_mass - (particle_mass*system->free_volume*system->pressure*ATM2REDUCED)/system->temperature)/frozen_mass;
			avg_observables->excess_ratio_error = 1000.0*avg_observables->N_error*particle_mass/frozen_mass;


			/* pore density */
			avg_observables->pore_density = avg_observables->N*particle_mass/(system->free_volume*NA*A32CM3);
			avg_observables->pore_density_error = avg_observables->N_error*particle_mass/(system->free_volume*NA*A32CM3);

		}

		/* calculate the isosteric heat */
		avg_observables->qst = -(avg_observables->NU - avg_observables->N*avg_observables->energy);
		avg_observables->qst /= (avg_observables->N_sq - avg_observables->N*avg_observables->N);
		avg_observables->qst += system->temperature;
		avg_observables->qst *= KB*NA/1000.0;	/* convert to kJ/mol */

		avg_observables->qst_sq = factor*avg_observables->qst_sq + (avg_observables->qst*avg_observables->qst) / ((double)counter);
		avg_observables->qst_error = 0.5*sqrt(avg_observables->qst_sq - avg_observables->qst*avg_observables->qst);

	}


}


