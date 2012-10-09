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

	if(ns->accept_volume + ns->reject_volume)
		ns->acceptance_rate_volume = ((double)ns->accept_volume) / ((double)(ns->accept_volume + ns->reject_volume));
	else
		ns->acceptance_rate_volume = 0;

}

/* update node statistics related to the processing */
void update_nodestats(nodestats_t *nodestats, avg_nodestats_t *avg_nodestats) {

	static int counter = 0;
	double factor;

	++counter;
	factor = ((double)(counter - 1))/((double)(counter));

	avg_nodestats->boltzmann_factor = factor*avg_nodestats->boltzmann_factor 
		+ nodestats->boltzmann_factor / ((double)counter);
	avg_nodestats->boltzmann_factor_sq = factor*avg_nodestats->boltzmann_factor_sq 
		+ nodestats->boltzmann_factor*nodestats->boltzmann_factor / ((double)counter);

	/* these really aren't averages, but accumulative values */
	avg_nodestats->acceptance_rate = nodestats->acceptance_rate;
	avg_nodestats->acceptance_rate_insert = nodestats->acceptance_rate_insert;
	avg_nodestats->acceptance_rate_remove = nodestats->acceptance_rate_remove;
	avg_nodestats->acceptance_rate_displace = nodestats->acceptance_rate_displace;
	avg_nodestats->acceptance_rate_adiabatic = nodestats->acceptance_rate_adiabatic;
	avg_nodestats->acceptance_rate_spinflip = nodestats->acceptance_rate_spinflip;
	avg_nodestats->acceptance_rate_volume = nodestats->acceptance_rate_volume;

	avg_nodestats->cavity_bias_probability = factor*avg_nodestats->cavity_bias_probability 
		+ nodestats->cavity_bias_probability / ((double)counter);
	avg_nodestats->cavity_bias_probability_sq = factor*avg_nodestats->cavity_bias_probability_sq 
		+ nodestats->cavity_bias_probability*nodestats->cavity_bias_probability / ((double)counter);

	avg_nodestats->polarization_iterations = factor*avg_nodestats->polarization_iterations 
		+ nodestats->polarization_iterations / ((double)counter);
	avg_nodestats->polarization_iterations_sq = factor*avg_nodestats->polarization_iterations_sq 
		+ nodestats->polarization_iterations*nodestats->polarization_iterations / ((double)counter);

}

// update the stats for the individual sorbates
void update_sorbate_stats( system_t *system )
{

	static int counter = 0;
	int i;
	counter++;
	double factor = ((double)(counter-1))/((double)(counter));
	sorbateAverages_t *sorbate_ptr;
	sorbateAverages_t *sorb_ptr;
	molecule_t        *molecule_ptr;

	// Zero every count (N) for each sorbate in the averages list
	for( sorbate_ptr = system->sorbateStats.next; sorbate_ptr; sorbate_ptr = sorbate_ptr->next )
		sorbate_ptr->currentN = 0;

	// Count each sorbate in the system and record the total in
	// the corresponding entry in the sorbate averages list.
	for( molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) {
		for( sorbate_ptr = system->sorbateStats.next; sorbate_ptr; sorbate_ptr = sorbate_ptr->next ) {
			if( !strcasecmp( sorbate_ptr->id, molecule_ptr->moleculetype )) {
				sorbate_ptr->currentN++;
				break;
			}
		}	
	}
	
	// Using the current sorbate count (currentN), the number of readings
	// contributing to the average N (count) and the old average N (avgN)
	// (via "factor"), calculate the new average N (avgN). 
	for( sorbate_ptr = system->sorbateStats.next; sorbate_ptr; sorbate_ptr = sorbate_ptr->next ) {
		sorbate_ptr->avgN = factor * sorbate_ptr->avgN + ((double) sorbate_ptr->currentN)/((double)(counter));
	}


	// Calculate the weight percent and sorbed mass of each sorbate

	//     First, calculate the frozen mass of the system
	double system_mass = 0.0,
	       frozen_mass = 0.0;
	for( molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) { 
		// Frozen Mass:
		if(molecule_ptr->frozen || molecule_ptr->adiabatic)
			frozen_mass += molecule_ptr->mass;
	}

	for( sorbate_ptr = system->sorbateStats.next; sorbate_ptr; sorbate_ptr = sorbate_ptr->next ) {
		
		// sorbed_mass 
		double sorbed_mass =
		       sorbate_ptr->sorbed_mass = 
                       sorbate_ptr->avgN * sorbate_ptr->mass; 
		double R          = 0.8205746; //  L*atm/(mol*K)
		double A3PerLiter = 1.0e27;    //  A^3/L  (cubic angstroms per liter)
		
		// weight% 
		sorbate_ptr->percent_wt    = 100.0 * sorbed_mass/(frozen_mass + sorbed_mass);
		
		// weight% (ME) 
		sorbate_ptr->percent_wt_me = 100.0 * sorbed_mass/frozen_mass;

		// bulk mass -- mass in present in system-sized container under ideal conditions.
//		sorbate_ptr->bulk_mass =   (system->free_volume * system->pressure * sorbate_ptr->mass)  
//		                         / (system->temperature * R * A3PerLiter );

		// excess ratio
		sorbate_ptr->excess_ratio=1000*sorbate_ptr->mass * 
			(sorbate_ptr->avgN - sorbate_ptr->mass*system->free_volume*system->pressure*ATM2REDUCED/system->temperature);

		// density & pore density
		sorbate_ptr->density = sorbate_ptr->sorbed_mass/(system->avg_observables->volume*NA*A32CM3);
		sorbate_ptr->pore_density = sorbate_ptr->sorbed_mass/(system->free_volume*NA*A32CM3);



		// Selectivity
		// The denominator is temporarily collected in sorbate_ptr->selectivity, and is a sum of all the
		// avgN values for every sorbate in the system. Then, the avgN value for THIS sorbate is divided
		// by said total to yield selectivity for THIS sorbate. If the denominator is 0, selectivity is 
		// assigned a value of positive infinity.

		sorbate_ptr->selectivity = 0;
		for( sorb_ptr = system->sorbateStats.next; sorb_ptr; sorb_ptr = sorb_ptr->next ) {
			if(  strcasecmp(sorbate_ptr->id, sorb_ptr->id)  &&  (sorb_ptr->avgN != 0)  ) {
				sorbate_ptr->selectivity += sorb_ptr->avgN;
			} 
		}
		if( sorbate_ptr->selectivity == 0 )
			sorbate_ptr->selectivity = atof( "+Inf" );
		else 
			sorbate_ptr->selectivity = sorbate_ptr->avgN / sorbate_ptr->selectivity;
	}
}

void update_root_averages(system_t *system, observables_t *observables, avg_nodestats_t *avg_nodestats, avg_observables_t *avg_observables) {

	double particle_mass, frozen_mass, curr_density;
	sorbateAverages_t *sptr;

	molecule_t *molecule_ptr;
	static int counter = 0;
	double m, factor, gammaratio;

	++counter;
	m = (double)counter;
	factor = (m - 1.0)/m;

	/* the physical observables */
	avg_observables->energy = factor*avg_observables->energy 
		+ observables->energy / m;
	avg_observables->energy_sq = factor*avg_observables->energy_sq 
		+ (observables->energy*observables->energy) / m;
	avg_observables->energy_error = 0.5*sqrt(avg_observables->energy_sq  
		- avg_observables->energy*avg_observables->energy);

	avg_observables->energy_sq_sq = factor*avg_observables->energy_sq_sq 
		+ pow(observables->energy, 4) / m;
	avg_observables->energy_sq_error = 0.5*sqrt(avg_observables->energy_sq_sq  
		- pow(avg_observables->energy, 4));

	avg_observables->coulombic_energy = factor*avg_observables->coulombic_energy 
		+ observables->coulombic_energy / m;
	avg_observables->coulombic_energy_sq = factor*avg_observables->coulombic_energy_sq 
		+ (observables->coulombic_energy*observables->coulombic_energy) / m;
	avg_observables->coulombic_energy_error = 0.5*sqrt(avg_observables->coulombic_energy_sq  
		- avg_observables->coulombic_energy*avg_observables->coulombic_energy);

	avg_observables->rd_energy = factor*avg_observables->rd_energy 
		+ observables->rd_energy / m;
	avg_observables->rd_energy_sq = factor*avg_observables->rd_energy_sq 
		+ (observables->rd_energy*observables->rd_energy) / m;
	avg_observables->rd_energy_error = 0.5*sqrt(avg_observables->rd_energy_sq  
		- avg_observables->rd_energy*avg_observables->rd_energy);

	avg_observables->polarization_energy = factor*avg_observables->polarization_energy 
		+ observables->polarization_energy / m;
	avg_observables->polarization_energy_sq = factor*avg_observables->polarization_energy_sq 
		+ (observables->polarization_energy*observables->polarization_energy) / m;
	avg_observables->polarization_energy_error = 0.5*sqrt(avg_observables->polarization_energy_sq  
		- avg_observables->polarization_energy*avg_observables->polarization_energy);

	avg_observables->vdw_energy = factor*avg_observables->vdw_energy 
		+ observables->vdw_energy / m;
	avg_observables->vdw_energy_sq = factor*avg_observables->vdw_energy_sq 
		+ (observables->vdw_energy*observables->vdw_energy) / m;
	avg_observables->vdw_energy_error = 0.5*sqrt(avg_observables->vdw_energy_sq  
		- avg_observables->vdw_energy*avg_observables->vdw_energy);

	avg_observables->dipole_rrms = factor*avg_observables->dipole_rrms 
		+ observables->dipole_rrms / m;
	avg_observables->dipole_rrms_sq = factor*avg_observables->dipole_rrms_sq 
		+ (observables->dipole_rrms*observables->dipole_rrms) / m;
	avg_observables->dipole_rrms_error = 0.5*sqrt(avg_observables->dipole_rrms_sq  
		- avg_observables->dipole_rrms*avg_observables->dipole_rrms);

	avg_observables->kinetic_energy = factor*avg_observables->kinetic_energy 
		+ observables->kinetic_energy / m;
	avg_observables->kinetic_energy_sq = factor*avg_observables->kinetic_energy_sq 
		+ (observables->kinetic_energy*observables->kinetic_energy) / m;
	avg_observables->kinetic_energy_error = 0.5*sqrt(avg_observables->kinetic_energy_sq  
		- avg_observables->kinetic_energy*avg_observables->kinetic_energy);

	avg_observables->temperature = factor*avg_observables->temperature 
		+ observables->temperature / m;
	avg_observables->temperature_sq = factor*avg_observables->temperature_sq 
		+ (observables->temperature*observables->temperature) / m;
	avg_observables->temperature_error = 0.5*sqrt(avg_observables->temperature_sq 
		- avg_observables->temperature*avg_observables->temperature);

	avg_observables->volume = factor*avg_observables->volume 
		+ observables->volume / m;
	avg_observables->volume_sq = factor*avg_observables->volume_sq 
		+ (observables->volume*observables->volume) / m;
	avg_observables->volume_error = 0.5*sqrt(avg_observables->volume_sq  
		- avg_observables->volume*avg_observables->volume);

	avg_observables->N = factor*avg_observables->N 
		+ observables->N / m;
	avg_observables->N_sq = factor*avg_observables->N_sq 
		+ (observables->N*observables->N) / m;
	avg_observables->N_error = 0.5*sqrt(avg_observables->N_sq  
		- avg_observables->N*avg_observables->N);

	avg_observables->spin_ratio = factor*avg_observables->spin_ratio 
		+ observables->spin_ratio / m;
	avg_observables->spin_ratio_sq = factor*avg_observables->spin_ratio_sq 
		+ (observables->spin_ratio*observables->spin_ratio) / m;
	avg_observables->spin_ratio_error = 0.5*sqrt(avg_observables->spin_ratio_sq  
		- avg_observables->spin_ratio*avg_observables->spin_ratio);

	avg_observables->NU = factor*avg_observables->NU 
		+ observables->NU / m;

	/* avg in nodestats */
	avg_observables->boltzmann_factor = factor*avg_observables->boltzmann_factor 
		+ avg_nodestats->boltzmann_factor / m;
	avg_observables->boltzmann_factor_sq = factor*avg_observables->boltzmann_factor_sq 
		+ avg_nodestats->boltzmann_factor_sq / m;
	avg_observables->boltzmann_factor_error = 0.5*sqrt(avg_observables->boltzmann_factor_sq 
		- avg_observables->boltzmann_factor*avg_observables->boltzmann_factor);

	avg_observables->acceptance_rate = factor*avg_observables->acceptance_rate 
		+ avg_nodestats->acceptance_rate / m;
	avg_observables->acceptance_rate_insert = factor*avg_observables->acceptance_rate_insert 
		+ avg_nodestats->acceptance_rate_insert / m;
	avg_observables->acceptance_rate_remove = factor*avg_observables->acceptance_rate_remove 
		+ avg_nodestats->acceptance_rate_remove / m;
	avg_observables->acceptance_rate_displace = factor*avg_observables->acceptance_rate_displace 
		+ avg_nodestats->acceptance_rate_displace / m;
	avg_observables->acceptance_rate_adiabatic = factor*avg_observables->acceptance_rate_adiabatic 
		+ avg_nodestats->acceptance_rate_adiabatic / m;
	avg_observables->acceptance_rate_spinflip = factor*avg_observables->acceptance_rate_spinflip 
		+ avg_nodestats->acceptance_rate_spinflip / m;
	avg_observables->acceptance_rate_volume = factor*avg_observables->acceptance_rate_volume 
		+ avg_nodestats->acceptance_rate_volume / m;

	avg_observables->cavity_bias_probability = factor*avg_observables->cavity_bias_probability 
		+ avg_nodestats->cavity_bias_probability / m;
	avg_observables->cavity_bias_probability_sq = factor*avg_observables->cavity_bias_probability_sq 
		+ avg_nodestats->cavity_bias_probability_sq / m;
	avg_observables->cavity_bias_probability_error = 0.5*sqrt(avg_observables->cavity_bias_probability_sq 
		- avg_observables->cavity_bias_probability*avg_observables->cavity_bias_probability);

	avg_observables->polarization_iterations = factor*avg_observables->polarization_iterations 
		+ avg_nodestats->polarization_iterations / m;
	avg_observables->polarization_iterations_sq = factor*avg_observables->polarization_iterations_sq 
		+ avg_nodestats->polarization_iterations_sq / m;
	avg_observables->polarization_iterations_error = 0.5*sqrt(avg_observables->polarization_iterations_sq 
		- avg_observables->polarization_iterations*avg_observables->polarization_iterations);


	/* get the mass of the two phases */
	system->observables->total_mass = 0;
	for(molecule_ptr = system->molecules, particle_mass = 0, frozen_mass = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		if(molecule_ptr->frozen || molecule_ptr->adiabatic)
			frozen_mass += molecule_ptr->mass;
		else { 
			particle_mass = molecule_ptr->mass;
			system->observables->total_mass += particle_mass;
		}
	}
	system->observables->frozen_mass = frozen_mass; //store this for sorbate calculations
	system->observables->total_mass += frozen_mass;

	/* density in g/cm^3 */ //had to modify since density isn't neccessarily constant since adding NPT
	curr_density = observables->N * particle_mass / (system->pbc->volume*NA*A32CM3);

	avg_observables->density = factor*avg_observables->density 
		+ curr_density/m;
	avg_observables->density_sq = factor * avg_observables->density_sq 
		+ (curr_density*curr_density)/m;
	avg_observables->density_error = 0.5*sqrt(avg_observables->density_sq 
		- (avg_observables->density)*(avg_observables->density) );

	/* needed for calculating sstdev (stdev of stdev) */
	/*
	gammaratio = tgamma(0.5*(double)counter) / tgamma(0.5*((double)counter-1));
	gammaratio = sqrt(1.0/counter * ((double)counter - 1.0 - 2.0*gammaratio*gammaratio) );
	*/
	//stirling approx to avoid numerical overflow
	gammaratio = pow((m-2.0)/(m-1.0),0.5*m-1.0)*sqrt(0.5*(m-2.0))*exp(0.5);
	gammaratio = sqrt(1.0/counter * (m - 1.0 - 2.0*gammaratio*gammaratio) );
	
	/* heat capacity in kJ/mol K */
	avg_observables->heat_capacity = (KB*NA/1000.0)*(avg_observables->energy_sq 
		- avg_observables->energy*avg_observables->energy)/(system->temperature*system->temperature);
	/* error in heat capacity is the standard deviation of the variance, = 2*sstdev */
	avg_observables->heat_capacity_error = 2.0 * gammaratio * avg_observables->heat_capacity;

	/* compressibility */
	if ( system->ensemble != ENSEMBLE_NPT )
		avg_observables->compressibility = ATM2PASCALS*(system->pbc->volume/pow(METER2ANGSTROM, 3))*(avg_observables->N_sq 
			- avg_observables->N*avg_observables->N)/(KB*system->temperature*avg_observables->N*avg_observables->N);
	else 
		avg_observables->compressibility = ATM2PASCALS * pow(METER2ANGSTROM,-3) *
			( avg_observables->volume_sq - avg_observables->volume * avg_observables->volume ) / 
			( KB * system->temperature * avg_observables->volume );
	avg_observables->compressibility_error = 2.0 * gammaratio * avg_observables->compressibility;

	/* we have a solid phase */
	if(frozen_mass > 0.0) {

		/* percent weight */
		avg_observables->percent_wt = 100.0*avg_observables->N*particle_mass/(frozen_mass 
			+ avg_observables->N*particle_mass);
		avg_observables->percent_wt_error = 100.0*avg_observables->N_error*particle_mass/(frozen_mass 
			+ avg_observables->N_error*particle_mass);

		/* percent weight like ME*/
		avg_observables->percent_wt_me = 100.0*avg_observables->N*particle_mass/frozen_mass;
		avg_observables->percent_wt_me_error = 100.0*avg_observables->N_error*particle_mass/frozen_mass;

		/* excess weight mg/g */
		if(system->free_volume > 0.0) {

			if(system->fugacities)
				avg_observables->excess_ratio = 1000.0*(avg_observables->N*particle_mass 
					- (particle_mass*system->free_volume*system->fugacities[0]*ATM2REDUCED)/system->temperature)/frozen_mass;
			else
				avg_observables->excess_ratio = 1000.0*(avg_observables->N*particle_mass 
					- (particle_mass*system->free_volume*system->pressure*ATM2REDUCED)/system->temperature)/frozen_mass;
			avg_observables->excess_ratio_error = 1000.0*avg_observables->N_error*particle_mass/frozen_mass;


			/* pore density */ //only valid for constant V, pure systems
			avg_observables->pore_density = curr_density * system->pbc->volume / system->free_volume;
			avg_observables->pore_density_error = avg_observables->N_error*particle_mass/(system->free_volume*NA*A32CM3);

		}

		/* calculate the isosteric heat */
		avg_observables->qst = -(avg_observables->NU 
			- avg_observables->N*avg_observables->energy);
		avg_observables->qst /= (avg_observables->N_sq 
			- avg_observables->N*avg_observables->N);
		avg_observables->qst += system->temperature;
		avg_observables->qst *= KB*NA/1000.0;	/* convert to kJ/mol */

	}


}


