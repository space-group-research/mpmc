/* 

Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

//convert string *a to a double and store at the address of f.
int safe_atof ( char * a, double * f ) {
	if ( sscanf(a,"%lf",f) == 0 ) return 1; //failure
	return 0; //success
}
	
//convert string *a to an int and store at the address of f.
int safe_atoi ( char * a, int * i ) {
	if ( sscanf(a,"%d",i) == 0 ) return 1; //failure
	return 0; //success
}

//convert string *a to an unsigned int and store at the address of f.
int safe_atou ( char * a, uint32_t * i ) {
	if ( sscanf(a,"%u",i) == 0 ) return 1; //failure
	return 0; //success
}

//converts string *a to a long
int safe_atol ( char * a, long unsigned int * l ) {
	if ( sscanf(a,"%lu",l) == 0 ) return 1; //fail
	return 0;
}


/* check each input command and set system flags */
int do_command (system_t * system, char ** token ) {
	int i;

	// check for comment/blanks
	if (!strncasecmp(token[0], "!", 1)) return 0;
	else if (!strncasecmp(token[0],"#", 1)) return 0;
	else if (!strcasecmp(token[0],"")) return 0;

	// ensemble options
	else if (!strcasecmp(token[0],"ensemble")) {
		if (!strcasecmp(token[1],"nvt"))
			system->ensemble = ENSEMBLE_NVT;
		else if (!strcasecmp(token[1],"uvt"))
			system->ensemble = ENSEMBLE_UVT;
		else if (!strcasecmp(token[1],"surf"))
			system->ensemble = ENSEMBLE_SURF;
		else if (!strcasecmp(token[1],"surf_fit"))
			system->ensemble = ENSEMBLE_SURF_FIT;
		else if (!strcasecmp(token[1],"nve"))
			system->ensemble = ENSEMBLE_NVE;
		else if (!strcasecmp(token[1],"total_energy")) 
			system->ensemble = ENSEMBLE_TE;
		else if (!strcasecmp(token[1],"npt"))
			system->ensemble = ENSEMBLE_NPT;
		else if (!strcasecmp(token[1],"replay"))
			system->ensemble = ENSEMBLE_REPLAY;
	}

	// random seed options
	else if (!strcasecmp(token[0],"preset_seeds")) {
		{ if ( safe_atou(token[1],&(system->preset_seeds[0])) ) return 1; }
		{ if ( safe_atou(token[2],&(system->preset_seeds[1])) ) return 1; }
		{ if ( safe_atou(token[3],&(system->preset_seeds[2])) ) return 1; }
		{ if ( safe_atou(token[4],&(system->preset_seeds[3])) ) return 1; }
		system->preset_seeds_on = 1;
	}

	//deprecated seed option
	else if (!strcasecmp(token[0],"seed")) {
		error( "INPUT: seed is deprecated\n" );
		error( "INPUT: use \"preset_seeds <int> <int> <int> <int>\" to *manually* seed the rng.\n" );
		die(-1);
	}



	// surf options
	////////////////////////
	

	// Option for fitting against arbitrary configurations, VS the default behavior of fitting
	// against a small set of orientations, while only varying their separation distance.
	else if(!strcasecmp(token[0], "fit_arbitrary_configs")) {
	    
		if(!strcasecmp(token[1], "on")) {
			system->surf_fit_arbitrary_configs = 1; 
		}
		else if (!strcasecmp(token[1], "off")) {
			system->surf_fit_arbitrary_configs = 0; 
		}
		else return 1;
	}
	
	else if(!strcasecmp(token[0],"surf_decomp")) {
		if(!strcasecmp(token[1], "on" ))
			system->surf_decomp = 1;
		else if (!strcasecmp(token[1], "off" ))
			system->surf_decomp = 0;
		else return 1; //unrecognized argument
	}
	else if(!strcasecmp(token[0], "surf_min" ))
		{ if ( safe_atof(token[1],&(system->surf_min)) ) return 1; }
	else if(!strcasecmp(token[0], "surf_max" ))
		{ if ( safe_atof(token[1],&(system->surf_max)) ) return 1; }
	else if(!strcasecmp(token[0], "surf_inc" ))
		{ if ( safe_atof(token[1],&(system->surf_inc)) ) return 1; }
	else if(!strcasecmp(token[0], "surf_ang" ))
		{ if ( safe_atof(token[1],&(system->surf_ang)) ) return 1; }
	else if(!strcasecmp(token[0], "surf_print_level" ))
		{ if ( safe_atoi(token[1],&(system->surf_print_level)) ) return 1; }
	//allows us to specify the surf-fit scales in the input file
	else if(!strcasecmp(token[0], "surf_weight_constant")) {
		{ if ( safe_atof(token[1],&(system->surf_weight_constant)) ) return 1; }
		system->surf_weight_constant_on = 1;
	}
	else if(!strcasecmp(token[0], "surf_scale_q")) {
		{ if ( safe_atof(token[1],&(system->surf_scale_q)) ) return 1; }
		system->surf_scale_q_on = 1;
	}
	else if(!strcasecmp(token[0], "surf_scale_r")) {
		{ if ( safe_atof(token[1],&(system->surf_scale_r)) ) return 1; }
		system->surf_scale_r_on = 1;
	}
	else if(!strcasecmp(token[0], "surf_scale_epsilon")) {
		{ if ( safe_atof(token[1],&(system->surf_scale_epsilon)) ) return 1; }
		system->surf_scale_epsilon_on = 1;
	}
	else if(!strcasecmp(token[0], "surf_scale_sigma")) {
		{ if ( safe_atof(token[1],&(system->surf_scale_sigma)) ) return 1; }
		system->surf_scale_sigma_on = 1;
	}
	else if(!strcasecmp(token[0], "surf_scale_omega")) {
		{ if ( safe_atof(token[1],&(system->surf_scale_omega)) ) return 1; }
		system->surf_scale_omega_on = 1;
	}
	else if(!strcasecmp(token[0], "surf_scale_alpha")) {
		{ if ( safe_atof(token[1],&(system->surf_scale_alpha)) ) return 1; }
		system->surf_scale_alpha_on = 1;
	}
	else if(!strcasecmp(token[0], "surf_scale_pol")) {
		{ if ( safe_atof(token[1],&(system->surf_scale_pol)) ) return 1; }
		system->surf_scale_pol_on = 1;
	}
	else if(!strcasecmp(token[0], "surf_qshift")) {
		if (!strcasecmp(token[1],"on")) {
			system->surf_qshift_on = 1;
			printf("INPUT: surf_qshift is ON.\n");
			printf("INPUT: only use qshift with x-axis aligned linear molecules.\n");
		}
		else if (!strcasecmp(token[1],"off")) {
			system->surf_qshift_on = 0;
			printf("INPUT: surf_qshift is OFF.\n");
		}
		else return 1;
	}
	else if(!strcasecmp(token[0], "surf_preserve")) {
		if(!strcasecmp(token[1],"on"))
			system->surf_preserve = 1;
		else if(!strcasecmp(token[1],"off"))
			system->surf_preserve = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "surf_preserve_rotation")) {
		if(system->surf_preserve_rotation_on != NULL) {
			error("INPUT: surf_preserve_rotationalready set.\n");
			return 1;
		}
		system->surf_preserve_rotation_on = malloc(6*sizeof(double));
		{ if ( safe_atof(token[1],&(system->surf_preserve_rotation_on->alpha1)) ) return 1; }
		{ if ( safe_atof(token[2],&(system->surf_preserve_rotation_on->beta1)) ) return 1; }
		{ if ( safe_atof(token[3],&(system->surf_preserve_rotation_on->gamma1)) ) return 1; }
		{ if ( safe_atof(token[4],&(system->surf_preserve_rotation_on->alpha2)) ) return 1; }
		{ if ( safe_atof(token[5],&(system->surf_preserve_rotation_on->beta2)) ) return 1; }
		{ if ( safe_atof(token[6],&(system->surf_preserve_rotation_on->gamma2)) ) return 1; }
	}
	else if (!strcasecmp(token[0], "surf_output")) {
		if(!system->surf_output) {
			system->surf_output = calloc(MAXLINE,sizeof(char));
			memnullcheck(system->surf_output,MAXLINE*sizeof(char),__LINE__-1, __FILE__);
			strcpy(system->surf_output,token[1]);
		} else return 1;
	}
	else if(!strcasecmp(token[0], "surf_global_axis")) {
		if (!strcasecmp(token[1],"on")) {
			system->surf_global_axis_on = 1;
			printf("INPUT: surf_global_axis is ON.\n");
		}
		else if (!strcasecmp(token[1],"off")) {
			system->surf_global_axis_on = 0;
			printf("INPUT: surf_global_axis is OFF.\n");
		}
		else return 1;
	}
	else if(!strcasecmp(token[0], "surf_descent")) {
		if (!strcasecmp(token[1],"on")) {
			system->surf_descent = 1;
			printf("INPUT: surf_descent is ON.\n");
		}
		else if (!strcasecmp(token[1],"off")) {
			system->surf_descent = 0;
			printf("INPUT: surf_descent is OFF.\n");
		}
		else return 1;
	}

	//spectre options
	else if(!strcasecmp(token[0], "spectre")) {
		if(!strcasecmp(token[1], "on" ))
			system->spectre = 1;
		else if(!strcasecmp(token[1], "off" )) 
			system->spectre = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "spectre_max_charge"))
		system->spectre_max_charge = fabs(atof(token[1]));
	else if(!strcasecmp(token[0], "spectre_max_target"))
		system->spectre_max_target = fabs(atof(token[1]));

	//cavity options
	else if (!strcasecmp(token[0], "cavity_bias")) {
		if (!strcasecmp(token[1], "on"))
			system->cavity_bias = 1;
		else if (!strcasecmp(token[1], "off"))
			system->cavity_bias = 0;
		else return 1; //no match
	}
	else if (!strcasecmp(token[0], "cavity_grid"))
		{ if ( safe_atoi(token[1],&(system->cavity_grid_size)) ) return 1; }
	else if (!strcasecmp(token[0],"cavity_radius"))
		{ if ( safe_atof(token[1],&(system->cavity_radius)) ) return 1; }
	else if (!strcasecmp(token[0],"cavity_autoreject")) {
		if (!strcasecmp(token[1], "on"))
			system->cavity_autoreject = 1;
		else if (!strcasecmp(token[1], "off"))
			system->cavity_autoreject = 0;
		else return 1; //no match
	}
	else if (!strcasecmp(token[0],"cavity_autoreject_absolute")) {
		if (!strcasecmp(token[1], "on"))
			system->cavity_autoreject_absolute = 1;
		else if (!strcasecmp(token[1], "off"))
			system->cavity_autoreject_absolute = 0;
		else return 1; //no match
	}
	else if (!strcasecmp(token[0],"cavity_autoreject_scale")) {
		{ if ( safe_atof(token[1],&(system->cavity_autoreject_scale)) ) return 1; }
	}

	//polar options
	else if(!strcasecmp(token[0], "polarization")) {
		if(!strcasecmp(token[1], "on"))
			system->polarization = 1; 
		else if (!strcasecmp(token[1], "off"))
			system->polarization = 0; 
		else return 1;
	}
	else if(!strcasecmp(token[0], "polarvdw") || !strcasecmp(token[0], "cdvdw")) {
		if(!strcasecmp(token[1], "on")) {
			system->polarvdw = 1; 
			system->polarization=1;
			system->polar_iterative=1; //matrix inversion destroys A_matrix before vdw can use it.
			output("INPUT: Forcing polar_iterative ON for CP-VdW.\n");
		}
		else if (!strcasecmp(token[1], "evects")) {
			system->polarvdw = 2; //calculate eigenvectors
			system->polarization=1;
			system->polar_iterative=1; //matrix inversion destroys A_matrix before vdw can use it.
			output("INPUT: Forcing polar_iterative ON for CP-VdW.\n");
		}
		else if ( !strcasecmp(token[1], "comp")) {
			system->polarvdw = 3; //calculate eigenvectors
			system->polarization=1;
			system->polar_iterative=1; //matrix inversion destroys A_matrix before vdw can use it.
			output("INPUT: Forcing polar_iterative ON for CP-VdW.\n");
		}
		else if (!strcasecmp(token[1], "off")) 
			system->polarvdw = 0;
		else return 1;
	}
	else if (!strcasecmp(token[0], "cdvdw_9th_repulsion")) {
		if (!strcasecmp(token[1], "on"))
			system->cdvdw_9th_repulsion = 1;
		else if (!strcasecmp(token[1], "off"))
			system->cdvdw_9th_repulsion = 0;
		else return 1;
	}
	else if (!strcasecmp(token[0], "cdvdw_exp_repulsion")) {
		if (!strcasecmp(token[1], "on"))
			system->cdvdw_exp_repulsion = 1;
		else if (!strcasecmp(token[1], "off"))
			system->cdvdw_exp_repulsion = 0;
		else return 1;
	}
	else if (!strcasecmp(token[0], "cdvdw_sig_repulsion")) {
		if (!strcasecmp(token[1], "on"))
			system->cdvdw_sig_repulsion = 1;
		else if (!strcasecmp(token[1], "off"))
			system->cdvdw_sig_repulsion = 0;
		else return 1;
	}
	
	else if (!strcasecmp(token[0], "polar_ewald_full")) {
		if (!strcasecmp(token[1], "on"))
			system->polar_ewald_full = 1;
		else if (!strcasecmp(token[1], "off"))
			system->polar_ewald_full = 0;
		else return 1;
	}
	else if (!strcasecmp(token[0], "polar_ewald")) {
		if (!strcasecmp(token[1], "on"))
			system->polar_ewald = 1;
		else if (!strcasecmp(token[1], "off"))
			system->polar_ewald = 0;
		else return 1;
	}
	//polar wolf shiz
	else if (!strcasecmp(token[0], "polar_wolf_full")) {
		if (!strcasecmp(token[1], "on"))
			system->polar_wolf_full = 1;
		else if (!strcasecmp(token[1], "off"))
			system->polar_wolf_full = 0;
		else return 1;
	}
	else if (!strcasecmp(token[0], "polar_wolf")) {
		if (!strcasecmp(token[1], "on"))
			system->polar_wolf = 1;
		else if (!strcasecmp(token[1], "off"))
			system->polar_wolf = 0;
		else return 1;
	}
	else if (!strcasecmp(token[0], "polar_wolf_alpha_lookup")) {
		if (!strcasecmp(token[1], "on"))
			system->polar_wolf_alpha_lookup = 1;
		else if (!strcasecmp(token[1], "off"))
			system->polar_wolf_alpha_lookup = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "polar_wolf_damp")) 
		{ if ( safe_atof(token[1],&(system->polar_wolf_alpha)) ) return 1; }
	else if(!strcasecmp(token[0], "polar_wolf_alpha"))  //same as polar_wolf_damp
		{ if ( safe_atof(token[1],&(system->polar_wolf_alpha)) ) return 1; }
	else if(!strcasecmp(token[0], "polar_wolf_alpha_lookup_cutoff"))  //store shit in lookup table
		{ if ( safe_atof(token[1],&(system->polar_wolf_alpha_lookup_cutoff)) ) return 1; }

	// replay options
	else if (!strcasecmp(token[0], "calc_pressure")) {
		if (!strcasecmp(token[1], "on"))
			system->calc_pressure = 1;
		else if (!strcasecmp(token[1], "off"))
			system->calc_pressure = 0;
		else return 1;
	}
	else if (!strcasecmp(token[0], "calc_pressure_dv")) 
		{ if ( safe_atof(token[1],&(system->calc_pressure_dv)) ) return 1; }

	/*set total energy for NVE*/
	else if(!strcasecmp(token[0], "total_energy"))
		{ if ( safe_atof(token[1],&(system->total_energy)) ) return 1; }

	else if(!strcasecmp(token[0], "numsteps"))
		{ if ( safe_atoi(token[1],&(system->numsteps)) ) return 1; }

	else if(!strcasecmp(token[0], "corrtime"))
		{ if ( safe_atoi(token[1],&(system->corrtime)) ) return 1; }
	
	/* set Monte Carlo options */
	else if(!strcasecmp(token[0], "move_probability")) 
		{ if ( safe_atof(token[1],&(system->move_probability)) ) return 1; }
	else if(!strcasecmp(token[0], "gwp_probability")) 
		{ if ( safe_atof(token[1],&(system->gwp_probability)) ) return 1; }
	else if(!strcasecmp(token[0], "rot_probability")) 
		{ if ( safe_atof(token[1],&(system->rot_probability)) ) return 1; }
	else if(!strcasecmp(token[0], "insert_probability")) 
		{ if ( safe_atof(token[1],&(system->insert_probability)) ) return 1; }
	else if(!strcasecmp(token[0], "adiabatic_probability")) 
		{ if ( safe_atof(token[1],&(system->adiabatic_probability)) ) return 1; }
	else if(!strcasecmp(token[0], "spinflip_probability")) 
		{ if ( safe_atof(token[1],&(system->spinflip_probability)) ) return 1; }
	else if(!strcasecmp(token[0], "volume_probability")) 
		{ if ( safe_atof(token[1],&(system->volume_probability)) ) return 1; }
	else if(!strcasecmp(token[0], "volume_change_factor")) 
		{ if ( safe_atof(token[1],&(system->volume_change_factor)) ) return 1; }
	else if(!strcasecmp(token[0], "ptemp_freq")) 
		{ if ( safe_atoi(token[1],&(system->ptemp_freq)) ) return 1; }
	/*end setting MC options*/

	/* parallel tempering options */
	else if(!strcasecmp(token[0], "parallel_tempering")) {
		if(!strcasecmp(token[1], "on")) 
			system->parallel_tempering = 1;
		else if(!strcasecmp(token[1], "off")) 
			system->parallel_tempering = 0;
		else return 1; //no match
	}
	else if(!strcasecmp(token[0], "max_temperature")) 
		{ if ( safe_atof(token[1],&(system->max_temperature)) ) return 1; }

	else if(!strcasecmp(token[0], "temperature")) 
		{ if ( safe_atof(token[1],&(system->temperature)) ) return 1; }

	else if(!strcasecmp(token[0], "simulated_annealing")) {
		if(!strcasecmp(token[1], "on")) 
			system->simulated_annealing = 1;
		else if(!strcasecmp(token[1], "off")) 
			system->simulated_annealing = 0;
		else return 1; //no match
	}

	else if(!strcasecmp(token[0], "simulated_annealing_schedule")) 
		{ if ( safe_atof(token[1],&(system->simulated_annealing_schedule)) ) return 1; }

	else if(!strcasecmp(token[0], "simulated_annealing_target")) 
		{ if ( safe_atof(token[1],&(system->simulated_annealing_target)) ) return 1; }

	else if(!strcasecmp(token[0], "pressure")) 
		{ if ( safe_atof(token[1],&(system->pressure)) ) return 1; }

/* fugacity shits */
	else if(!strcasecmp(token[0], "h2_fugacity")) {
		if(!strcasecmp(token[1], "on"))
			system->h2_fugacity = 1;
		else if(!strcasecmp(token[1], "off"))
			system->h2_fugacity = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "co2_fugacity")) {
		if(!strcasecmp(token[1], "on"))
			system->co2_fugacity = 1;
		else if(!strcasecmp(token[1], "off"))
			system->co2_fugacity = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "ch4_fugacity")) {
		if(!strcasecmp(token[1], "on"))
			system->ch4_fugacity = 1;
		else if(!strcasecmp(token[1], "off"))
			system->ch4_fugacity = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "n2_fugacity")) {
		if(!strcasecmp(token[1], "on"))
			system->n2_fugacity = 1;
		else if(!strcasecmp(token[1],"off"))
			system->n2_fugacity = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "user_fugacities")) { //read in a list of user-defined fugacities
		system->user_fugacities = 1;
		for ( i=0; strlen(token[i+1]) > 0; i++ );
		if ( i == 0 ) return 1;
		system->fugacitiesCount = i;
		system->fugacities = calloc(i,sizeof(double));
		memnullcheck(system->fugacities,i*sizeof(double),__LINE__-1, __FILE__);
		for ( i=0; strlen(token[i+1]) > 0; i++ )
			if ( safe_atof(token[i+1], &(system->fugacities[i])) )
				return 1;
	}

	else if(!strcasecmp(token[0], "free_volume"))
		{ if ( safe_atof(token[1],&(system->free_volume)) ) return 1; }

	else if(!strcasecmp(token[0], "rd_only")) {
		if(!strcasecmp(token[1], "on"))
			system->rd_only = 1;
		else if (!strcasecmp(token[1], "off"))
			system->rd_only = 0;
		else return 1;
	}

	else if(!strcasecmp(token[0], "gwp")) {
		if(!strcasecmp(token[1], "on"))
			system->gwp = 1;
		else if(!strcasecmp(token[1], "off"))
			system->gwp = 0;
		else return 1;
	}

	else if(!strcasecmp(token[0], "wolf")) {
		if(!strcasecmp(token[1], "on"))
			system->wolf = 1;
		else if(!strcasecmp(token[1], "off"))
			system->wolf = 0;
		else return 1;
	}

	// rd options
	else if(!strcasecmp(token[0], "rd_lrc")) {
		if(!strcasecmp(token[1], "on"))
			system->rd_lrc = 1;
		else if(!strcasecmp(token[1], "off"))
			system->rd_lrc = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "rd_crystal")) {
		if(!strcasecmp(token[1], "on"))
			system->rd_crystal = 1;
		else if(!strcasecmp(token[1], "off"))
			system->rd_crystal = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "rd_crystal_order"))
		{ if ( safe_atoi(token[1],&(system->rd_crystal_order)) ) return 1; }
	else if(!strcasecmp(token[0], "rd_anharmonic")) {
		if(!strcasecmp(token[1], "on"))
			system->rd_anharmonic = 1;
		else if(!strcasecmp(token[1], "off"))
			system->rd_anharmonic = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "rd_anharmonic_k"))
		{ if ( safe_atof(token[1],&(system->rd_anharmonic_k)) ) return 1; }
	else if(!strcasecmp(token[0], "rd_anharmonic_g"))
		{ if ( safe_atof(token[1],&(system->rd_anharmonic_g)) ) return 1; }

	else if(!strcasecmp(token[0], "feynman_hibbs")) {
		if(!strcasecmp(token[1],"on"))
			system->feynman_hibbs = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->feynman_hibbs = 0;
		else return 1;
	}
	//shitty feynman-hibbs correction for polarvdw (the default is better)
	else if(!strcasecmp(token[0], "vdw_fh_2be")) {
		if(!strcasecmp(token[1],"on"))
			system->vdw_fh_2be = 1;
		else if (!strcasecmp(token[1],"off"))
			system->vdw_fh_2be = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "feynman_kleinert")) {
		if(!strcasecmp(token[1],"on"))
			system->feynman_kleinert = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->feynman_kleinert = 0;
		else return 1;
	}
	
	else if(!strcasecmp(token[0], "feynman_hibbs_order"))
		{ if ( safe_atoi(token[1],&(system->feynman_hibbs_order)) ) return 1; }

	else if(!strcasecmp(token[0], "sg")) {
		if(!strcasecmp(token[1],"on"))
			system->sg = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->sg = 0;
		else return 1;
	}

	else if(!strcasecmp(token[0], "waldmanhagler")) {
		if(!strcasecmp(token[1],"on"))
			system->waldmanhagler = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->waldmanhagler = 0;
		else return 1;
	}

	else if(!strcasecmp(token[0], "halgren_mixing")) {
		if(!strcasecmp(token[1],"on"))
			system->halgren_mixing = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->halgren_mixing = 0;
		else return 1;
	}
	
	else if(!strcasecmp(token[0], "dreiding")) {
		if(!strcasecmp(token[1],"on"))
			system->dreiding = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->dreiding = 0;
		else return 1;
	}

	else if(!strcasecmp(token[0], "lj_buffered_14_7")) {
		if(!strcasecmp(token[1],"on"))
			system->lj_buffered_14_7 = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->lj_buffered_14_7 = 0;
		else return 1;
	}
	
	else if(!strcasecmp(token[0], "wrapall")) {
		if(!strcasecmp(token[1],"on"))
			system->wrapall = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->wrapall = 0;
		else return 1;
	}
	
	else if(!strcasecmp(token[0], "scale_charge"))
		{ if ( safe_atof(token[1],&(system->scale_charge)) ) return 1; }

	else if(!strcasecmp(token[0], "ewald_alpha")) { 
		if ( safe_atof(token[1],&(system->ewald_alpha)) ) return 1; 
		system->ewald_alpha_set = 1;
	}
	
	else if(!strcasecmp(token[0], "ewald_kmax"))
		{ if ( safe_atoi(token[1],&(system->ewald_kmax)) ) return 1; }
	
	else if(!strcasecmp(token[0], "pbc_cutoff"))
		{ if ( safe_atof(token[1],&(system->pbc->cutoff)) ) return 1; }

//polar options
	else if(!strcasecmp(token[0], "polar_ewald")) {
		if(!strcasecmp(token[1],"on"))
			system->polar_ewald = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->polar_ewald = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "polar_ewald_alpha")) {
		if ( safe_atof(token[1],&(system->polar_ewald_alpha)) ) return 1;
		system->polar_ewald_alpha_set = 1;
	}

	else if(!strcasecmp(token[0], "polarizability_tensor")) {
		if(!strcasecmp(token[1],"on"))
			system->polarizability_tensor = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->polarizability_tensor = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "polar_zodid")) {
		if(!strcasecmp(token[1],"on"))
			system->polar_zodid = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->polar_zodid = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "polar_iterative")) {
		if(!strcasecmp(token[1],"on"))
			system->polar_iterative = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->polar_iterative = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "polar_palmo")) {
		if(!strcasecmp(token[1],"on"))
			system->polar_palmo = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->polar_palmo = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "polar_gs")) {
		if(!strcasecmp(token[1],"on"))
			system->polar_gs = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->polar_gs = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "polar_gs_ranked")) {
		if(!strcasecmp(token[1],"on"))
			system->polar_gs_ranked = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->polar_gs_ranked = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "polar_sor")) {
		if(!strcasecmp(token[1],"on"))
			system->polar_sor = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->polar_sor = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "polar_esor")) {
		if(!strcasecmp(token[1],"on"))
			system->polar_esor = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->polar_esor = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "polar_gamma"))
		{ if ( safe_atof(token[1],&(system->polar_gamma)) ) return 1; }
	else if(!strcasecmp(token[0], "polar_damp"))
		{ if ( safe_atof(token[1],&(system->polar_damp)) ) return 1; }
	else if(!strcasecmp(token[0], "polar_precision"))
		{ if ( safe_atof(token[1],&(system->polar_precision)) ) return 1; }
	else if(!strcasecmp(token[0], "polar_max_iter"))
		{ if ( safe_atoi(token[1],&(system->polar_max_iter)) ) return 1; }
	else if(!strcasecmp(token[0], "polar_damp_type")) {
		if(!strcasecmp(token[1],"none"))
			system->damp_type = DAMPING_OFF;
		else if(!strcasecmp(token[1],"off"))
			system->damp_type = DAMPING_OFF;
		else if(!strcasecmp(token[1],"linear"))
			system->damp_type = DAMPING_LINEAR;
		else if (!strcasecmp(token[1],"exponential")) 
			system->damp_type = DAMPING_EXPONENTIAL;
		else return 1;
	}
	else if(!strcasecmp(token[0], "polar_rrms")) {
		if(!strcasecmp(token[1],"on"))
			system->polar_rrms = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->polar_rrms = 0;
		else return 1;
	}
		
	else if(!strcasecmp(token[0], "cuda")) {
		if(!strcasecmp(token[1],"on"))
			system->cuda = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->cuda = 0;
		else return 1;
	}
	
	else if(!strcasecmp(token[0], "independent_particle")) {
		if(!strcasecmp(token[1],"on"))
			system->independent_particle = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->independent_particle = 0;
		else return 1;
	}

//quantum rotation stuff
#ifdef QM_ROTATION
	else if(!strcasecmp(token[0], "quantum_rotation")) {
		if(!strcasecmp(token[1],"on"))
			system->quantum_rotation = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->quantum_rotation = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "quantum_rotation_hindered")) {
		if(!strcasecmp(token[1],"on"))
			system->quantum_rotation_hindered = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->quantum_rotation_hindered = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "quantum_rotation_hindered_barrier"))
		{ if ( safe_atof(token[1],&(system->quantum_rotation_hindered_barrier)) ) return 1; }
	else if(!strcasecmp(token[0], "quantum_rotation_B"))
		{ if ( safe_atof(token[1],&(system->quantum_rotation_B)) ) return 1; }
	else if(!strcasecmp(token[0], "quantum_rotation_level_max"))
		{ if ( safe_atoi(token[1],&(system->quantum_rotation_level_max)) ) return 1; }
	else if(!strcasecmp(token[0], "quantum_rotation_l_max"))
		{ if ( safe_atoi(token[1],&(system->quantum_rotation_l_max)) ) return 1; }
	else if(!strcasecmp(token[0], "quantum_rotation_sum"))
		{ if ( safe_atoi(token[1],&(system->quantum_rotation_sum)) ) return 1; }
#endif //end QM rotation

/* #ifdef XXX
	else if(!strcasecmp(token[0], "quantum_vibration")) {
		if(!strcasecmp(token[1],"on"))
			system->quantum_vibration = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->quantum_vibration = 0;
		else return 1;
	}
#endif */

	// set job name (CRC)
	else if (!strcasecmp(token[0], "job_name")) {
		//already allocated
		strcpy(system->job_name,token[1]);
	}

	else if (!strcasecmp(token[0], "pqr_input")) {
		if(!system->pqr_input) {
			system->pqr_input = calloc(MAXLINE,sizeof(char));
			memnullcheck(system->pqr_input,MAXLINE*sizeof(char),__LINE__-1, __FILE__);
			strcpy(system->pqr_input,token[1]);
		} else return 1;
	}
	else if (!strcasecmp(token[0], "pqr_output")) {
		if(!system->pqr_output) {
			system->pqr_output = calloc(MAXLINE,sizeof(char));
			memnullcheck(system->pqr_output,MAXLINE*sizeof(char),__LINE__-1, __FILE__);
			strcpy(system->pqr_output,token[1]);
		} else return 1;
	}
	else if (!strcasecmp(token[0], "pqr_restart")) {
		if(!system->pqr_restart) {
			system->pqr_restart = calloc(MAXLINE,sizeof(char));
			memnullcheck(system->pqr_restart,MAXLINE*sizeof(char),__LINE__-1, __FILE__);
			strcpy(system->pqr_restart,token[1]);
		} else return 1;
	}
	else if (!strcasecmp(token[0], "traj_output")) {
		if(!system->traj_output) {
			system->traj_output = calloc(MAXLINE,sizeof(char));
			memnullcheck(system->traj_output,MAXLINE*sizeof(char),__LINE__-1, __FILE__);
			strcpy(system->traj_output,token[1]);
		} else return 1;
	}
	else if (!strcasecmp(token[0], "traj_input")) {
		if(!system->traj_input) { 
			system->traj_input = calloc(MAXLINE,sizeof(char));
			memnullcheck(system->traj_input,MAXLINE*sizeof(char),__LINE__-1, __FILE__);
			strcpy(system->traj_input,token[1]);
		} else return 1;
	}
	else if (!strcasecmp(token[0], "energy_output")) {
		if(!system->energy_output) {
			system->energy_output = calloc(MAXLINE,sizeof(char));
			memnullcheck(system->energy_output,MAXLINE*sizeof(char),__LINE__-1, __FILE__);
			strcpy(system->energy_output,token[1]);
		} else return 1;
	}
	else if (!strcasecmp(token[0], "energy_output_csv")) {
		if(!system->energy_output_csv) {
			system->energy_output_csv = calloc(MAXLINE,sizeof(char));
			memnullcheck(system->energy_output_csv,MAXLINE*sizeof(char),__LINE__-1, __FILE__);
			strcpy(system->energy_output_csv,token[1]);
		} else return 1;
	}
	else if (!strcasecmp(token[0], "pop_histogram_output")) {
		if(!system->histogram_output) {
			system->histogram_output = calloc(MAXLINE,sizeof(char));
			memnullcheck(system->histogram_output,MAXLINE*sizeof(char),__LINE__-1, __FILE__);
			strcpy(system->histogram_output,token[1]);
		} else return 1;
	}
	else if (!strcasecmp(token[0], "dipole_output")) {
		if(!system->dipole_output) {
			system->dipole_output = calloc(MAXLINE,sizeof(char));
			memnullcheck(system->dipole_output,MAXLINE*sizeof(char),__LINE__-1, __FILE__);
			strcpy(system->dipole_output,token[1]);
		} else return 1;
	}
	else if (!strcasecmp(token[0], "field_output")) {
		if(!system->field_output) {
			system->field_output = calloc(MAXLINE,sizeof(char));
			memnullcheck(system->field_output,MAXLINE*sizeof(char),__LINE__-1, __FILE__);
			strcpy(system->field_output,token[1]);
		} else return 1;
	}
	else if (!strcasecmp(token[0], "frozen_output")) {
		if(!system->frozen_output) {
			system->frozen_output = calloc(MAXLINE,sizeof(char));
			memnullcheck(system->frozen_output,MAXLINE*sizeof(char),__LINE__-1, __FILE__);
			strcpy(system->frozen_output,token[1]);
		} else return 1;
	}
	else if (!strcasecmp(token[0], "insert_input")) {
		if(!system->insert_input) {
			system->insert_input = calloc(MAXLINE,sizeof(char));
			memnullcheck(system->insert_input,MAXLINE*sizeof(char),__LINE__-1, __FILE__);
			strcpy(system->insert_input,token[1]);
		} else return 1;
	}

	// Force long (non-PDB compliant, %11.6f) output of coordinates
	else if(!strcasecmp(token[0], "long_output")) {
		if(!strcasecmp(token[1],"on")) {
			system->long_output = 1;
			output("INPUT: Long coordinate output requested.\n");
		} else if (!strcasecmp(token[1],"off")) {
			system->long_output = 0;
		} else return 1;
	}

	// read box limits from pqr input
	else if(!strcasecmp(token[0], "read_pqr_box")) {
		if(!strcasecmp(token[1],"on"))
			system->read_pqr_box_on = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->read_pqr_box_on = 0;
		else return 1;
	}

	// surface fit input parameters
	
	else if( !strcasecmp( token[0], "fit_schedule")) {
		register double temp = atof(token[1]);
		if( temp<= 0.0 || temp>=1.0 ) {
			error( "INPUT: Invalid schedule.\n");
			return 1;
		}
		else
			system->fit_schedule = temp;
	}
	else if( !strcasecmp( token[0], "fit_max_energy")) {
		register double temp = atof(token[1]);
		if( temp<= 0.0 ) {
			error( "INPUT: fit_max_energy parameter must be greater than zero.\n");
			return 1;
		}
		else
			system->fit_max_energy = temp;
	}
	else if( !strcasecmp( token[0], "fit_start_temp" )) {
		register double temp = atof(token[1]);
		if( temp<= 0.0 ) {
			error( "INPUT: fit_start_temp parameter must be greater than zero.\n");
			return 1;
		}
		else
			{ if ( safe_atof(token[1],&(system->fit_start_temp)) ) return 1; }
	}
	else if(!strcasecmp(token[0], "fit_boltzmann_weight")) {
		if(!strcasecmp(token[1],"on"))
			system->fit_boltzmann_weight = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->fit_boltzmann_weight = 0;
		else return 1;
	}
	else if(!strcasecmp(token[0], "fit_input")) {
		// navigate to the end of the input file linked list
		fileNode_t *node = &(system->fit_input_list);
		while( node->next )
			node = node->next;
		// allocate a new node
		if( !(node->next = malloc( sizeof(fileNode_t)   ))   ) {
			error( "INPUT: Exhausted memory during input file node allocation.\n");
			return (1);
		}
		// advance to new node and initialize
		node = node->next;
		node->next = 0; // terminate list
		if(   !(node->data.filename = calloc(MAXLINE, sizeof(char)))   ) {
			error( "INPUT: Exhausted memory during string allocation for fit input filename.\n");
			return (1);
		}
		// copy filename to node and increment list count
		strcpy(node->data.filename, token[1]);
		system->fit_input_list.data.count++;
	}	

	// set basis
	else if(!strcasecmp(token[0], "basis1")) {
		{ if ( safe_atof(token[1],&(system->pbc->basis[0][0])) ) return 1; }
		{ if ( safe_atof(token[2],&(system->pbc->basis[0][1])) ) return 1; }
		{ if ( safe_atof(token[3],&(system->pbc->basis[0][2])) ) return 1; }
	}
	else if(!strcasecmp(token[0], "basis2")) {
		{ if ( safe_atof(token[1],&(system->pbc->basis[1][0])) ) return 1; }
		{ if ( safe_atof(token[2],&(system->pbc->basis[1][1])) ) return 1; }
		{ if ( safe_atof(token[3],&(system->pbc->basis[1][2])) ) return 1; }
	}
	else if(!strcasecmp(token[0], "basis3")) {
		{ if ( safe_atof(token[1],&(system->pbc->basis[2][0])) ) return 1; }
		{ if ( safe_atof(token[2],&(system->pbc->basis[2][1])) ) return 1; }
		{ if ( safe_atof(token[3],&(system->pbc->basis[2][2])) ) return 1; }
	}

	else if(!strcasecmp(token[0], "max_bondlength"))
		{ if ( safe_atof(token[1],&(system->max_bondlength)) ) return 1; }

	else if(!strcasecmp(token[0], "pop_histogram")) {
		if(!strcasecmp(token[1],"on"))
			system->calc_hist = 1;
		else if (!strcasecmp(token[1],"off")) 
			system->calc_hist = 0;
		else return 1;
	}

	else if (!strcasecmp(token[0], "pop_hist_resolution"))
		{ if ( safe_atof(token[1],&(system->hist_resolution)) ) return 1; }
	
	else return 1; //no match

	return 0;
}		



void setdefaults(system_t * system) {

	/* set the default scaling to 1 */
	system->scale_charge = 1.0;
	system->rot_probability = 1.0;
	system->spinflip_probability = 0.0;
	system->volume_probability = 0.0;

	/* set default volume change factor (for NPT) to 0.25 */
	system->volume_change_factor = 0.25;

	/* set histogram flag default: off */
	system->calc_hist=0;
	system->hist_resolution=0.0;
	system->histogram_output=NULL;

	/* default ewald parameters */
	system->ewald_alpha = EWALD_ALPHA;
	system->ewald_kmax = EWALD_KMAX;
	system->polar_ewald_alpha = EWALD_ALPHA;
	system->polar_wolf_alpha_lookup_cutoff = 30.0; //angstroms

	/* default polarization parameters */
	system->polar_gamma = 1.0;

	/* default rd LRC flag */
	system->rd_lrc = 1;

	// Initialize fit_input_list to reflect an empty list
	system->fit_input_list.next       = 0;
	system->fit_input_list.data.count = 0;

	// Initialize surface fitting parameters
	system->surf_fit_arbitrary_configs = 0;   // default is to use a small set of orientations, varying their separation
	system->fit_schedule               = 0;
	system->fit_start_temp             = 0;
	system->fit_max_energy             = 0;

	// Initialize insertion parameters
	system->num_insertion_molecules   = 0;
	system->insertion_molecules       = (molecule_t  *) 0;
	system->insertion_molecules_array = (molecule_t **) 0;

#ifdef QM_ROTATION
	/* default QR parameters */
	system->quantum_rotation_level_max = QUANTUM_ROTATION_LEVEL_MAX;
	system->quantum_rotation_l_max = QUANTUM_ROTATION_L_MAX;
	system->quantum_rotation_theta_max = QUANTUM_ROTATION_THETA_MAX;
	system->quantum_rotation_phi_max = QUANTUM_ROTATION_PHI_MAX;
	system->quantum_rotation_sum = QUANTUM_ROTATION_SUM;
#endif /* QM_ROTATION */

	//set default jobname
	system->job_name = calloc(MAXLINE,sizeof(char));
	memnullcheck(system->job_name,MAXLINE*sizeof(char),__LINE__-1, __FILE__);
	sprintf(system->job_name,"untitled");

	return;
}

system_t *read_config(char *input_file) {

	system_t *system;
	char linebuffer[MAXLINE], *n;
	char errormsg[MAXLINE];
	char ** token;
	FILE *fp;
	int i, linenum;

	system = calloc(1, sizeof(system_t));
	memnullcheck(system,sizeof(system_t),__LINE__-1, __FILE__);

	system->pbc = calloc(1, sizeof(pbc_t));
	memnullcheck(system->pbc, sizeof(pbc_t),__LINE__-1, __FILE__);

	/* open the config file or error */
	fp = fopen(input_file, "r");
	filecheck(fp, input_file, READ);

	/* allocate space for tokens */
	token = calloc(10, sizeof(char *));
	memnullcheck(token, 10*sizeof(char *), __LINE__-1, __FILE__);
	for ( i=0; i<10; i++ ) {
		token[i] = calloc(MAXLINE, sizeof(char));
		memnullcheck(token[i], MAXLINE*sizeof(char), __LINE__-1, __FILE__);
	}

	/* set default vaules */
	setdefaults(system);

	/* loop over each line */
	memset(linebuffer, 0, MAXLINE);
	n = fgets(linebuffer, MAXLINE, fp);
	
	linenum=0;
	while(n) {

		linenum++;
		/* grab a line and parse it out */
		for ( i=0; i<10; i++ )
			memset(token[i], 0, MAXLINE); //clear a token
			sscanf(linebuffer, "%s %s %s %s %s %s %s %s %s %s", 
			token[0], token[1], token[2], token[3], token[4], 
			token[5], token[6], token[7], token[8], token[9]);

		//parse and apply a command
		if ( do_command(system, token) != 0 ) {
			sprintf(errormsg,"INPUT: invalid command on line %d.\n", linenum);
			error(errormsg);
			sprintf(errormsg,"> %s\n", linebuffer);
			error(errormsg);
			return(NULL);
		}	

		memset(linebuffer, 0, MAXLINE);
		n = fgets(linebuffer, MAXLINE, fp);

	}

	/* close the config file */
	fclose(fp);

	for (i=0; i<10; i++)
		free(token[i]);
	free(token);

	return(system);

}


system_t *setup_system(char *input_file) {

	system_t *system;
	char linebuf[MAXLINE];
	FILE * finput;

	//read the main input file and sets values to flags and sets options
	system = read_config(input_file); 
	if(!system) {
		// error message generated in read_config()
		return(NULL);
	} else
		output("INPUT: finished reading config file\n");
	
	/* validate configuration parameters */
	if(check_system(system) < 0) {
		error("INPUT: invalid config parameters specified\n");
		return(NULL);
	} else
		output("INPUT: config file validated\n");

	/* set up the simulation box: pbc and read in molecules */
	if ( system->ensemble == ENSEMBLE_REPLAY ) {
		finput = fopen(system->traj_input,"r");
		filecheck(finput,system->traj_input,READ);
	}
	else {
		finput = fopen(system->pqr_input,"r");
		filecheck(finput,system->pqr_input,READ);
	}
	rewind(finput);
	if(setup_simulation_box(finput,system) <0) {
		error("INPUT: error setting up simulation box.\n");
		return(NULL);
	}	else
		output("INPUT: simulation box configured.\n");
	fclose(finput);

	/* allocate the necessary pairs */
	setup_pairs(system);
	output("INPUT: finished allocating pair lists\n");

	/* get all of the pairwise interactions, exclusions, etc. */
	if(system->cavity_bias) setup_cavity_grid(system);
	pairs(system);

	/* set all pairs to initially have their energies calculated */
	flag_all_pairs(system);
	output("INPUT: finished calculating pairwise interactions\n");

	if(!(system->sg || system->rd_only)) {
		sprintf(linebuf, "INPUT: Ewald gaussian width = %f A\n", system->ewald_alpha);
		output(linebuf);
		sprintf(linebuf, "INPUT: Ewald kmax = %d\n", system->ewald_kmax);
		output(linebuf);
	}

	/* ensure that all SPECTRE charges lie within the restricted domain */
	if(system->spectre) spectre_wrapall(system);

	return(system);

}






// readFitInputFiles() stores all of data points found in the fit input files
// linked to by system->fit_input_list. It returns an array of type curveData_t,
// each element of which provides specific data about the curve it describes, as
// well as an array (a *double, actually) containing the actual point data.

curveData_t *readFitInputFiles( system_t *system, int nCurves )
{
	// General Algorithm:
	// This function traverses the list of fit_input files (@ system->fit_input_list), reading
	// each one a line at a time. It parses the line, determining whether it is a line with
	// data or metadata. If metadata, the value is stored in the appropriate curveData_t
	// field. If it is a point value, it is temporarily stored in a linked list--this is
	// done because A) the number of data points is not initially known, and B) the list can
	// be constructed in sorted order--just in case the points are stored out of order
	// in the input file. Once the entire file has been processed, an array is allocated and
	// the list data is transferred into the array (the function that processes the data
	// expects an array). The function then moves on to the next input file...
	// Before returning, the function verifies that each curve has the same number of points
	// and that the r-values of the curves all correspond to each other, and that the increment
	// between r-values is identical (to within threshold_for_identical_DeltaR). If everything
	// checks out, the normalized curve-weights are calculated and the function returns.
	
	// Precision to use when comparing delta-r values for equality:
	const double threshold_for_identical_DeltaR = 1.0E-12;
	
	
	FILE *fp_fit;
	char linebuf[MAXLINE ],  // Line buffer for reading file input
	     token1[ MAXLINE ],  // Variables for parsing parameters in input files.
	     token2[ MAXLINE ],
	     token3[ MAXLINE ],
	     errMsg[ MAXLINE ];  // Output buffer for error messages
	
	
	
	// Allocate memory for list of curve data/metadata
	curveData_t *curve_input = calloc( nCurves, sizeof(curveData_t));
	if( !curve_input )
	{
		error( "INPUT: Exhausted memory while allocating array for curve data.\n" );
		return (NULL);
	}





	// Read in point data/metadata array for each input curve and store
	// the data in one of the curveData_t array elements we just allocated
	///////////////////////////////////////////////////////////////////////
	
	fileNode_t *fNode = system->fit_input_list.next;
	int currentCurve = 0;
	
	while( fNode )
	{

		// orderedPairNode_t's will be used in making linked lists of data points.
		// When all points have been read, an array of the appropriate size will
		// be allocated and the list items will be transferred into the array.
		// The linked list will be maintained in order, sorted by r-value from
		// least to greatest. The list is never explicitly sorted, but when each
		// element is added, it is inserted at its in-order location.
		typedef struct _orderedPairNode {
			double r; // Separation distance
			double E; // Energy value
			struct _orderedPairNode *next;
		} orderedPairNode_t;

		// The list head, although a node, contains no data itself
		// and only serves as a pointer to the first element.
		orderedPairNode_t pointList; // The list head
		pointList.next = 0;          // The list starts out empty
		
		orderedPairNode_t *nodePtr = &pointList; // Our list iterator
		
		
		
		// Initialize metadata for the current curve
		curve_input[ currentCurve ].nPoints = 0;
		curve_input[ currentCurve ].alpha1  = 0;
		curve_input[ currentCurve ].beta1   = 0;
		curve_input[ currentCurve ].gamma1  = 0;
		curve_input[ currentCurve ].alpha2  = 0;
		curve_input[ currentCurve ].beta2   = 0;
		curve_input[ currentCurve ].gamma2  = 0;
		curve_input[ currentCurve ].r       = 0;
		curve_input[ currentCurve ].input   = 0;
		curve_input[ currentCurve ].weight  = 1;
		
		
		
		
		
		// Open file for reading
		////////////////////////////
		
		char *filename = fNode->data.filename;
		printf( "INPUT: Loading %s\n", filename );
		fp_fit = fopen( filename, "r" );
				filecheck(fp_fit,filename,READ);



		// Strip path from filename and save the remainder
		// for use as an id in the curveData_t structure
		////////////////////////////////////////////////////
		{

			char *i  = filename; // iterator for raw filename
			char *fn = filename; // temp pointer saving last good place in raw filename
			
			// Define the appropriate file separator character for the OS
			#ifndef _WIN32
			const char file_separator = '/';
			#else
			const char file_separator = '\\';
			#endif
			
			// Use iterator to look for file separator characters in the filename
			while( *i != '\0' )
			{
				if( (*i) == file_separator )
					fn = i+1;
				i++;
			}
			
			curve_input[currentCurve].filename = malloc( (strlen(fn) + 1 ) * sizeof(char) );
			// Use the filename as the id, at least until a proper id is found in the file
			curve_input[currentCurve].id       = malloc( (strlen(fn) + 1 ) * sizeof(char) );
			if( !(curve_input[currentCurve].filename && curve_input[currentCurve].filename) )
			{
				error( "INPUT: Exhausted memory while allocating curve id/filename fields.\n" );
				if( curve_input[currentCurve].filename )
					free( curve_input[currentCurve].filename );
				if( curve_input[currentCurve].id )
					free( curve_input[currentCurve].id );
				return (NULL);
			}
			strcpy( curve_input[currentCurve].filename, fn );
			strcpy( curve_input[currentCurve].id,       fn );
		}





		// Parse and store data found in the input file
		//////////////////////////////////////////////////

		while( fgets(linebuf, MAXLINE, fp_fit) )
		{
			memset(token1, 0, MAXLINE);
			memset(token2, 0, MAXLINE);
			memset(token3, 0, MAXLINE);
			sscanf(linebuf, "%s %s %s", token1, token2, token3 );
			
			if( !strcasecmp(token1, "id")) {
				// Free memory for previous id
				free( curve_input[currentCurve].id );
				// Allocate memory for new id
				curve_input[currentCurve].id = malloc( (strlen(token2) + 1) * sizeof(char) );
				// If allocation fails, free all memory and exit
				if( !curve_input[currentCurve].id ) {
					sprintf( errMsg, "INPUT: Exhausted memory while reading data from %s\n", filename );
					error( errMsg );
					fclose( fp_fit );
					// Free list memory
					nodePtr = pointList.next;
					while( nodePtr ){orderedPairNode_t *temp = nodePtr; nodePtr = nodePtr->next; free(temp);}
					// Free id memory
					int i;
					for( i=0; i<=currentCurve; i++ ) {
						if( curve_input[i].filename )
							free( curve_input[i].filename );
						if( curve_input[i].id )
							free( curve_input[i].id);
					}
					return( NULL );
				}
				// If allocation was successful, copy new id into curve_input
				strcpy( curve_input[ currentCurve ].id, token2 );
			}

			else if( !strcasecmp(token1, "alpha1" )) {
				curve_input[ currentCurve ].alpha1 = atof( token2 );
				if( !strcasecmp( token3, "pi" ))
					curve_input[ currentCurve ].alpha1 *= M_PI;
			}

			else if( !strcasecmp(token1, "alpha2" )) {
				curve_input[ currentCurve ].alpha2 = atof( token2 );
				if( !strcasecmp( token3, "pi" ))
					curve_input[ currentCurve ].alpha2 *= M_PI;
			}

			else if( !strcasecmp(token1, "beta1" )) {
				curve_input[ currentCurve ].beta1 = atof( token2 );
				if( !strcasecmp( token3, "pi" ))
					curve_input[ currentCurve ].beta1 *= M_PI;
			}
			
			else if( !strcasecmp(token1, "beta2" )) {
				curve_input[ currentCurve ].beta2 = atof( token2 );
				if( !strcasecmp( token3, "pi" ))
					curve_input[ currentCurve ].beta2 *= M_PI;
			}
			
			else if( !strcasecmp(token1, "gamma1" )) {
				curve_input[ currentCurve ].gamma1 = atof( token2 );
				if( !strcasecmp( token3, "pi" ))
					curve_input[ currentCurve ].gamma1 *= M_PI;
			}
			
			else if( !strcasecmp(token1, "gamma2" ))  {
				curve_input[ currentCurve ].gamma2 = atof( token2 );
				if( !strcasecmp( token3, "pi" ))
					curve_input[ currentCurve ].gamma2 *= M_PI;
			}
			
			else if( !strcasecmp(token1, "weight" ))
				curve_input[ currentCurve ].weight = atof( token2 );
			
			else if(  (token1[0]=='*') || (token1[0]==0)  )
				; // Do nothing for comment lines and blank lines
			else
			{

				// If we get here, this line SHOULD be a data-point.
				// We will now store the point in our linked list.
				//////////////////////////////////////////////////////


				// Create a new node
				orderedPairNode_t *newNodePtr = malloc( sizeof(orderedPairNode_t));
				if( !newNodePtr )
				{
					sprintf( errMsg, "INPUT: Exhausted memory while reading data from %s\n", filename );
					error( errMsg );
					fclose( fp_fit );
					// Free list memory
					nodePtr = pointList.next;
					while( nodePtr ){orderedPairNode_t *temp = nodePtr; nodePtr = nodePtr->next; free(temp);}
					// Free id memory
					int i; for( i=0; i<=currentCurve; i++ ) { free( curve_input[i].filename ); free( curve_input[i].id);}
					return (NULL);
				}

				// Transfer data to new node
				double r =                      // for convenient referencing of current r-value
				newNodePtr->r = atof( token1 ); // molecule-molecule distance
				newNodePtr->E = atof( token2 ); // energy value at distance r
				
				// Check for an invalid r-value
				// An r-value of 0.0 may be indicative of a corrupt input line or an invalid keyword
				if( r <= 0.0 )
				{
					sprintf( errMsg, "INPUT: Invalid r-value (%.2lf) in %s\n", r, filename );
					error( errMsg );
					if( r==0.0 )
					{
						error( "       Verify that all keywords in input file are spelled correctly.\n" );
						error( "       Keywords are not case sensitive.\n" );
					}
					fclose( fp_fit );
					// Free list memory
					nodePtr = pointList.next;
					while( nodePtr ){orderedPairNode_t *temp = nodePtr; nodePtr = nodePtr->next; free(temp);}
					// Free id memory
					int i; for( i=0; i<=currentCurve; i++ ) { free( curve_input[i].filename ); free( curve_input[i].id);}
					return (NULL);
				}

				// Increment point-count metadata for this curve
				curve_input[ currentCurve ].nPoints++;
				
				// Reset nodePtr and advance it to the correct list location
				nodePtr = &pointList;
				while(   (nodePtr->next)   &&   ((nodePtr->next->r) < r)   ) {
					// Stop list-traversal if there is no next item, or if next r-value is >= to the current
					// r-value (consequently maintaining a sorted list).
					// Recall that for the && operator, the 2nd condition will not be evaluated if the 1st
					// fails, i.e., we will not check the next r-value if there is no next r-value...
					nodePtr = nodePtr->next;
				}
				
				// Link newNode into the list at the location determined by previous while loop
				newNodePtr->next = nodePtr->next;
				nodePtr->next = newNodePtr;


			} // done processing this line, proceed to the next...
		} // while (fgets), i.e.  while(!EOF)
		fclose(fp_fit);






		// Ensure some point data was actually read
		////////////////////////////////////////////
		if(!curve_input[ currentCurve ].nPoints)
		{
			sprintf( errMsg, "INPUT: No data points found in file %s\n", filename );
			error( errMsg );
			// Free list memory
			nodePtr = pointList.next;
			while( nodePtr ){orderedPairNode_t *temp = nodePtr; nodePtr = nodePtr->next; free(temp);}
			// Free memory for curve identification
			int i; for( i=0; i<=currentCurve; i++ ) {free( curve_input[i].filename ); free( curve_input[i].id);}
			return(NULL);
		}
		
		
		
		
		
		
		// Allocate array space to hold data points and copy
		// data from the list to the array.
		///////////////////////////////////////////////////////
		
		curve_input[ currentCurve ].r      = calloc( curve_input[ currentCurve ].nPoints, sizeof(double) );
		curve_input[ currentCurve ].input  = calloc( curve_input[ currentCurve ].nPoints, sizeof(double) );
		
		if( !curve_input[ currentCurve ].r || !curve_input[ currentCurve ].input )
		{
			sprintf( errMsg, "INPUT: Exhausted memory transferring curve data to array, for curve: %s\n", filename );
			error( errMsg );
			// Free list memory
			nodePtr = pointList.next;
			while( nodePtr ){orderedPairNode_t *temp = nodePtr; nodePtr = nodePtr->next; free(temp);}
			// Free curveData_t memory
			int i;
			for( i=0; i<=currentCurve; i++ ) {
				if( curve_input[i].filename ) free( curve_input[i].filename );
				if( curve_input[i].id       ) free( curve_input[i].id );
				if( curve_input[i].r        ) free( curve_input[i].r );
				if( curve_input[i].input    ) free( curve_input[i].input    );
			}
			return (NULL);
		}

		nodePtr = &pointList;  // reset nodePtr to list head
		int N=0;               // reset N to first array element
		// Copy point data from the linked list into the array
		while( nodePtr->next )
		{
			nodePtr = nodePtr->next;
			curve_input[ currentCurve ].r[N]     = nodePtr->r;
			curve_input[ currentCurve ].input[N] = nodePtr->E;
			++N;
		}

		// Free memory used by the list
		nodePtr = pointList.next;
		while( nodePtr )
		{
			orderedPairNode_t *temp = nodePtr;
			nodePtr = nodePtr->next;
			free(temp);
		}



		fNode = fNode->next;
		++currentCurve;
	} // while(fNode) i.e., Process next file in the list...






	 // Normalize curve weights 
	////////////////////////////
	
	int totalWeight = 0;
	// Get the combined sum of all curve weights
	for( currentCurve=0; currentCurve < nCurves; currentCurve++ )
		totalWeight += curve_input[currentCurve].weight;
	// Divide each curves individual weight by the total weight
	for( currentCurve=0; currentCurve < nCurves; currentCurve++ )
		curve_input[currentCurve].normalized_weight = (double) curve_input[currentCurve].weight / totalWeight;







	// Check data validity...
	///////////////////////////
	
	// Check that all curves have the same number of points
	int lastPointCount = curve_input[0].nPoints;
	for( currentCurve=1; currentCurve<nCurves; ++currentCurve )
	if( lastPointCount == curve_input[currentCurve].nPoints )
		lastPointCount =  curve_input[currentCurve].nPoints;
	else
	{
		error( "INPUT: All curves must have the same number of data points.\n" );
		// Free curveData_t memory
		int i;
		for( i=0; i<=currentCurve; i++ ) {
			if( curve_input[i].filename ) free( curve_input[i].filename );
			if( curve_input[i].id       ) free( curve_input[i].id );
			if( curve_input[i].r        ) free( curve_input[i].r );
			if( curve_input[i].input    ) free( curve_input[i].input    );
		}
		return( NULL );
	}


	// Check that all corresponding points on all the curves have identical r-values
	int currentPoint;
	for( currentPoint=0; currentPoint<curve_input[0].nPoints; ++currentPoint )
	{
		double lastRValue = curve_input[0].r[currentPoint];
		for( currentCurve=1; currentCurve<nCurves; ++currentCurve ) {
			if( lastRValue == curve_input[currentCurve].r[currentPoint] )
				lastRValue =  curve_input[currentCurve].r[currentPoint];
			else
			{
				error( "INPUT: Every curve must have identical r-values.\n" );
				// Free curveData_t memory
				int i;
				for( i=0; i<=currentCurve; ++i ) {
					if( curve_input[i].filename ) free( curve_input[i].filename );
					if( curve_input[i].id       ) free( curve_input[i].id );
					if( curve_input[i].r        ) free( curve_input[i].r );
					if( curve_input[i].input    ) free( curve_input[i].input    );
				}
				return( NULL );
			}
		}
	}


	// Ensure that delta-r values are uniform
	// Since we know the magnitude and quantity of r-values are identical across all
	// curves, we will only check the delta-r's of the first curve.
	
	
	// The standard by which all other delta-r values will be judged:
	double reference_deltaR = curve_input[0].r[1] - curve_input[0].r[0];
	
	for( currentPoint=2; currentPoint<curve_input[0].nPoints; ++currentPoint )
	{
		double deltaR = curve_input[0].r[currentPoint] - curve_input[0].r[currentPoint-1];
		if( fabs(deltaR - reference_deltaR) > threshold_for_identical_DeltaR )
		{
			error( "INPUT: Data points on curve must be evenly spaced.\n         " );
			sprintf( errMsg, "       Tolerance currently set at: +/- %.17lf\n", threshold_for_identical_DeltaR );
			error( errMsg );
			// Free curveData_t memory
			int i;
			for( i=0; i<=currentCurve; i++ ) {
				if( curve_input[i].filename ) free( curve_input[i].filename );
				if( curve_input[i].filename ) free( curve_input[i].id       );
				if( curve_input[i].r        ) free( curve_input[i].r        );
				if( curve_input[i].input    ) free( curve_input[i].input    );
			}
			return( NULL );
		}
	}

	return curve_input;
}



