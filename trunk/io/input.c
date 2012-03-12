/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>
#define au2invseconds 4.13412763705666648752113572754445220741745180640e16 

system_t *read_config(char *input_file) {

	system_t *system;
	char linebuffer[MAXLINE], *n;
	char token1[MAXLINE], token2[MAXLINE], token3[MAXLINE], token4[MAXLINE], token5[MAXLINE], token6[MAXLINE], token7[MAXLINE];
	FILE *fp;

	system = calloc(1, sizeof(system_t));
	if(!system) {
		error("INPUT: couldn't allocate system data structure\n");
		return(NULL);
	}

	system->pbc = calloc(1, sizeof(pbc_t));
	if(!system->pbc) {
		error("INPUT: couldn't allocate the PBC data structure\n");
		return(NULL);
	}

	/* open the config file or error */
	fp = fopen(input_file, "r");
	if(!fp) {
		error("INPUT: couldn't open config file\n");
		return(NULL);
	}

	/* set the default scaling to 1 */
	system->scale_charge = 1.0;
	system->scale_rd = 1.0;
	system->rot_probability = 1.0;
	system->spinflip_probability = 0.0;

	/* set histogram flag default: off */
	system->calc_hist=0;
	system->hist_resolution=0.0;
	system->histogram_output=NULL;

	/* default ewald parameters */
	system->ewald_alpha = EWALD_ALPHA;
	system->ewald_kmax = EWALD_KMAX;

	/* default polarization parameters */
	system->polar_gamma = 1.0;

	/* default rd LRC flag */
	system->rd_lrc = 1;

	/*default vdw feynman-hibbs correction cap*/
	system->VDW_FH_cap = -1000;

        // Initialize fit_input_list to reflect an empty list
        system->fit_input_list.next       = 0;
        system->fit_input_list.data.count = 0;

        // Initialize surface fitting parameters
        system->fit_schedule         = 0;
        system->fit_start_temp       = 0;
        system->fit_max_energy       = 0;

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

	/* loop over each line */
	memset(linebuffer, 0, MAXLINE);
	n = fgets(linebuffer, MAXLINE, fp);
	while(n) {

		/* grab a line and parse it out */
		memset(token1, 0, MAXLINE);
		memset(token2, 0, MAXLINE);
		memset(token3, 0, MAXLINE);
		memset(token4, 0, MAXLINE);
		memset(token5, 0, MAXLINE);
		memset(token6, 0, MAXLINE);
		memset(token7, 0, MAXLINE);
		sscanf(linebuffer, "%s %s %s %s %s %s %s", token1, token2, token3, token4, token5, token6, token7);

		if(!strcasecmp(token1, "ensemble")) {
			if(!strcasecmp(token2, "nvt"))
				system->ensemble = ENSEMBLE_NVT;
			else if(!strcasecmp(token2, "uvt"))
				system->ensemble = ENSEMBLE_UVT;
			else if(!strcasecmp(token2, "surf"))
				system->ensemble = ENSEMBLE_SURF;
			else if(!strcasecmp(token2, "surf_fit"))
				system->ensemble = ENSEMBLE_SURF_FIT;
			else if(!strcasecmp(token2, "nve"))
				system->ensemble = ENSEMBLE_NVE;
			else if(!strcasecmp(token2, "total_energy")) {
				system->ensemble = ENSEMBLE_TE;
				system->numsteps = 0; //do no mc steps
			}

		}

		/* for NVE only */
		if(!strcasecmp(token1, "total_energy"))
                        system->total_energy = atof(token2);

		if(!strcasecmp(token1, "surf_decomp")) {
			if(!strcasecmp(token2, "on"))
				system->surf_decomp = 1;
			else
				system->surf_decomp = 0;
		}

		if(!strcasecmp(token1, "surf_min"))
			system->surf_min = atof(token2);

		if(!strcasecmp(token1, "surf_max"))
			system->surf_max = atof(token2);

		if(!strcasecmp(token1, "surf_inc"))
			system->surf_inc = atof(token2);

		if(!strcasecmp(token1, "surf_ang"))
			system->surf_ang = atof(token2);

		//allows us to specify the surf-fit scales in the input file
		if(!strcasecmp(token1, "surf_weight_constant")) {
			system->surf_weight_constant = atof(token2);
			system->surf_weight_constant_on = 1;
		}

		if(!strcasecmp(token1, "surf_scale_q")) {
			system->surf_scale_q = atof(token2);
			system->surf_scale_q_on = 1;
		}

		if(!strcasecmp(token1, "surf_scale_r")) {
			system->surf_scale_r = atof(token2);
			system->surf_scale_r_on = 1;
		}

		if(!strcasecmp(token1, "surf_scale_epsilon")) {
			system->surf_scale_epsilon = atof(token2);
			system->surf_scale_epsilon_on = 1;
		}

		if(!strcasecmp(token1, "surf_scale_sigma")) {
			system->surf_scale_sigma = atof(token2);
			system->surf_scale_sigma_on = 1;
		}

		if(!strcasecmp(token1, "surf_scale_omega")) {
			system->surf_scale_omega = atof(token2);
			system->surf_scale_omega_on = 1;
		}

		//let the charge site shift during fitting
		if(!strcasecmp(token1, "surf_qshift")) {
			if (!strcasecmp(token2, "on")) {
				system->surf_qshift_on = 1;
				fprintf(stdout,"INPUT: SURF_QSHIFT IS ON. "
					"CHARGED SITES ARE ALLOWED TO MIGRATE.\n"
					"INPUT: THE QUADRUPOLE MOMENT WILL BE PRESERVED.\n"
					"INPUT: *** ONLY USE WITH LINEAR MOLECULES ALIGNED ON THE X-AXIS. ***\n"
					"INPUT: Label the charged sites, H2G and H2Q.\n");
			}
			else {
				fprintf(stderr,"INPUT: ERROR: Invalid argument to qshift option.\n");
				return(NULL);
			}
		}

		if(!strcasecmp(token1, "surf_preserve")) {
			if(!strcasecmp(token2, "on"))
				system->surf_preserve = 1;
			else
				system->surf_preserve = 0;
		}

	//set rotation in input file for surf_preserve (calc dimer geometry curves)
		if(!strcasecmp(token1, "surf_preserve_rotation")) {
			if(system->surf_preserve_rotation_on != NULL) {
				fprintf(stderr,"ERROR: surf_preserve_rotation already set.\n");
				return(NULL);
			}
			system->surf_preserve_rotation_on = malloc(6*sizeof(double));
			system->surf_preserve_rotation_on->alpha1=atof(token2);
			system->surf_preserve_rotation_on->beta1=atof(token3);
			system->surf_preserve_rotation_on->gamma1=atof(token4);
			system->surf_preserve_rotation_on->alpha2=atof(token5);
			system->surf_preserve_rotation_on->beta2=atof(token6);
			system->surf_preserve_rotation_on->gamma2=atof(token7);
		}

		if(!strcasecmp(token1, "spectre")) {
			if(!strcasecmp(token2, "on"))
				system->spectre = 1;
			else
				system->spectre = 0;
		}

		if(!strcasecmp(token1, "spectre_max_charge"))
			system->spectre_max_charge = fabs(atof(token2));

		if(!strcasecmp(token1, "spectre_max_target"))
			system->spectre_max_target = fabs(atof(token2));

		if(!strcasecmp(token1, "seed"))
			system->seed = atol(token2);

		if(!strcasecmp(token1, "numsteps"))
			system->numsteps = atoi(token2);

		if(!strcasecmp(token1, "corrtime"))
			system->corrtime = atoi(token2);

		if(!strcasecmp(token1, "move_probability"))
			system->move_probability = atof(token2);

		if(!strcasecmp(token1, "gwp_probability"))
			system->gwp_probability = atof(token2);

		if(!strcasecmp(token1, "rot_probability"))
			system->rot_probability = atof(token2);

		if(!strcasecmp(token1, "insert_probability"))
			system->insert_probability = atof(token2);

		if(!strcasecmp(token1, "adiabatic_probability"))
			system->adiabatic_probability = atof(token2);

		if(!strcasecmp(token1, "spinflip_probability"))
			system->spinflip_probability = atof(token2);

		if(!strcasecmp(token1, "cavity_bias")) {
			if(!strcasecmp(token2, "on"))
				system->cavity_bias = 1;
			else
				system->cavity_bias = 0;
		}

		if(!strcasecmp(token1, "cavity_grid"))
			system->cavity_grid_size = atoi(token2);

		if(!strcasecmp(token1, "cavity_radius"))
			system->cavity_radius = atof(token2);

		if(!strcasecmp(token1, "temperature"))
			system->temperature = atof(token2);

		if(!strcasecmp(token1, "simulated_annealing")) {
			if(!strcasecmp(token2, "on"))
				system->simulated_annealing = 1;
			else
				system->simulated_annealing = 0;
		}

		if(!strcasecmp(token1, "simulated_annealing_schedule"))
			system->simulated_annealing_schedule = atof(token2);

		if(!strcasecmp(token1, "pressure"))
			system->pressure = atof(token2);

		if(!strcasecmp(token1, "h2_fugacity")) {
			if(!strcasecmp(token2, "on"))
				system->h2_fugacity = 1;
			else
				system->h2_fugacity = 0;
		}

		if(!strcasecmp(token1, "co2_fugacity")) {
			if(!strcasecmp(token2, "on"))
				system->co2_fugacity = 1;
			else
				system->co2_fugacity = 0;
		}
		if(!strcasecmp(token1, "co2fug"))
		        system->fugacity=atof(token2);

                if(!strcasecmp(token1, "ch4_fugacity")) {
                        if(!strcasecmp(token2, "on"))
                                system->ch4_fugacity = 1; 
                        else
                                system->ch4_fugacity = 0; 
                }

		if(!strcasecmp(token1, "n2_fugacity")) {
			if(!strcasecmp(token2, "on"))
				system->n2_fugacity = 1;
			else
				system->n2_fugacity = 0;
		}

		if(!strcasecmp(token1, "free_volume"))
			system->free_volume = atof(token2);

		if(!strcasecmp(token1, "rd_only")) {
			if(!strcasecmp(token2, "on"))
				system->rd_only = 1;
			else
				system->rd_only = 0;
		}

		if(!strcasecmp(token1, "gwp")) {
			if(!strcasecmp(token2, "on"))
				system->gwp = 1;
			else
				system->gwp = 0;
		}

		if(!strcasecmp(token1, "wolf")) {
			if(!strcasecmp(token2, "on"))
				system->wolf = 1;
			else
				system->wolf = 0;
		}

		if(!strcasecmp(token1, "rd_lrc")) {
			if(!strcasecmp(token2, "on"))
				system->rd_lrc = 1;
			else
				system->rd_lrc = 0;
		}

		if(!strcasecmp(token1, "rd_anharmonic")) {
			if(!strcasecmp(token2, "on"))
				system->rd_anharmonic = 1;
			else
				system->rd_anharmonic = 0;
		}

		if(!strcasecmp(token1, "rd_anharmonic_k"))
			system->rd_anharmonic_k = atof(token2);

		if(!strcasecmp(token1, "rd_anharmonic_g"))
			system->rd_anharmonic_g = atof(token2);

		if(!strcasecmp(token1, "feynman_hibbs")) {
			if(!strcasecmp(token2, "on"))
				system->feynman_hibbs = 1;
			else
				system->feynman_hibbs = 0;
		}

		if(!strcasecmp(token1, "feynman_kleinert")) {
			if(!strcasecmp(token2, "on"))
				system->feynman_kleinert = 1;
			else
				system->feynman_kleinert = 0;
		}

		if(!strcasecmp(token1, "feynman_hibbs_order")) {
			system->feynman_hibbs_order = atoi(token2);
		}

		if(!strcasecmp(token1, "wpi")) {
			if(!strcasecmp(token2, "on"))
				system->wpi = 1;
			else
				system->wpi = 0;
		}

		if(!strcasecmp(token1, "wpi_grid"))
			system->wpi_grid = atoi(token2);

		if(!strcasecmp(token1, "fvm")) {
			if(!strcasecmp(token2, "on"))
				system->fvm = 1;
			else
				system->wpi = 0;
		}

		if(!strcasecmp(token1, "sg")) {
			if(!strcasecmp(token2, "on"))
				system->sg = 1;
			else if(!strcasecmp(token2, "off"))
				system->sg = 0;
		}

		if(!strcasecmp(token1, "dreiding")) {
			if(!strcasecmp(token2, "on"))
				system->dreiding = 1;
			else if(!strcasecmp(token2, "off"))
				system->dreiding = 0;
		}
		if(!strcasecmp(token1, "wrapall")) {
			if(!strcasecmp(token2, "on"))
				system->wrapall = 1;
			else
				system->wrapall = 0;
		}

		if(!strcasecmp(token1, "scale_charge")) {
			system->scale_charge = atof(token2);
		}

		if(!strcasecmp(token1, "scale_rd")) {
			system->scale_rd = atof(token2);
		}

		if(!strcasecmp(token1, "ewald_alpha")) {
			system->ewald_alpha = atof(token2);
		}

		if(!strcasecmp(token1, "ewald_kmax")) {
			system->ewald_kmax = atoi(token2);
		}

		if(!strcasecmp(token1, "pbc_cutoff")) {
			system->pbc->cutoff = atof(token2);
		}

 		//cap the value of FH correction to VDW
 		if (!strcasecmp(token1,"vdw_fh_cap")) {
 			system->VDW_FH_cap = atof(token2);
 			if ( system->VDW_FH_cap > 0 ) {
 				fprintf(stderr,"INPUT: Error: VDW_FH_cap must be negative.\n");
 				return(NULL);
 			}
 		}

		if(!strcasecmp(token1, "polarization")) {
			if(!strcasecmp(token2, "on"))
				system->polarization = 1;
			else
				system->polarization = 0;
		}

		if(!strcasecmp(token1, "polarvdw")) {
			if(!strcasecmp(token2, "on")) {
				system->polarvdw = 1;
				//polarvdw requires polarization matricies to be built
				system->polarization=1;
				system->polar_iterative=1; //matrix inversion destroys A_matrix before vdw can use it.
				fprintf(stdout,"INPUT: Forcing polar_iterative ON for CP-VdW.\n");
			}
			else if (!strcasecmp(token2, "evects")) {
				system->polarvdw = 2; //calculate eigenvectors
				system->polarization=1;
				system->polar_iterative=1; //matrix inversion destroys A_matrix before vdw can use it.
				fprintf(stdout,"INPUT: Forcing polar_iterative ON for CP-VdW.\n");
			}
			else if ( !strcasecmp(token2, "comp")) {
				system->polarvdw = 3; //calculate eigenvectors
				system->polarization=1;
				system->polar_iterative=1; //matrix inversion destroys A_matrix before vdw can use it.
				fprintf(stdout,"INPUT: Forcing polar_iterative ON for CP-VdW.\n");
			}
			else 
				system->polarvdw = 0;
		}

		if(!strcasecmp(token1, "cavity_autoreject")) {
			if(!strcasecmp(token2, "on"))
				system->cavity_autoreject = 1;
			else
				system->cavity_autoreject = 0;
		}

		if(!strcasecmp(token1, "cavity_autoreject_scale")) {
			system->cavity_autoreject_scale = atof(token2);
		}

		if(!strcasecmp(token1, "polar_ewald")) {
			if(!strcasecmp(token2, "on"))
				system->polar_ewald = 1;
			else
				system->polar_ewald = 0;
		}

		if(!strcasecmp(token1, "polarizability_tensor")) {
			if(!strcasecmp(token2, "on"))
				system->polarizability_tensor = 1;
			else
				system->polarizability_tensor = 0;
		}

		if(!strcasecmp(token1, "polar_zodid")) {
			if(!strcasecmp(token2, "on"))
				system->polar_zodid = 1;
			else
				system->polar_zodid = 0;
		}

		if(!strcasecmp(token1, "polar_iterative")) {
			if(!strcasecmp(token2, "on"))
				system->polar_iterative = 1;
			else
				system->polar_iterative = 0;
		}

		if(!strcasecmp(token1, "polar_palmo")) {
			if(!strcasecmp(token2, "on"))
				system->polar_palmo = 1;
			else
				system->polar_palmo = 0;
		}

		if(!strcasecmp(token1, "polar_gs")) {
			if(!strcasecmp(token2, "on"))
				system->polar_gs = 1;
			else
				system->polar_gs = 0;
		}

		if(!strcasecmp(token1, "polar_gs_ranked")) {
			if(!strcasecmp(token2, "on"))
				system->polar_gs_ranked = 1;
			else
				system->polar_gs_ranked = 0;
		}

		if(!strcasecmp(token1, "polar_sor")) {
			if(!strcasecmp(token2, "on"))
				system->polar_sor = 1;
			else
				system->polar_sor = 0;
		}

		if(!strcasecmp(token1, "polar_esor")) {
			if(!strcasecmp(token2, "on"))
				system->polar_esor = 1;
			else
				system->polar_esor = 0;
		}

		if(!strcasecmp(token1, "polar_gamma")) {
			system->polar_gamma = atof(token2);
		}

		if(!strcasecmp(token1, "polar_damp")) {
			system->polar_damp = atof(token2);
		}

		if(!strcasecmp(token1, "field_damp")) {
			system->field_damp = atof(token2);
		}

		if(!strcasecmp(token1, "polar_precision")) {
			system->polar_precision = atof(token2);
		}

		if(!strcasecmp(token1, "polar_max_iter")) {
			system->polar_max_iter = atoi(token2);
		}

		if(!strcasecmp(token1, "polar_damp_type")) {
			if(!strcasecmp(token2, "linear"))
				system->damp_type = DAMPING_LINEAR;
			else if(!strcasecmp(token2, "exponential"))
				system->damp_type = DAMPING_EXPONENTIAL;
			else
				system->damp_type = 0;
		}

		if(!strcasecmp(token1, "polar_self")) {
			if(!strcasecmp(token2, "on"))
				system->polar_self = 1;
			else
				system->polar_self = 0;
		}

		if(!strcasecmp(token1, "cuda")) {
			if(!strcasecmp(token2, "on"))
				system->cuda = 1;
			else
				system->cuda = 0;
		}

		if(!strcasecmp(token1, "independent_particle")) {
			if(!strcasecmp(token2, "on"))
				system->independent_particle = 1;
			else
				system->independent_particle = 0;
		}

#ifdef QM_ROTATION
		if(!strcasecmp(token1, "quantum_rotation")) {
			if(!strcasecmp(token2, "on"))
				system->quantum_rotation = 1;
			else
				system->quantum_rotation = 0;
		}

		if(!strcasecmp(token1, "quantum_rotation_hindered")) {
			if(!strcasecmp(token2, "on"))
				system->quantum_rotation_hindered = 1;
			else
				system->quantum_rotation_hindered = 0;
		}

		if(!strcasecmp(token1, "quantum_rotation_hindered_barrier")) {
			system->quantum_rotation_hindered_barrier = atof(token2);
		}

		if(!strcasecmp(token1, "quantum_rotation_B")) {
			system->quantum_rotation_B = atof(token2);
		}

		if(!strcasecmp(token1, "quantum_rotation_level_max")) {
			system->quantum_rotation_level_max = atof(token2);
		}

		if(!strcasecmp(token1, "quantum_rotation_l_max")) {
			system->quantum_rotation_l_max = atof(token2);
		}

		if(!strcasecmp(token1, "quantum_rotation_sum")) {
			system->quantum_rotation_sum = atof(token2);
		}

#endif /* QM_ROTATION */

#ifdef XXX
		if(!strcasecmp(token1, "quantum_vibration")) {
			if(!strcasecmp(token2, "on"))
				system->quantum_vibration = 1;
			else
				system->quantum_vibration = 0;
		}
#endif /* XXX */


		if(!strcasecmp(token1, "pdb_input")) {
			if(!system->pdb_input) {
				system->pdb_input = calloc(MAXLINE, sizeof(char));
				memnullcheck(system->pdb_input,MAXLINE*sizeof(char),93);
				strcpy(system->pdb_input, token2);
			}
		}

		if(!strcasecmp(token1, "pdb_output")) {
			if(!system->pdb_output) {
				system->pdb_output = calloc(MAXLINE, sizeof(char));
				memnullcheck(system->pdb_output,MAXLINE*sizeof(char),94);
				strcpy(system->pdb_output, token2);
			}
		}

		if(!strcasecmp(token1, "pdb_restart")) {
			if(!system->pdb_restart) {
				system->pdb_restart = calloc(MAXLINE, sizeof(char));
				memnullcheck(system->pdb_restart,MAXLINE*sizeof(char),95);
				strcpy(system->pdb_restart, token2);
			}
		}

		if(!strcasecmp(token1, "traj_output")) {
			if(!system->traj_output) {
				system->traj_output = calloc(MAXLINE, sizeof(char));
				memnullcheck(system->traj_output,MAXLINE*sizeof(char),96);
				strcpy(system->traj_output, token2);
			}
		}

		if(!strcasecmp(token1, "energy_output")) {
			if(!system->energy_output) {
				system->energy_output = calloc(MAXLINE, sizeof(char));
				memnullcheck(system->energy_output,MAXLINE*sizeof(char),97);
				strcpy(system->energy_output, token2);
			}
		}

		if(!strcasecmp(token1, "pop_histogram_output")) {
			if(!system->histogram_output) {
				system->histogram_output = calloc(MAXLINE, sizeof(char));
				memnullcheck(system->histogram_output,MAXLINE*sizeof(char),98);
				strcpy(system->histogram_output, token2);
			}
		}

		if(system->polarization && !strcasecmp(token1, "dipole_output")) {
			if(!system->dipole_output) {
				system->dipole_output = calloc(MAXLINE, sizeof(char));
				memnullcheck(system->dipole_output,MAXLINE*sizeof(char),99);
				strcpy(system->dipole_output, token2);
			}
		}

		if(system->polarization && !strcasecmp(token1, "field_output")) {
			if(!system->field_output) {
				system->field_output = calloc(MAXLINE, sizeof(char));
				memnullcheck(system->field_output,MAXLINE*sizeof(char),100);
				strcpy(system->field_output, token2);
			}
		}

		if(!strcasecmp(token1, "frozen_output")) {
			if(!system->frozen_output) {
				system->frozen_output = calloc(MAXLINE, sizeof(char));
				memnullcheck(system->frozen_output,MAXLINE*sizeof(char),101);
				strcpy(system->frozen_output, token2);
			}
		}

		if(!strcasecmp(token1, "insert_input")) {
			if(!system->insert_input) {
				system->insert_input = calloc(MAXLINE, sizeof(char));
				memnullcheck(system->insert_input,MAXLINE*sizeof(char),102);
				strcpy(system->insert_input, token2);
			}
		}






                // Surface fit input parameters
                //////////////////////////////////////////////////

                if( !strcasecmp( token1, "fit_schedule")) {
                    register double temp = atof(token2);
                    if( temp<= 0.0 || temp>=1.0 )
                        error( "INPUT: Invalid schedule. Using default.\n");
                    else
                        system->fit_schedule = temp;
                }
                else if( !strcasecmp( token1, "fit_max_energy")) {
                    register double temp = atof(token2);
                    if( temp<= 0.0 )
                        error( "INPUT: fit_max_energy parameter must be greater than zero. Using default.\n");
                    else
                        system->fit_max_energy = temp;
                }
                else if( !strcasecmp( token1, "fit_start_temp" )) {
                    register double temp = atof(token2);
                    if( temp<= 0.0 )
                        error( "INPUT: fit_start_temp parameter must be greater than zero. Using default.\n");
                    else
                        system->fit_start_temp = atof( token2 );
                }
                
                else if(!strcasecmp(token1, "fit_input")) {

                    // navigate to the end of the input file linked list
                    fileNode_t *node = &(system->fit_input_list);
                    while( node->next )
                        node = node->next;

                    // allocate a new node
                    if( !(node->next = malloc( sizeof(fileNode_t)   ))   )
                    {
                        error( "INPUT: Exhausted memory during input file node allocation.\n");
                        return (NULL);
                    }

                    // advance to new node and initialize
                    node = node->next;
                    node->next = 0; // terminate list
                    if(   !(node->data.filename = calloc(MAXLINE, sizeof(char)))   )
                    {
                        error( "INPUT: Exhausted memory during string allocation for fit input filename.\n");
                        return (NULL);
                    }

                    // copy filename to node and increment list count
                    strcpy(node->data.filename, token2);
                    system->fit_input_list.data.count++;
                }




		if(!strcasecmp(token1, "basis1")) {
			system->pbc->basis[0][0] = atof(token2);
			system->pbc->basis[0][1] = atof(token3);
			system->pbc->basis[0][2] = atof(token4);
		}

		if(!strcasecmp(token1, "basis2")) {
			system->pbc->basis[1][0] = atof(token2);
			system->pbc->basis[1][1] = atof(token3);
			system->pbc->basis[1][2] = atof(token4);
		}

		if(!strcasecmp(token1, "basis3")) {
			system->pbc->basis[2][0] = atof(token2);
			system->pbc->basis[2][1] = atof(token3);
			system->pbc->basis[2][2] = atof(token4);
		}

		if(!strcasecmp(token1, "max_bondlength")) {
			system->max_bondlength = atof(token2);
		}

		if(!strcasecmp(token1, "pop_histogram")) {
			if(!strcasecmp(token2, "on")) 
				system->calc_hist=1;
			else 
				system->calc_hist=0;
		}

		if(!strcasecmp(token1, "pop_hist_resolution"))
			system->hist_resolution=atof(token2);

		

		memset(linebuffer, 0, MAXLINE);
		n = fgets(linebuffer, MAXLINE, fp);
	}

	/* close the config file */
	fclose(fp);

	/* calculate things related to the periodic boundary conditions */
	pbc(system->pbc);

	return(system);

}

double h2_fugacity(double temperature, double pressure) {

	if((temperature == 77.0) && (pressure <= 200.0)) {

		output("INPUT: fugacity calculation using Zhou function\n");
		return(h2_fugacity_zhou(temperature, pressure));

	}  else if(temperature >= 273.15) {

		output("INPUT: fugacity calculation using Shaw function\n");
		return(h2_fugacity_shaw(temperature, pressure));

	} else {

		output("INPUT: fugacity calculation using BACK EoS\n");
		return(h2_fugacity_back(temperature, pressure));

	}

	return(0); /* NOT REACHED */

}

/* use the semi-empirical BACK equation of state */
/* Tomas Boublik, "The BACK equation of state for hydrogen and related compounds", Fluid Phase Equilibria, 240, 96-100 (2005) */

double h2_fugacity_back(double temperature, double pressure) {

	double fugacity_coefficient, fugacity;
	double comp_factor;
	double P, dP;
	char linebuf[MAXLINE];

	/* integrate (z-1)/P from 0 to P */
	fugacity_coefficient = 0;
	for(P = 0.001, dP = 0.001; P <= pressure; P += dP) {

		comp_factor = h2_comp_back(temperature, P);
		fugacity_coefficient += dP*(comp_factor - 1.0)/P;

	}
	fugacity_coefficient = exp(fugacity_coefficient);

	comp_factor = h2_comp_back(temperature, pressure);
	sprintf(linebuf, "INPUT: BACK compressibility factor at %.3f atm is %.3f\n", pressure, comp_factor);
	output(linebuf);

	fugacity = pressure*fugacity_coefficient;
	return(fugacity);

}


#define BACK_H2_ALPHA	1.033
#define BACK_H2_U0	38.488
#define BACK_H2_V00	9.746
#define BACK_H2_N	0.00
#define BACK_C		0.12

#define BACK_MAX_M	9
#define BACK_MAX_N	4

double h2_comp_back(double temperature, double pressure) {

	double alpha, y;				/* repulsive part of the compressibility factor */
	double V, V0, u, D[BACK_MAX_M][BACK_MAX_N];	/* attractive part */
	int n, m;					/* indices for double sum of attractive part */
	double comp_factor;
	double comp_factor_repulsive;
	double comp_factor_attractive;
	double fugacity;


	/* setup the BACK universal D constants */
	D[0][0] = -8.8043;	D[0][1] = 2.9396;	D[0][2] = -2.8225;	D[0][3] = 0.34;
	D[1][0] = 4.164627;	D[1][1] = -6.0865383;	D[1][2] = 4.7600148;	D[1][3] = -3.1875014;
	D[2][0] = -48.203555;	D[2][1] = 40.137956;	D[2][2] = 11.257177;	D[2][3] = 12.231796;
	D[3][0] = 140.4362;	D[3][1] = -76.230797;	D[3][2] = -66.382743;	D[3][3] = -12.110681;
	D[4][0] = -195.23339;	D[4][1] = -133.70055;	D[4][2] = 69.248785;	D[4][3] = 0.0;
	D[5][0] = 113.515;	D[5][1] = 860.25349;	D[5][2] = 0.0;		D[5][3] = 0.0;
	D[6][0] = 0.0;		D[6][1] = -1535.3224;	D[6][2] = 0.0;		D[6][3] = 0.0;
	D[7][0] = 0.0;		D[7][1] = 1221.4261;	D[7][2] = 0.0;		D[7][3] = 0.0;
	D[8][0] = 0.0;		D[8][1] = -409.10539;	D[8][2] = 0.0;		D[8][3] = 0.0;

	/* calculate attractive part */
	V0 = BACK_H2_V00*(1.0 - BACK_C*exp(-3.0*BACK_H2_U0/temperature));
	V = NA*KB*temperature/(pressure*ATM2PASCALS*1.0e-6);
	u = BACK_H2_U0*(1.0 + BACK_H2_N/temperature);

	comp_factor_attractive = 0;
	for(n = 0; n < BACK_MAX_N; n++)
		for(m = 0; m < BACK_MAX_M; m++)
			comp_factor_attractive += ((double)(m+1))*D[m][n]*pow(u/temperature, ((double)(n+1)))*pow(V0/V, ((double)(m+1)));

	/* calculate repulsive part */
	alpha = BACK_H2_ALPHA;
	y = (M_PI*sqrt(2.0)/6.0)*(pressure*ATM2PASCALS*1.0e-6)/(NA*KB*temperature)*V0;
	comp_factor_repulsive = 1.0 + (3.0*alpha - 2.0)*y;
	comp_factor_repulsive += (3.0*pow(alpha, 2.0) - 3.0*alpha + 1.0)*pow(y, 2.0);
	comp_factor_repulsive -= pow(alpha, 2.0)*pow(y, 3.0);
	comp_factor_repulsive /= pow((1.0 - y), 3.0);

	comp_factor = comp_factor_repulsive + comp_factor_attractive;
	return(comp_factor);

}


/* calculate the fugacity correction for H2 for 0 C and higher */
/* this empirical relation follows from: */
/* H.R. Shaw, D.F. Wones, American Journal of Science, 262, 918-929 (1964) */
double h2_fugacity_shaw(double temperature, double pressure) {

	double C1, C2, C3;
	double fugacity, fugacity_coefficient;

	C1 = -3.8402*pow(temperature, 1.0/8.0) + 0.5410;
	C1 = exp(C1);

	C2 = -0.1263*pow(temperature, 1.0/2.0) - 15.980;
	C2 = exp(C2);

	C3 = -0.11901*temperature - 5.941;
	C3 = exp(C3);
	C3 *= 300.0;

	fugacity_coefficient = C1*pressure - C2*pow(pressure, 2.0) + C3*exp(-pressure/300.0 - 1.0);
	fugacity_coefficient = exp(fugacity_coefficient);
	fugacity = fugacity_coefficient*pressure;

	return(fugacity);

}

/* fugacity for low temperature and up to 200 atm */
/* Zhou, Zhou, Int. J. Hydrogen Energy, 26, 597-601 (2001) */
double h2_fugacity_zhou(double temperature, double pressure) {

	double fugacity, fugacity_coefficient;

	pressure *= ATM2PSI;

	fugacity_coefficient = -1.38130e-4*pressure;
	fugacity_coefficient += 4.67096e-8*pow(pressure, 2.0)/2;
	fugacity_coefficient += 5.93690e-12*pow(pressure, 3.0)/3;
	fugacity_coefficient += -3.24527e-15*pow(pressure, 4.0)/4;
	fugacity_coefficient += 3.54211e-19*pow(pressure, 5.0)/5;

	pressure /= ATM2PSI;

	fugacity_coefficient = exp(fugacity_coefficient);
	fugacity = pressure*fugacity_coefficient;

	return(fugacity);

}

/* ***************************** CH4 EQUATION OF STATE *************************************** */
double ch4_fugacity(double temperature, double pressure) {

        if((temperature >= 298.0) && (temperature <= 300.0) && (pressure <= 500.0)) {

                output("INPUT: CH4 fugacity calculation using BACK EoS\n");
                return(ch4_fugacity_back(temperature, pressure));

        } else if((temperature == 150.0) && (pressure <= 200.0)) {

                output("INPUT: CH4 fugacity calculation using Peng-Robinson EoS\n");
                return(ch4_fugacity_PR(temperature, pressure));

        } else {

                output("INPUT: Unknown if CH4 fugacity will be correct at the requested temperature & pressure...defaulting to use the BACK EoS.\n");
                return(ch4_fugacity_back(temperature, pressure));

        }

        return(0); /* NOT REACHED */

}

/* Incorporate BACK EOS */
double ch4_fugacity_back(double temperature, double pressure) {

        double fugacity_coefficient, fugacity;
        double comp_factor;
        double P, dP;
        char linebuf[MAXLINE];

        /* integrate (z-1)/P from 0 to P */
        fugacity_coefficient = 0;
        for(P = 0.001, dP = 0.001; P <= pressure; P += dP) {

                comp_factor = ch4_comp_back(temperature, P);
                fugacity_coefficient += dP*(comp_factor - 1.0)/P;

        }
        fugacity_coefficient = exp(fugacity_coefficient);

        comp_factor = ch4_comp_back(temperature, pressure);
        sprintf(linebuf, "INPUT: CH4 BACK compressibility factor at %.3f atm is %.3f\n", pressure, comp_factor);
        output(linebuf);

        fugacity = pressure*fugacity_coefficient;
        return(fugacity);
}

#define MWCH4 16.043
#define BACK_CH4_ALPHA  1.000
#define BACK_CH4_U0     188.047
#define BACK_CH4_V00    21.532
#define BACK_CH4_N      2.40
#define BACK_C          0.12

#define BACK_MAX_M      9
#define BACK_MAX_N      4

double ch4_comp_back(double temperature, double pressure) {

        double alpha, y;                                /* repulsive part of the compressibility factor */
        double V, V0, u, D[BACK_MAX_M][BACK_MAX_N];     /* attractive part */
        int n, m;                                       /* indices for double sum of attractive part */
        double comp_factor;
        double comp_factor_repulsive;
        double comp_factor_attractive;
        double fugacity;

        /* setup the BACK universal D constants */
        D[0][0] = -8.8043;      D[0][1] = 2.9396;       D[0][2] = -2.8225;      D[0][3] = 0.34;
        D[1][0] = 4.164627;     D[1][1] = -6.0865383;   D[1][2] = 4.7600148;    D[1][3] = -3.1875014;
        D[2][0] = -48.203555;   D[2][1] = 40.137956;    D[2][2] = 11.257177;    D[2][3] = 12.231796;
        D[3][0] = 140.4362;     D[3][1] = -76.230797;   D[3][2] = -66.382743;   D[3][3] = -12.110681;
        D[4][0] = -195.23339;   D[4][1] = -133.70055;   D[4][2] = 69.248785;    D[4][3] = 0.0;
        D[5][0] = 113.515;      D[5][1] = 860.25349;    D[5][2] = 0.0;          D[5][3] = 0.0;
        D[6][0] = 0.0;          D[6][1] = -1535.3224;   D[6][2] = 0.0;          D[6][3] = 0.0;
        D[7][0] = 0.0;          D[7][1] = 1221.4261;    D[7][2] = 0.0;          D[7][3] = 0.0;
        D[8][0] = 0.0;          D[8][1] = -409.10539;   D[8][2] = 0.0;          D[8][3] = 0.0;

        /* calculate attractive part */
        V0 = BACK_CH4_V00*(1.0 - BACK_C*exp(-3.0*BACK_CH4_U0/temperature));
        V = NA*KB*temperature/(pressure*ATM2PASCALS*1.0e-6);
        u = BACK_CH4_U0*(1.0 + BACK_CH4_N/temperature);

        comp_factor_attractive = 0;
        for(n = 0; n < BACK_MAX_N; n++)
                for(m = 0; m < BACK_MAX_M; m++)
                        comp_factor_attractive += ((double)(m+1))*D[m][n]*pow(u/temperature, ((double)(n+1)))*pow(V0/V, ((double)(m+1)));

        /* calculate repulsive part */
        alpha = BACK_CH4_ALPHA;
        y = (M_PI*sqrt(2.0)/6.0)*(pressure*ATM2PASCALS*1.0e-6)/(NA*KB*temperature)*V0;
        comp_factor_repulsive = 1.0 + (3.0*alpha - 2.0)*y;
        comp_factor_repulsive += (3.0*pow(alpha, 2.0) - 3.0*alpha + 1.0)*pow(y, 2.0);
        comp_factor_repulsive -= pow(alpha, 2.0)*pow(y, 3.0);
        comp_factor_repulsive /= pow((1.0 - y), 3.0);

        comp_factor = comp_factor_repulsive + comp_factor_attractive;
        return(comp_factor);

}

/* Apply the Peng-Robinson EoS to methane */
double ch4_fugacity_PR(double temperature, double pressure) {

  double Z, A, B, aa, bb, TcCH4, PcCH4, Tr;
  double alpha, alpha2, wCH4, R, Q, X, j, k, l;
  double theta, X2, Q3;
  double uu, U, V, root1, root2, root3, answer, stuff1, stuff2, stuff3;
  double f1, f2, f3, f4, fugacity, lnfoverp;
  double pi=acos(-1.0);

  /*Peng Robinson variables and equations for CH4 units K,atm, L, mole*/
     TcCH4 = 190.564;    /* K */
     PcCH4 = 45.391;    /* atm */
     wCH4  = 0.01142;
     R    = 0.08206;   /* gas constant atmL/moleK */

     aa = (0.45724*R*R*TcCH4*TcCH4) / PcCH4;
     bb = (0.07780*R*TcCH4) / PcCH4;
     Tr = temperature / TcCH4;
     stuff1 = 0.37464 + 1.54226*wCH4 - 0.26992*wCH4*wCH4;
     stuff2=1.0-sqrt(Tr);
     alpha=1.0+stuff1*stuff2;
     alpha2=alpha*alpha;
     A=alpha2*aa*pressure/(R*R*temperature*temperature);
     B=bb*pressure/(R*temperature);

     /* solving a cubic equation part */
     j=-1.0*(1-B);
     k=A-3.0*B*B-2.0*B;
     l= -1*(A*B- B*B -B*B*B);
     Q=(j*j-3.0*k)/9.0;
     X=(2.0*j*j*j -9.0*j*k+27.0*l)/54.0;
     Q3=Q*Q*Q;
     X2=X*X;

     /* Need to check X^2 < Q^3 */
     if((X*X)<(Q*Q*Q)){    /* THREE REAL ROOTS  */
       theta=acos((X/sqrt(Q3)));
       root1=-2.0*sqrt(Q)*cos(theta/3.0)-j/3.0;
       root2=-2.0*sqrt(Q)*cos((theta+2.0*pi)/3.0)-j/3.0;
       root3=-2.0*sqrt(Q)*cos((theta-2.0*pi)/3.0)-j/3.0;

       /*Choose the root closest to 1, which is "ideal gas law" */
       if((1.0-root1)<(1.0-root2) && (1.0-root1)<(1.0-root3))
         Z=root1;
       else if((1.0-root2)<(1.0-root3) && (1.0-root2)<(1.0-root1))
         Z=root2;
       else
         Z=root3;
     }
     else{   /* ONLY ONE real root */
       stuff3= X*X-Q*Q*Q;
       uu=X-sqrt(stuff3);
       /*Power function must have uu a positive number*/
       if(uu<0.0)
         uu=-1.0*uu;
       U=pow(uu,(1.0/3.0));
       V=Q/U;
       root1=U+V-j/3.0;
       Z=root1;
     }

     /* using Z calculate the fugacity */
     f1=(Z-1.0)-log(Z-B);
     f2=A/(2.0*sqrt(2.0)*B);
     f3=Z+(1.0+sqrt(2.0))*B;
     f4=Z+(1.0-sqrt(2.0))*B;
     lnfoverp=f1-f2*log(f3/f4);
     fugacity=exp(lnfoverp)*pressure;

  return(fugacity);
}
/* ******************************* END CH4 FUGACITY ****************************************** */

/* *************************** N2 BACK EQUATION OF STATE ************************************* */
double n2_fugacity(double temperature, double pressure) {

        if((temperature == 78.0) && (pressure <= 1.0)) {

                output("INPUT: N2 fugacity calculation using Zhou\n");
                return(n2_fugacity_zhou(temperature, pressure));

        } else if((temperature == 78.0) && (pressure >= 10.0) && (pressure <= 300.0)) {

                output("INPUT: N2 fugacity calculation using Peng-Robinson EoS\n");
                return(n2_fugacity_PR(temperature, pressure));

        } else if((temperature == 150.0) && (pressure < 175.0)) {

                output("INPUT: N2 fugacity calculation using Peng-Robinson EoS\n");
                return(n2_fugacity_PR(temperature, pressure));

        } else if((temperature == 150.0) && (pressure >= 175.0) && (pressure <= 325.0)) {

                output("INPUT: N2 fugacity calculation using BACK EoS\n");
                return(n2_fugacity_back(temperature, pressure));

        } else if((temperature >= 298.0) && (temperature <= 300.0) && (pressure <= 350.0)) {

                output("INPUT: N2 fugacity calculation using Peng-Robinson EoS\n");
                return(n2_fugacity_PR(temperature, pressure));

        } else {
                output("INPUT: Unknown if N2 fugacity will be correct at the requested temperature & pressure...defaulting to use the PR EoS.\n");
                return(n2_fugacity_PR(temperature, pressure));

        }

        return(0); /* NOT REACHED */

}

/* Incorporate BACK EOS */
double n2_fugacity_back(double temperature, double pressure) {

        double fugacity_coefficient, fugacity;
        double comp_factor;
        double P, dP;
        char linebuf[MAXLINE];

        /* integrate (z-1)/P from 0 to P */
        fugacity_coefficient = 0;
        for(P = 0.001, dP = 0.001; P <= pressure; P += dP) {

                comp_factor = n2_comp_back(temperature, P);
                fugacity_coefficient += dP*(comp_factor - 1.0)/P;

        }
        fugacity_coefficient = exp(fugacity_coefficient);

        comp_factor = n2_comp_back(temperature, pressure);
        sprintf(linebuf, "INPUT: BACK compressibility factor at %.3f atm is %.3f\n", pressure, comp_factor);
        output(linebuf);

        fugacity = pressure*fugacity_coefficient;
        return(fugacity);

}

#define BACK_N2_ALPHA   1.048
#define BACK_N2_U0      120.489
#define BACK_N2_V00     18.955
#define BACK_N2_N       10.81
#define BACK_C          0.12

#define BACK_MAX_M      9
#define BACK_MAX_N      4

double n2_comp_back(double temperature, double pressure) {

        double alpha, y;                                /* repulsive part of the compressibility factor */
        double V, V0, u, D[BACK_MAX_M][BACK_MAX_N];     /* attractive part */
        int n, m;                                       /* indices for double sum of attractive part */
        double comp_factor;
        double comp_factor_repulsive;
        double comp_factor_attractive;
        double fugacity;


        /* setup the BACK universal D constants */
        D[0][0] = -8.8043;      D[0][1] = 2.9396;       D[0][2] = -2.8225;      D[0][3] = 0.34;
        D[1][0] = 4.164627;     D[1][1] = -6.0865383;   D[1][2] = 4.7600148;    D[1][3] = -3.1875014;
        D[2][0] = -48.203555;   D[2][1] = 40.137956;    D[2][2] = 11.257177;    D[2][3] = 12.231796;
        D[3][0] = 140.4362;     D[3][1] = -76.230797;   D[3][2] = -66.382743;   D[3][3] = -12.110681;
        D[4][0] = -195.23339;   D[4][1] = -133.70055;   D[4][2] = 69.248785;    D[4][3] = 0.0;
        D[5][0] = 113.515;      D[5][1] = 860.25349;    D[5][2] = 0.0;          D[5][3] = 0.0;
        D[6][0] = 0.0;          D[6][1] = -1535.3224;   D[6][2] = 0.0;          D[6][3] = 0.0;
        D[7][0] = 0.0;          D[7][1] = 1221.4261;    D[7][2] = 0.0;          D[7][3] = 0.0;
        D[8][0] = 0.0;          D[8][1] = -409.10539;   D[8][2] = 0.0;          D[8][3] = 0.0;

        /* calculate attractive part */
        V0 = BACK_N2_V00*(1.0 - BACK_C*exp(-3.0*BACK_N2_U0/temperature));
        V = NA*KB*temperature/(pressure*ATM2PASCALS*1.0e-6);
        u = BACK_N2_U0*(1.0 + BACK_N2_N/temperature);

        comp_factor_attractive = 0;
        for(n = 0; n < BACK_MAX_N; n++)
                for(m = 0; m < BACK_MAX_M; m++)
                        comp_factor_attractive += ((double)(m+1))*D[m][n]*pow(u/temperature, ((double)(n+1)))*pow(V0/V, ((double)(m+1)));

        /* calculate repulsive part */
        alpha = BACK_N2_ALPHA;
        y = (M_PI*sqrt(2.0)/6.0)*(pressure*ATM2PASCALS*1.0e-6)/(NA*KB*temperature)*V0;
        comp_factor_repulsive = 1.0 + (3.0*alpha - 2.0)*y;
        comp_factor_repulsive += (3.0*pow(alpha, 2.0) - 3.0*alpha + 1.0)*pow(y, 2.0);
        comp_factor_repulsive -= pow(alpha, 2.0)*pow(y, 3.0);
        comp_factor_repulsive /= pow((1.0 - y), 3.0);

        comp_factor = comp_factor_repulsive + comp_factor_attractive;
        return(comp_factor);
}

/* Apply the Peng-Robinson EoS to N2 */
double n2_fugacity_PR(double temperature, double pressure) {

  double Z, A, B, aa, bb, TcN2, PcN2, Tr;
  double alpha, alpha2, wN2, R, Q, X, j, k, l;
  double theta, X2, Q3;
  double uu, U, V, root1, root2, root3, answer, stuff1, stuff2, stuff3;
  double f1, f2, f3, f4, fugacity, lnfoverp;
  double pi=acos(-1.0);

  /*Peng Robinson variables and equations for N2 units K,atm, L, mole*/
     TcN2 = 126.192;   /* K */
     PcN2 = 33.514;    /* atm */
     wN2  = 0.037;
     R    = 0.08206;   /* gas constant atmL/moleK */

     aa = (0.45724*R*R*TcN2*TcN2) / PcN2;
     bb = (0.07780*R*TcN2) / PcN2;
     Tr = temperature / TcN2;
     stuff1 = 0.37464 + 1.54226*wN2 - 0.26992*wN2*wN2;
     stuff2=1.0-sqrt(Tr);
     alpha=1.0+stuff1*stuff2;
     alpha2=alpha*alpha;
     A=alpha2*aa*pressure/(R*R*temperature*temperature);
     B=bb*pressure/(R*temperature);

     /* solving a cubic equation part */
     j=-1.0*(1-B);
     k=A-3.0*B*B-2.0*B;
     l= -1*(A*B- B*B -B*B*B);
     Q=(j*j-3.0*k)/9.0;
     X=(2.0*j*j*j -9.0*j*k+27.0*l)/54.0;
     Q3=Q*Q*Q;
     X2=X*X;

     /* Need to check X^2 < Q^3 */
     if((X*X)<(Q*Q*Q)){    /* THREE REAL ROOTS  */
       theta=acos((X/sqrt(Q3)));
       root1=-2.0*sqrt(Q)*cos(theta/3.0)-j/3.0;
       root2=-2.0*sqrt(Q)*cos((theta+2.0*pi)/3.0)-j/3.0;
       root3=-2.0*sqrt(Q)*cos((theta-2.0*pi)/3.0)-j/3.0;

       /*Choose the root closest to 1, which is "ideal gas law" */
       if((1.0-root1)<(1.0-root2) && (1.0-root1)<(1.0-root3))
         Z=root1;
       else if((1.0-root2)<(1.0-root3) && (1.0-root2)<(1.0-root1))
         Z=root2;
       else
         Z=root3;
     }
     else{   /* ONLY ONE real root */
       stuff3= X*X-Q*Q*Q;
       uu=X-sqrt(stuff3);
       /*Power function must have uu a positive number*/
       if(uu<0.0)
         uu=-1.0*uu;
       U=pow(uu,(1.0/3.0));
       V=Q/U;
       root1=U+V-j/3.0;
       Z=root1;
     }

     /* using Z calculate the fugacity */
     f1=(Z-1.0)-log(Z-B);
     f2=A/(2.0*sqrt(2.0)*B);
     f3=Z+(1.0+sqrt(2.0))*B;
     f4=Z+(1.0-sqrt(2.0))*B;
     lnfoverp=f1-f2*log(f3/f4);
     fugacity=exp(lnfoverp)*pressure;

  return(fugacity);
}

/* Apply the Zhou function to N2 */
double n2_fugacity_zhou(double temperature, double pressure) {

        double fugacity_coefficient, fugacity;

        output("INPUT: N2 fugacity calculation using Zhou function\n");

        pressure *= ATM2PSI;

        fugacity_coefficient = -1.38130e-4*pressure;
        fugacity_coefficient += 4.67096e-8*pow(pressure, 2.0)/2;
        fugacity_coefficient += 5.93690e-12*pow(pressure, 3.0)/3;
        fugacity_coefficient += -3.24527e-15*pow(pressure, 4.0)/4;
        fugacity_coefficient += 3.54211e-19*pow(pressure, 5.0)/5;

        pressure /= ATM2PSI;

        fugacity_coefficient = exp(fugacity_coefficient);
        fugacity = pressure*fugacity_coefficient;

        return(fugacity);
}
/* ********************************** END N2 ****************************************************** */


/* reads in temperature in K, and pressure (of the ideal gas in the resevoir) in atm */
/* return the CO2 fugacity via the Peng-Robinson equation of state */
/* else return 0.0 on error - I don't have an error statement*/
/* units are  atm, K,  */

double co2_fugacity(double temperature, double pressure) {

  double Z, A, B, aa, bb, Tc, Pc, Tr;
  double alpha,alpha2, w, R, Q, X, j,k,l;
  double theta, X2, Q3;
  double uu,U,V,root1,root2,root3, answer, stuff1, stuff2, stuff3;
  double f1, f2, f3, f4, fugacity, lnfoverp;
  double pi=acos(-1.0);

  /*Peng Robinson variables and equations for CO2 units K,atm, L, mole*/
     Tc=304.12;   /*K  */
     Pc=73.74/1.01325;    /*bar to atm*/
     w=0.225;
     R=0.08206;   /* gas constant atmL/moleK */
  
     aa=0.45724*R*R*Tc*Tc/Pc;
     bb=0.07780*R*Tc/Pc;
     Tr=temperature/Tc;
     stuff1=0.37464+1.54226*w -0.26992*w*w;
     stuff2=1.0-sqrt(Tr);
     alpha=1.0+stuff1*stuff2;
     alpha2=alpha*alpha;
     A=alpha2*aa*pressure/(R*R*temperature*temperature);
     B=bb*pressure/(R*temperature);
        
     /* solving a cubic equation part */
     j=-1.0*(1-B);
     k=A-3.0*B*B-2.0*B;
     l= -1*(A*B- B*B -B*B*B);
     Q=(j*j-3.0*k)/9.0;
     X=(2.0*j*j*j -9.0*j*k+27.0*l)/54.0;
     Q3=Q*Q*Q;
     X2=X*X;
     
     /* Need to check X^2 < Q^3 */
     if((X*X)<(Q*Q*Q)){    /* THREE REAL ROOTS  */
       theta=acos((X/sqrt(Q3)));
       root1=-2.0*sqrt(Q)*cos(theta/3.0)-j/3.0;
       root2=-2.0*sqrt(Q)*cos((theta+2.0*pi)/3.0)-j/3.0;
       root3=-2.0*sqrt(Q)*cos((theta-2.0*pi)/3.0)-j/3.0;
       
       /*Choose the root closest to 1, which is "ideal gas law" */
       if((1.0-root1)<(1.0-root2) && (1.0-root1)<(1.0-root3))
	 Z=root1;
       else if((1.0-root2)<(1.0-root3) && (1.0-root2)<(1.0-root1))
	 Z=root2;
       else
	 Z=root3;
     }
     else{   /* ONLY ONE real root */
       stuff3= X*X-Q*Q*Q;
       uu=X-sqrt(stuff3);
       /*Power function must have uu a positive number*/
       if(uu<0.0)
	 uu=-1.0*uu;
       U=pow(uu,(1.0/3.0));
       V=Q/U;
       root1=U+V-j/3.0;
       Z=root1;
     }

     /* using Z calculate the fugacity */
     f1=(Z-1.0)-log(Z-B);
     f2=A/(2.0*sqrt(2.0)*B);
     f3=Z+(1.0+sqrt(2.0))*B;
     f4=Z+(1.0-sqrt(2.0))*B;
     lnfoverp=f1-f2*log(f3/f4);
     fugacity=exp(lnfoverp)*pressure;

  return(fugacity);

}


int check_system(system_t *system) {

	char linebuf[MAXLINE];

	switch(system->ensemble) {

		case ENSEMBLE_UVT:
			output("INPUT: Grand canonical ensemble\n");
			break;
		case ENSEMBLE_NVT:
			output("INPUT: Canonical ensemble\n");
			break;
		case ENSEMBLE_SURF:
			output("INPUT: Potential energy surface\n");
			break;
		case ENSEMBLE_SURF_FIT:
			output("INPUT: Potential energy surface fitting\n");
			break;
		case ENSEMBLE_NVE:
			output("INPUT: Microcanonical ensemble\n");
			break;
		case ENSEMBLE_TE:
			output("INPUT: Single-point energy calculation\n");
			break;
		default:
			error("INPUT: improper ensemble specified\n");
			return(-1);
	}

	if(system->ensemble == ENSEMBLE_SURF_FIT) {
		// Record number of curves for convenient reference
		int nCurves = system->fit_input_list.data.count;

		if( !nCurves ) {
			error( "INPUT: There were no fit_input files specified in the main input file.\n" );
			return (-1);
		}

		if( nCurves < 2 ) {
			error( "INPUT: There were less than two fit_input files specified in the main input file.\n" );
			error( "       A minmum of two fit_input files are required for surface-fit calculations.\n" );
			return (-1);
		}
	}

	if(system->ensemble == ENSEMBLE_SURF) {

		output("INPUT: surface module activated\n");
		if(system->surf_max < system->surf_min) {
			error("INPUT: surf_max is greater than surf_min\n");
			return(-1);
		} else {
			sprintf(linebuf, "INPUT: minimum surface coordinate is %.3f\n", system->surf_min);
			output(linebuf);
			sprintf(linebuf, "INPUT: maximum surface coordinate is %.3f\n", system->surf_max);
			output(linebuf);
		}

		if(system->surf_inc <= 0.0) {
			error("INPUT: surf_inc is less than or equal to 0\n");
			return(-1);
		} else {
			sprintf(linebuf, "INPUT: incremental surface displacement coordinate is %.3f\n", system->surf_inc);
			output(linebuf);
		}

		if(!system->surf_preserve && (system->surf_ang <= 0.0)) {
			error("INPUT: surf_ang is less than or equal to 0\n");
			return(-1);
		} else {
			sprintf(linebuf, "INPUT: incremental surface angle coordinate is %.3f\n", system->surf_ang);
			output(linebuf);
		}


	}

	if(system->rd_only) output("INPUT: calculating repulsion/dispersion only\n");

	if(system->gwp) {
		output("INPUT: Gaussian wavepacket code active\n");

		if(system->gwp_probability == 0.) {

			output("INPUT: GWP move scaling not input - setting equal to move_probability\n");
			system->gwp_probability = system->move_probability;
		}
	}


	if(system->wolf) output("INPUT: Wolf summation active\n");

	if(system->rd_lrc)
		output("INPUT: rd long-range corrections are ON\n");
	else
		output("INPUT: rd long-range corrections are OFF\n");

	if(system->sg) output("INPUT: Molecular potential is Silvera-Goldman\n");

	if(system->dreiding) output("INPUT: Molecular potential is DREIDING\n");

	sprintf(linebuf, "INPUT: frozen atom charges scaled by %.2f\n", system->scale_charge);
	output(linebuf);

	sprintf(linebuf, "INPUT: frozen atom rd scaled by %.2f\n", system->scale_rd);
	output(linebuf);

	if(system->spectre) {

		if(system->ensemble != ENSEMBLE_NVT) {
			error("INPUT: SPECTRE algorithm requires canonical ensemble\n");
			return(-1);
		} else {

			output("INPUT: SPECTRE algorithm activated\n");

			sprintf(linebuf, "INPUT: SPECTRE max charge = %.3f\n", system->spectre_max_charge);
			output(linebuf);

			sprintf(linebuf, "INPUT: SPECTRE max target = %.3f\n", system->spectre_max_target);
			output(linebuf);

		}

	}

	if(system->feynman_hibbs) {

		output("INPUT: Feynman-Hibbs effective potential activated\n");

		if(system->feynman_kleinert) {
			output("INPUT: Feynman-Kleinert iteration method activated\n");

			if(!system->rd_anharmonic) {
				error("INPUT: Feynman-Kleinert iteration only implemented for anharmonic oscillator\n");
				return(-1);
			}

		} else {

			switch(system->feynman_hibbs_order) {

				case 2:
					sprintf(linebuf, "INPUT: Feynman-Hibbs second-order quantum correction activated\n");
					output(linebuf);
					break;
				case 4:
					sprintf(linebuf, "INPUT: Feynman-Hibbs fourth-order quantum correction activated\n");
					output(linebuf);
					break;
				default:
					output("INPUT: Feynman-Hibbs order unspecified - defaulting to h^2\n");
					system->feynman_hibbs_order = 2;
					break;

			}
		}

	}

#ifdef QM_ROTATION
	if(system->quantum_rotation){
		output("INPUT: Quantum rotational eigenspectrum calculation enabled\n");
		if(system->quantum_rotation_B <= 0.0) {
			error("INPUT: invalid quantum rotational constant B specified\n");
			return(-1);
		} else {
			sprintf(linebuf, "INPUT: Quantum rotational constant B = %.3f K (%.3f cm^-1)\n", system->quantum_rotation_B, system->quantum_rotation_B*KB/(100.0*H*C));
			output(linebuf);
		}

		if(system->quantum_rotation_level_max <= 0) {
			error("INPUT: invalid quantum rotation level max\n");
			return(-1);
		} else {
			sprintf(linebuf, "INPUT: Quantum rotation level max = %d\n", system->quantum_rotation_level_max);
			output(linebuf);
		}

		if(system->quantum_rotation_l_max <= 0) {
			error("INPUT: invalid quantum rotation l_max\n");
			return(-1);
		} else {
			sprintf(linebuf, "INPUT: Quantum rotation l_max = %d\n", system->quantum_rotation_l_max);
			output(linebuf);
		}

		if(system->quantum_rotation_level_max > (system->quantum_rotation_l_max+1)*(system->quantum_rotation_l_max+1)) {
			error("INPUT: quantum rotational levels cannot exceed l_max + 1 X l_max +1\n");
			return(-1);
		}

		if(system->quantum_rotation_sum <= 0) {
			error("INPUT: quantum rotational sum for partition function invalid\n");
			return(-1);
		} else {
			sprintf(linebuf, "INPUT: Quantum rotation sum = %d\n", system->quantum_rotation_sum);
			output(linebuf);
		}

	}
#endif /* QM_ROTATION */

#ifdef XXX
	if(system->quantum_vibration) output("INPUT: Quantum vibrational eigenspectrum calculation enabled\n");
#endif /* XXX */

	if(system->simulated_annealing) {

		if(system->ensemble != ENSEMBLE_NVT) {
			error("INPUT: Simulated annealing only valid for canonical ensemble\n");
			return(-1);
		} else {
			output("INPUT: Simulated annealing active\n");
		}

		if((system->simulated_annealing_schedule < 0.0) || (system->simulated_annealing_schedule > 1.0)) {
			error("INPUT: invalid simulated annealing temperature schedule specified\n");
			return(-1);
		} else {
			sprintf(linebuf, "INPUT: Simulated annealing temperature schedule = %.3f\n", system->simulated_annealing_schedule);
			output(linebuf);
		}

	}


	if(!system->pdb_input) {
		error("INPUT: must specify an input PDB\n");
		return(-1);
	} else {
		sprintf(linebuf, "INPUT: molecular coordinates are in %s\n", system->pdb_input);
		output(linebuf);
	}


	if(system->cuda) {
#ifndef CUDA
		error("INPUT: cuda keyword enabled, but not compiled into this code\n");
		return(-1);
#else
		output("INPUT: CUDA GPU acceleration activated\n");
#endif /* CUDA */
	}

	if(system->polarization) {

		output("INPUT: Thole polarization activated\n");

		if(system->cuda) {
			if(!system->polar_iterative) {
				error("INPUT: CUDA GPU acceleration available for iterative Thole only\n");
				return(-1);
			} else if(system->damp_type != DAMPING_EXPONENTIAL) {
				error("INPUT: CUDA GPU accleration available for exponential Thole damping only\n");
				return(-1);
			} else if(!system->polar_max_iter) {
				/* XXX disable for 1 iter testing */
				//error("INPUT: Must set polar_max_iter for CUDA GPU acceleration\n");
				//return(-1);
			} else
				output("INPUT: CUDA GPU Thole SCF solver activated\n");
		}

		if(system->polar_iterative && system->polarizability_tensor) {
			error("INPUT: iterative polarizability tensor method not implemented\n");
			return(-1);
		}

		if(!system->polar_iterative && system->polar_zodid) {
			error("INPUT: ZODID and matrix inversion cannot both be set!\n");
			return(-1);
		}

		if(system->damp_type == DAMPING_LINEAR)
			output("INPUT: Thole linear damping activated\n");
		else if(system->damp_type == DAMPING_EXPONENTIAL)
			output("INPUT: Thole exponential damping activated\n");
		else {
			error("INPUT: Thole damping method not specified\n");
			return(-1);
		}

		if(system->polar_damp <= 0.0) {
			error("INPUT: damping factor must be specified\n");
			return(-1);
		} else {
			sprintf(linebuf, "INPUT: Thole damping parameter is %.4f\n", system->polar_damp);
			output(linebuf);

			sprintf(linebuf, "INPUT: Field damping parameter is %.4f\n", system->field_damp);
			output(linebuf);
		}

		if(system->polar_iterative) {

			output("INPUT: Thole iterative solver activated\n");
			if(system->polar_zodid) {
				output("INPUT: ZODID polarization enabled\n");
			}

			if((system->polar_precision > 0.0) && (system->polar_max_iter > 0)) {
				error("INPUT: cannot specify both polar_precision and polar_max_iter, must pick one\n");
				return(-1);
			}

			if(system->polar_precision < 0.0) {
				error("INPUT: invalid polarization iterative precision specified\n");
				return(-1);
			} else if(system->polar_precision > 0.0) {
				sprintf(linebuf, "INPUT: Thole iterative precision is %e A*sqrt(KA) (%e D)\n", system->polar_precision, system->polar_precision/DEBYE2SKA);
				output(linebuf);
			} else {
				sprintf(linebuf, "INPUT: using polar max SCF iterations = %d\n", system->polar_max_iter);
				output(linebuf);
			}

			if(system->polar_sor && system->polar_esor) {
				error("INPUT: cannot specify both SOR and ESOR SCF methods\n");
				return(-1);
			}

			if(system->polar_sor) output("INPUT: SOR SCF scheme active\n");
			else if(system->polar_esor) output("INPUT: ESOR SCF scheme active\n");

			if(system->polar_gamma < 0.0) {
				error("INPUT: invalid Pre-cond/SOR/ESOR gamma set\n");
				return(-1);
			} else {
				sprintf(linebuf, "INPUT: Pre-cond/SOR/ESOR gamma = %.3f\n", system->polar_gamma);
				output(linebuf);
			}


			if(system->polar_gs && system->polar_gs_ranked) {
				error("INPUT: both polar_gs and polar_gs_ranked cannot be set\n");
				return(-1);
			}

			if(system->polar_gs)
				output("INPUT: Gauss-Seidel iteration scheme active\n");
			else if(system->polar_gs_ranked)
				output("INPUT: Gauss-Seidel Ranked iteration scheme active\n");

			if(system->polar_palmo) output("INPUT: Polarization energy of Palmo and Krimm enabled\n");

		} else {
			output("INPUT: Matrix polarization activated\n");
			if(system->polarizability_tensor)
				output("INPUT: Polarizability tensor calculation activated\n");
		}

		if(system->polar_self) output("INPUT: Polarization self-induction is active\n");

	}

	if((system->pbc->volume < 0.0) || (system->pbc->cutoff < 0.0))
		error("INPUT: something is wrong with the periodic boundary conditions\n");
	else {

		sprintf(linebuf, "INPUT: pbc_cutoff set to %.5f A\n", system->pbc->cutoff);
		output(linebuf);
		sprintf(linebuf, "INPUT: basis vector 1 = %.5f %.5f %.5f\n", system->pbc->basis[0][0], system->pbc->basis[0][1], system->pbc->basis[0][2]);
		output(linebuf);
		sprintf(linebuf, "INPUT: basis vector 2 = %.5f %.5f %.5f\n", system->pbc->basis[1][0], system->pbc->basis[1][1], system->pbc->basis[1][2]);
		output(linebuf);
		sprintf(linebuf, "INPUT: basis vector 3 = %.5f %.5f %.5f\n", system->pbc->basis[2][0], system->pbc->basis[2][1], system->pbc->basis[2][2]);
		output(linebuf);
		sprintf(linebuf, "INPUT: unit cell volume = %.3f A^3 (cutoff = %.3f A)\n", system->pbc->volume, system->pbc->cutoff);
		output(linebuf);
		sprintf(linebuf, "INPUT: recip basis vector 1 = %.5f %.5f %.5f\n", system->pbc->reciprocal_basis[0][0], system->pbc->reciprocal_basis[0][1], system->pbc->reciprocal_basis[0][2]);
		output(linebuf);
		sprintf(linebuf, "INPUT: recip basis vector 2 = %.5f %.5f %.5f\n", system->pbc->reciprocal_basis[1][0], system->pbc->reciprocal_basis[1][1], system->pbc->reciprocal_basis[1][2]);
		output(linebuf);
		sprintf(linebuf, "INPUT: recip basis vector 3 = %.5f %.5f %.5f\n", system->pbc->reciprocal_basis[2][0], system->pbc->reciprocal_basis[2][1], system->pbc->reciprocal_basis[2][2]);
		output(linebuf);
	}

	if(system->rd_anharmonic) {

		if(!system->rd_only) {
			error("INPUT: rd_anharmonic being set requires rd_only\n");
			return(-1);
		} else {
			sprintf(linebuf, "INPUT: rd_anharmonic_k = %.3f K/A^2\n", system->rd_anharmonic_k);
			output(linebuf);
			sprintf(linebuf, "INPUT: rd_anharmonic_g = %.3f K/A^4\n", system->rd_anharmonic_g);
			output(linebuf);
		}

	}

	if((system->ensemble != ENSEMBLE_SURF) && (system->ensemble != ENSEMBLE_SURF_FIT)) {

		if(!system->seed) {
			error("INPUT: no seed specified\n");
			return(-1);
		} else {
			sprintf(linebuf, "INPUT: rng seed is %ld\n", system->seed);
			output(linebuf);
		}

		if((system->numsteps < 1) && (system->ensemble != ENSEMBLE_TE) ) {
			error("INPUT: improper numsteps specified\n");
			return(-1);
		} else {
			sprintf(linebuf, "INPUT: each core performing %d simulation steps\n", system->numsteps);
			output(linebuf);
		}
	
		if((system->corrtime < 1) && (system->ensemble != ENSEMBLE_TE) )  {
			error("INPUT: improper corrtime specified\n");
			return(-1);
		} else {
			sprintf(linebuf, "INPUT: system correlation time is %d steps\n", system->corrtime);
			output(linebuf);
		}

		if((system->ensemble != ENSEMBLE_NVE) && (system->ensemble != ENSEMBLE_TE) ) {
			if(system->temperature <= 0.0) {
				error("INPUT: invalid temperature specified\n");
				return(-1);
			} else {
				sprintf(linebuf, "INPUT: system temperature is %.3f K\n", system->temperature);
				output(linebuf);
			}
		} else if ( system->ensemble == ENSEMBLE_NVE ) {
			sprintf(linebuf, "INPUT: NVE energy is %.3f K\n", system->total_energy);
			output(linebuf);
		}

		if(system->free_volume > 0.0) {
			sprintf(linebuf, "INPUT: system free_volume is %.3f A^3\n", system->free_volume);
			output(linebuf);
		}

		if((system->ensemble == ENSEMBLE_UVT) && (system->pressure <= 0.0)) {
			error("INPUT: invalid pressure set for GCMC\n");
			return(-1);
		} else {
			if(system->ensemble == ENSEMBLE_UVT) {
				sprintf(linebuf, "INPUT: reservoir pressure is %.3f atm\n", system->pressure);
				sprintf(linebuf, "INPUT: fugacity is set to %.3f\n", system->fugacity);
				output(linebuf);
			}

			if(system->h2_fugacity) {

				system->fugacity = h2_fugacity(system->temperature, system->pressure);
				if(system->h2_fugacity == 0.0) {
					error("INPUT: error in H2 fugacity assignment\n");
					return(-1);
				}

				sprintf(linebuf, "INPUT: H2 fugacity = %.3f atm\n", system->fugacity);
				output(linebuf);
			}

			if(system->co2_fugacity) {

				system->fugacity = co2_fugacity(system->temperature, system->pressure);
				if(system->co2_fugacity == 0.0) {
					error("INPUT: error in CO2 fugacity assignment\n");
					return(-1);
				}

				sprintf(linebuf, "INPUT: CO2 fugacity = %.3f atm\n", system->fugacity);
				output(linebuf);
			}

                        if(system->ch4_fugacity) {

                                system->fugacity = ch4_fugacity(system->temperature, system->pressure);
                                if(system->ch4_fugacity == 0.0) {
                                        error("INPUT: error in CH4 fugacity assignment\n");
                                        return(-1);
                                }

                                sprintf(linebuf, "INPUT: CH4 fugacity = %.3f atm\n", system->fugacity);
                                output(linebuf);
                        }

                        if(system->n2_fugacity) {

                                system->fugacity = n2_fugacity(system->temperature, system->pressure);
                                if(system->n2_fugacity == 0.0) {
                                        error("INPUT: error in N2 fugacity assignment\n");
                                        return(-1);
                                }

                                sprintf(linebuf, "INPUT: N2 fugacity = %.3f atm\n", system->fugacity);
                                output(linebuf);
                        }
		}

		if ( system->ensemble != ENSEMBLE_TE ) {

			sprintf(linebuf, "INPUT: insert probability is %.3f\n", system->insert_probability);
			output(linebuf);

			sprintf(linebuf, "INPUT: move probability is %.3f\n", system->move_probability);
			output(linebuf);

			sprintf(linebuf, "INPUT: gwp probability is %.3f\n", system->gwp_probability);
			output(linebuf);

			sprintf(linebuf, "INPUT: rotation probability is %.3f\n", system->rot_probability);
			output(linebuf);

			sprintf(linebuf, "INPUT: spinflip probability is %.3f\n", system->spinflip_probability);
			output(linebuf);

		}

		/* autoreject insertions closer than some scaling factor of sigma */
		if(system->cavity_autoreject) {

			output("INPUT: cavity autorejection activated\n");

			if((system->cavity_autoreject_scale <= 0.0) || (system->cavity_autoreject_scale > 1.0))
				error("INPUT: cavity_autoreject_scale either not set or out of range\n");

		}

		if(system->cavity_bias) {

			if((system->cavity_grid_size <= 0) || (system->cavity_radius <= 0.0)) {
				error("INPUT: invalid cavity grid or radius specified\n");
			} else {
				output("INPUT: cavity-biased umbrella sampling activated\n");
				sprintf(linebuf, "INPUT: cavity grid size is %dx%dx%d points with a sphere radius of %.3f A\n", system->cavity_grid_size, system->cavity_grid_size, system->cavity_grid_size, system->cavity_radius);
				output(linebuf);
			}

		}

		if(!system->pdb_output) {
			error("INPUT: must specify an output PDB\n");
			return(-1);
		} else {
			sprintf(linebuf, "INPUT: will be writing final configuration to %s\n", system->pdb_output);
			output(linebuf);
		}

		if(!system->pdb_restart) {
			error("INPUT: must specify a restart PDB\n");
			return(-1);
		} else {
			sprintf(linebuf, "INPUT: will be writing restart configuration to %s\n", system->pdb_restart);
			output(linebuf);
		}

		if(!system->traj_output) {
			error("INPUT: must specify a trajectory PDB\n");
			return(-1);
		} else {
			sprintf(linebuf, "INPUT: will be writing trajectory to %s\n", system->traj_output);
			output(linebuf);
		}

		if(system->insert_input) {
			sprintf( linebuf, "INPUT: inserted molecules will be selected from: %s\n", system->insert_input );
			output( linebuf );
		} else {
			output( "INPUT: inserted molecules will be selected from the restart file.\n" );
		}

		if(system->polarization && (!system->dipole_output || !system->field_output)) {
			error("INPUT: must specify a dipole and field output file\n");
			return(-1);
		} else if(system->polarization) {
			sprintf(linebuf, "INPUT: dipole field will be written to %s\n", system->dipole_output);
			output(linebuf);
			sprintf(linebuf, "INPUT: electric field will be written to %s\n", system->field_output);
			output(linebuf);
		}

		if(system->wpi) {
			output("INPUT: Widom Particle Insertion is enabled\n");
		}

	}
		
	if(system->calc_hist){

		output("INPUT: histogram calculation will be performed\n");
		if(!system->hist_resolution){
			output("INPUT: no histogram resolution set but histogram calculation requested\n");
			output("INPUT: setting hist_resolution to default value of 0.7A\n");
			system->hist_resolution=0.7;
		}
		else if(system->hist_resolution<0.01 || system->hist_resolution>5.0){
			output("INPUT: histogram resolution out of bounds\n");
			output("INPUT: setting hist_resolution to default value of 0.7A\n");
			system->hist_resolution=0.7;
		}
		else if(!system->histogram_output){
			output("INPUT: no histogram outputfile selected, defaulting to histogram.dat\n");
			system->histogram_output=calloc(MAXLINE,sizeof(char));
			memnullcheck(system->histogram_output,MAXLINE*sizeof(char),103);
			sprintf(system->histogram_output,"histogram.dat");
		}
		else{
			sprintf(linebuf,"INPUT: histogram resolution set to %.3f A\n",system->hist_resolution);
			output(linebuf);
		}


		if(system->max_bondlength < .5){
			output("INPUT: max_bondlength either not set or out of bounds\n");
			output("INPUT: setting max_bondlength to default value of 1.8A\n");
			system->max_bondlength=1.8;
		}	

		if(!system->frozen_output){
			output("INPUT: no frozen_output set! setting frozen coordinate output file to frozen.dx\n");
			system->frozen_output = calloc(MAXLINE, sizeof(char));
			memnullcheck(system->frozen_output,MAXLINE*sizeof(char),104);
			sprintf(system->frozen_output,"frozen.dx");
		} else {

			sprintf(linebuf, "INPUT: will be writing frozen coordinates to %s\n", system->frozen_output);
			output(linebuf);

		}


	}


	return(0);

}


molecule_t *read_molecules(system_t *system) {

	int i;
	molecule_t *molecules, *molecule_ptr;
	atom_t     *atom_ptr, *prev_atom_ptr;
	char       linebuf[MAXLINE], *n;
	FILE       *fp;
	char       token_atom[MAXLINE],
                   token_atomid[MAXLINE],
                   token_atomtype[MAXLINE],
                   token_moleculetype[MAXLINE];
	char       token_frozen[MAXLINE],
                   token_moleculeid[MAXLINE],
                   token_x[MAXLINE],
                   token_y[MAXLINE],
                   token_z[MAXLINE];
	char token_mass[MAXLINE], token_charge[MAXLINE], token_alpha[MAXLINE], token_epsilon[MAXLINE], token_sigma[MAXLINE], token_omega[MAXLINE], token_gwp_alpha[MAXLINE];
	int current_frozen, current_adiabatic, current_spectre, current_target;
	int current_moleculeid, current_atomid;
	double current_x, current_y, current_z;
	double current_mass, current_charge, current_alpha, current_epsilon, current_sigma, current_omega, current_gwp_alpha;
	double current_molecule_mass;
        int current_site_neighbor;
	int moveable, spectres, targets;
	int atom_counter;

	/* allocate the start of the list */
	molecules = calloc(1, sizeof(molecule_t));
	memnullcheck(molecules,sizeof(molecule_t),105);
	molecule_ptr = molecules;
	molecule_ptr->id = 1;
	molecule_ptr->atoms = calloc(1, sizeof(atom_t));
	memnullcheck(molecule_ptr->atoms,sizeof(atom_t),106);
	atom_ptr = molecule_ptr->atoms;
	prev_atom_ptr = atom_ptr;


	/* open the molecule input file */
	fp = fopen(system->pdb_input, "r");
	if(!fp) {
		sprintf(linebuf, "INPUT: couldn't open PDB input file %s\n", system->pdb_input);
		error(linebuf);	
		return(NULL);
	}

	/* clear the linebuffer and read the tokens in */
	atom_counter = 0;
	memset(linebuf, 0, MAXLINE);
	n = fgets(linebuf, MAXLINE, fp);
	while(n) {

		/* clear the tokens */
		memset( token_atom,         0, MAXLINE);
		memset( token_atomid,       0, MAXLINE);
		memset( token_atomtype,     0, MAXLINE);
		memset( token_moleculetype, 0, MAXLINE);
		memset( token_frozen,       0, MAXLINE);
		memset( token_moleculeid,   0, MAXLINE);
		memset( token_x,            0, MAXLINE);
		memset( token_y,            0, MAXLINE);
		memset( token_z,            0, MAXLINE);
		memset( token_mass,         0, MAXLINE);
		memset( token_charge,       0, MAXLINE);
		memset( token_alpha,        0, MAXLINE);
		memset( token_epsilon,      0, MAXLINE);
		memset( token_sigma,        0, MAXLINE);
		memset(token_omega, 0, MAXLINE);
		memset(token_gwp_alpha, 0, MAXLINE);

		/* parse the line */
		sscanf(linebuf, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n", token_atom, token_atomid, token_atomtype, token_moleculetype, token_frozen, token_moleculeid, token_x, token_y, token_z, token_mass, token_charge, token_alpha, token_epsilon, token_sigma, token_omega, token_gwp_alpha);

		if(!strcasecmp(token_atom, "ATOM") && strcasecmp(token_moleculetype, "BOX")) {

			current_frozen = 0; current_adiabatic = 0; current_spectre = 0; current_target = 0;
			if(!strcasecmp(token_frozen, "F"))
				current_frozen = 1;
			if(!strcasecmp(token_frozen, "A"))
				current_adiabatic = 1;
			if(!strcasecmp(token_frozen, "S"))
				current_spectre = 1;
			if(!strcasecmp(token_frozen, "T"))
				current_target = 1;

			current_moleculeid    = atoi(token_moleculeid);
			current_atomid        = atoi(token_atomid);
			current_x             = atof(token_x);
			current_y             = atof(token_y);
			current_z             = atof(token_z);
			current_mass          = atof(token_mass);	/* mass in amu */
			current_charge        = atof(token_charge);
			current_charge       *= E2REDUCED;		/* convert charge into reduced units  */
			current_alpha         = atof(token_alpha);
			current_epsilon       = atof(token_epsilon);
			current_sigma         = atof(token_sigma);
			current_omega         = atof(token_omega);
			current_gwp_alpha     = atof(token_gwp_alpha);;
			// Functionality of site_neighbor disabled in favor of omega/gwp_alpha parameters
			// Current behavior is to default to atom 0, typically the center of mass for
			// the molecule.
			current_site_neighbor = 0; //atoi( token_site_neighbor );
                        
			if(current_frozen)
                            current_charge *= system->scale_charge;

			if(molecule_ptr->id != current_moleculeid) {
				molecule_ptr->next = calloc(1, sizeof(molecule_t));
				memnullcheck(molecule_ptr,sizeof(molecule_t),107);
				molecule_ptr = molecule_ptr->next;
				molecule_ptr->atoms = calloc(1, sizeof(atom_t));
				memnullcheck(molecule_ptr->atoms,sizeof(atom_t),108);
				prev_atom_ptr->next = NULL;
				free(atom_ptr);
				atom_ptr = molecule_ptr->atoms;
			}

                        strcpy(molecule_ptr->moleculetype, token_moleculetype);

			molecule_ptr->id        = current_moleculeid;
			molecule_ptr->frozen    = current_frozen;
			molecule_ptr->adiabatic = current_adiabatic;
			molecule_ptr->spectre   = current_spectre;
			molecule_ptr->target    = current_target;
			molecule_ptr->mass     += current_mass;

#ifdef QM_ROTATION
			/* if quantum rot calc. enabled, allocate the necessary structures */
			if(system->quantum_rotation && !molecule_ptr->frozen) {

				molecule_ptr->quantum_rotational_energies = calloc(system->quantum_rotation_level_max, sizeof(double));
				memnullcheck(molecule_ptr->quantum_rotational_energies,system->quantum_rotation_level_max*sizeof(double),109);
				molecule_ptr->quantum_rotational_eigenvectors = calloc(system->quantum_rotation_level_max, sizeof(complex_t *));
				memnullcheck(molecule_ptr->quantum_rotational_eigenvectors,system->quantum_rotation_level_max*sizeof(complex_t *),110);
				for(i = 0; i < system->quantum_rotation_level_max; i++){
					molecule_ptr->quantum_rotational_eigenvectors[i] = 
						calloc((system->quantum_rotation_l_max + 1)*(system->quantum_rotation_l_max + 1), sizeof(complex_t));
					memnullcheck(molecule_ptr->quantum_rotational_eigenvectors[i],
						(system->quantum_rotation_l_max+1)*(system->quantum_rotation_l_max+1)*sizeof(complex_t),111);
				}
				molecule_ptr->quantum_rotational_eigensymmetry = calloc(system->quantum_rotation_level_max, sizeof(int));
				memnullcheck(molecule_ptr->quantum_rotational_eigensymmetry,system->quantum_rotation_level_max*sizeof(int),112);

			}
#endif /* QM_ROTATION */

#ifdef XXX
			/* if quantum vib calc. enabled, allocate the necessary structures */
			if(system->quantum_vibration && !molecule_ptr->frozen) {

				molecule_ptr->quantum_vibrational_energies = calloc(system->quantum_vibration_level_max, sizeof(double));
				memnullcheck(molecule_ptr->quantum_vibrational_energies,system->quantum_vibration_level_max*sizeof(double),113);
				molecule_ptr->quantum_vibrational_eigenvectors = calloc(system->quantum_vibration_level_max, sizeof(complex_t *));
				memnullcheck(molecule_ptr->quantum_vibrational_eigenvectors,system->quantum_vibration_level_max*sizeof(complex_t *));
				for(i = 0; i < system->quantum_vibration_level_max; i++) {
					molecule_ptr->quantum_vibrational_eigenvectors[i] = 
						calloc((system->quantum_vibration_l_max+1)*(system->quantum_vibration_l_max+1), sizeof(complex_t));
					memnullcheck(molecule_ptr->quantum_vibrational_eigenvectors[i],
						(system->quantum_vibration_l_max+1)*(system->quantum_vibration_l_max+1)*sizeof(complex_t),114);
				}
				molecule_ptr->quantum_vibrational_eigensymmetry = calloc(system->quantum_vibration_level_max, sizeof(int));
				memnullcheck(molecule_ptr->quantum_vibrational_eigensymmetry, system->quantum_vibration_level_max*sizeof(int),115);

			}
#endif /* XXX */

			++atom_counter;
			atom_ptr->id        = atom_counter;
                        atom_ptr->bond_id   = current_atomid;
			memset(atom_ptr->atomtype, 0, MAXLINE);
			strcpy(atom_ptr->atomtype, token_atomtype);
			atom_ptr->frozen    = current_frozen;
			atom_ptr->adiabatic = current_adiabatic;
			atom_ptr->spectre = current_spectre;
			atom_ptr->target = current_target;
			atom_ptr->pos[0] = current_x;
			atom_ptr->pos[1] = current_y;
			atom_ptr->pos[2] = current_z;
			atom_ptr->mass = current_mass;
			atom_ptr->charge = current_charge;
			atom_ptr->polarizability = current_alpha;
			atom_ptr->epsilon = current_epsilon;
			atom_ptr->sigma = current_sigma;
			//if polarvdw is on, we want sigma = 1
			//it's an unneeded parameter, and complicates the mixing rules
			if ( system->polarvdw && current_sigma != 1 ) {
				fprintf(stderr,"INPUT: WARNING POLARVDW REQUIRES SIGMA = 1.\n"
					"INPUT: INPUT VALUE OF %lf IGNORED.\n", current_sigma);
				current_sigma = 1;
				atom_ptr->sigma = 1;
			}
			atom_ptr->omega = current_omega;
			atom_ptr->gwp_alpha = current_gwp_alpha;
			if(current_gwp_alpha != 0.)
				atom_ptr->gwp_spin = 1;
			else
				atom_ptr->gwp_spin = 0;

			atom_ptr->site_neighbor_id = current_site_neighbor;
			atom_ptr->next = calloc(1, sizeof(atom_t));
			memnullcheck(atom_ptr->next,sizeof(atom_t),116);
			prev_atom_ptr  = atom_ptr;
			atom_ptr       = atom_ptr->next;
		}

		memset(linebuf, 0, MAXLINE);
		n = fgets(linebuf, MAXLINE, fp);
	}

	/* terminate the atom list */
	prev_atom_ptr->next = NULL;
	free(atom_ptr);

	/* scan the list, make sure that there is at least one moveable molecule */
	moveable = 0;
	spectres = 0;
	targets = 0;
	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {

		if(!molecule_ptr->frozen ) ++moveable;
		if(molecule_ptr ->target ) ++targets;
		if(molecule_ptr ->spectre) ++spectres;

	}

	if(system->spectre) {

		if(!spectres || !targets) {
			error("INPUT: either no targets or spectres found\n");
			return(NULL);
		}

	} else {

		if(!moveable) {
			error("INPUT: no moveable molecules found, there must be at least one in your PDB file\n");
			return(NULL);
		}

	}

	//you forgot to close the file!
	fclose(fp);

	return(molecules);

}



//  read_insertion_molecules( system_t * ) was a cut/paste of read_molecules
//  modified to read the candidate insertion molecules from a separate file. 
///////////////////////////////////////////////////////////////////////////////

molecule_t *read_insertion_molecules(system_t *system) {

	int i;

	molecule_t *molecules, 
	           *molecule_ptr;

	atom_t     *atom_ptr, 
	           *prev_atom_ptr;

	char       linebuf[MAXLINE], *n;

	FILE       *fp;
	char       token_atom[MAXLINE], token_atomid[MAXLINE], token_atomtype[MAXLINE],
	           token_moleculetype[MAXLINE],
	           token_frozen[MAXLINE],
	           token_moleculeid[MAXLINE],
	           token_x[MAXLINE], token_y[MAXLINE], token_z[MAXLINE],
	           token_mass[MAXLINE],
		   token_charge[MAXLINE],
	           token_alpha[MAXLINE], token_epsilon[MAXLINE], token_sigma[MAXLINE], 
	           token_omega[MAXLINE], token_gwp_alpha[MAXLINE];
	
	int        current_frozen, 
	           current_adiabatic,
	           current_spectre, 
	           current_target,
	           current_moleculeid, 
	           current_atomid,
	           current_site_neighbor;
	double     current_x, current_y, current_z,
	           current_mass,  current_charge, 
                   current_alpha, current_epsilon, current_sigma, current_omega, current_gwp_alpha, 
	           current_molecule_mass;

	int        moveable, 
	           spectres, 
	           targets,

	           atom_counter;


	// allocate the start of the list 
	molecules           = calloc(1, sizeof(molecule_t));
	memnullcheck(molecules,sizeof(molecule_t),117);
	molecule_ptr        = molecules;
	molecule_ptr->id    = 1;
	molecule_ptr->atoms = calloc(1, sizeof(atom_t));
	memnullcheck(molecule_ptr->atoms,sizeof(atom_t),118);
	atom_ptr            = molecule_ptr->atoms;
	prev_atom_ptr       = atom_ptr;


	// open the molecule input file 
	fp = fopen(system->insert_input, "r");
	if(!fp) {
		sprintf(linebuf, "INPUT: couldn't open insertion input file %s\n", system->insert_input);
		error(linebuf);	
		return(NULL);
	}

	// clear the linebuffer and read the tokens in 
	atom_counter = 0;
	memset(linebuf, 0, MAXLINE);
	n = fgets(linebuf, MAXLINE, fp);
	while(n) {

		// clear the tokens 
		memset( token_atom,         0, MAXLINE);
		memset( token_atomid,       0, MAXLINE);
		memset( token_atomtype,     0, MAXLINE);
		memset( token_moleculetype, 0, MAXLINE);
		memset( token_frozen,       0, MAXLINE);
		memset( token_moleculeid,   0, MAXLINE);
		memset( token_x,            0, MAXLINE);
		memset( token_y,            0, MAXLINE);
		memset( token_z,            0, MAXLINE);
		memset( token_mass,         0, MAXLINE);
		memset( token_charge,       0, MAXLINE);
		memset( token_alpha,        0, MAXLINE);
		memset( token_epsilon,      0, MAXLINE);
		memset( token_sigma,        0, MAXLINE);
		memset( token_omega,        0, MAXLINE);
		memset( token_gwp_alpha,    0, MAXLINE);

		/* parse the line */
		sscanf(linebuf, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n", 
		       token_atom, token_atomid, token_atomtype, 
		       token_moleculetype, 
		       token_frozen, 
		       token_moleculeid, 
		       token_x, token_y, token_z, 
		       token_mass, token_charge, 
		       token_alpha, token_epsilon, token_sigma, token_omega, token_gwp_alpha
		);

		if(!strcasecmp(token_atom, "ATOM") && strcasecmp(token_moleculetype, "BOX")) {

			current_frozen = 0; current_adiabatic = 0; current_spectre = 0; current_target = 0;
			if(!strcasecmp(token_frozen, "F"))
				current_frozen = 1;
			if(!strcasecmp(token_frozen, "A"))
				current_adiabatic = 1;
			if(!strcasecmp(token_frozen, "S"))
				current_spectre = 1;
			if(!strcasecmp(token_frozen, "T"))
				current_target = 1;

			current_moleculeid    = atoi(token_moleculeid);
			current_atomid        = atoi(token_atomid);
			current_x             = atof(token_x);
			current_y             = atof(token_y);
			current_z             = atof(token_z);
			current_mass          = atof(token_mass);	// mass in amu 
			current_charge        = atof(token_charge);
			current_charge       *= E2REDUCED;		// convert charge into reduced units  
			current_alpha         = atof(token_alpha);
			current_epsilon       = atof(token_epsilon);
			current_sigma         = atof(token_sigma);
			current_omega         = atof(token_omega);
			current_gwp_alpha     = atof(token_gwp_alpha);
			// Functionality of site_neighbor disabled in favor of omega/gwp_alpha parameters
			// Current behavior is to default to atom 0, typically the center of mass for
			// the molecule.
			current_site_neighbor = 0;

                        if(current_frozen)
                            current_charge *= system->scale_charge;

			if(molecule_ptr->id != current_moleculeid) {
				molecule_ptr->next = calloc(1, sizeof(molecule_t));
				memnullcheck(molecule_ptr->next,sizeof(molecule_t),119);
				molecule_ptr = molecule_ptr->next;
				molecule_ptr->atoms = calloc(1, sizeof(atom_t));
				memnullcheck(molecule_ptr->atoms,sizeof(atom_t),120);
				prev_atom_ptr->next = NULL;
				free(atom_ptr);
				atom_ptr = molecule_ptr->atoms;
			}

                        strcpy(molecule_ptr->moleculetype, token_moleculetype);

			molecule_ptr->id        = current_moleculeid;
			molecule_ptr->frozen    = current_frozen;
			molecule_ptr->adiabatic = current_adiabatic;
			molecule_ptr->spectre   = current_spectre;
			molecule_ptr->target    = current_target;
			molecule_ptr->mass     += current_mass;

#ifdef QM_ROTATION
			/* if quantum rot calc. enabled, allocate the necessary structures */
			if(system->quantum_rotation && !molecule_ptr->frozen) {

				molecule_ptr->quantum_rotational_energies = calloc(system->quantum_rotation_level_max, sizeof(double));
				memnullcheck(molecule_ptr->quantum_rotational_energies, system->quantum_rotation_level_max*sizeof(double),121);
				molecule_ptr->quantum_rotational_eigenvectors = calloc(system->quantum_rotation_level_max, sizeof(complex_t *));
				memnullcheck(molecule_ptr->quantum_rotational_eigenvectors, system->quantum_rotation_level_max*sizeof(complex_t *), 122);
				for(i = 0; i < system->quantum_rotation_level_max; i++) {
					molecule_ptr->quantum_rotational_eigenvectors[i] = 
						calloc((system->quantum_rotation_l_max + 1)*(system->quantum_rotation_l_max + 1), sizeof(complex_t));
					memnullcheck(molecule_ptr->quantum_rotational_eigenvectors[i],
						(system->quantum_rotation_l_max + 1)*(system->quantum_rotation_l_max + 1)*sizeof(complex_t),123);
				}
				molecule_ptr->quantum_rotational_eigensymmetry = calloc(system->quantum_rotation_level_max, sizeof(int));
				memnullcheck(molecule_ptr->quantum_rotational_eigensymmetry,system->quantum_rotation_level_max*sizeof(int),124);

			}
#endif /* QM_ROTATION */

#ifdef XXX
			/* if quantum vib calc. enabled, allocate the necessary structures */
			if(system->quantum_vibration && !molecule_ptr->frozen) {

				molecule_ptr->quantum_vibrational_energies = calloc(system->quantum_vibration_level_max, sizeof(double));
				memnullcheck(molecule_ptr->quantum_vibrational_energies, system->quantum_vibration_level_max*sizeof(double), 125);
				molecule_ptr->quantum_vibrational_eigenvectors = calloc(system->quantum_vibration_level_max, sizeof(complex_t *));
				memnullcheck(molecule_ptr->quantum_vibrational_eigenvectors, system->quantum_vibration_level_max*sizeof(complex_t *), 126);
				for(i = 0; i < system->quantum_vibration_level_max; i++) {
					molecule_ptr->quantum_vibrational_eigenvectors[i] = 
						calloc((system->quantum_vibration_l_max + 1)*(system->quantum_vibration_l_max + 1), sizeof(complex_t));
					memnullcheck(molecule_ptr->quantum_vibrational_eigenvectors[i],
						(system->quantum_vibration_l_max + 1)*(system->quantum_vibration_l_max + 1)*sizeof(complex_t), 127);
				}
				molecule_ptr->quantum_vibrational_eigensymmetry = calloc(system->quantum_vibration_level_max, sizeof(int));
				memnullcheck(molecule_ptr->quantum_vibrational_eigensymmetry, system->quantum_vibration_level_max*sizeof(int), 128);

			}
#endif /* XXX */

			++atom_counter;
			atom_ptr->id        = atom_counter;
                        atom_ptr->bond_id   = current_atomid;
			memset(atom_ptr->atomtype, 0, MAXLINE);
			strcpy(atom_ptr->atomtype, token_atomtype);
			atom_ptr->frozen    = current_frozen;
			atom_ptr->adiabatic = current_adiabatic;
			atom_ptr->spectre   = current_spectre;
			atom_ptr->target    = current_target;
			atom_ptr->pos[0]    = current_x;
			atom_ptr->pos[1]    = current_y;
			atom_ptr->pos[2]    = current_z;
			atom_ptr->mass      = current_mass;
			atom_ptr->charge    = current_charge;
			atom_ptr->polarizability = current_alpha;
			atom_ptr->epsilon   = current_epsilon;
			atom_ptr->sigma     = current_sigma;
			//if polarvdw is on, we want sigma = 1
			//it's an unneeded parameter, and complicates the mixing rules
			if ( system->polarvdw && current_sigma != 1 ) {
				fprintf(stderr,"INPUT: WARNING POLARVDW REQUIRES SIGMA = 1.\n"
					"INPUT: INPUT VALUE OF %lf IGNORED.\n", current_sigma);
				current_sigma   = 1;
				atom_ptr->sigma = 1;
			}
			atom_ptr->omega     = current_omega;
			atom_ptr->gwp_alpha = current_gwp_alpha;
			if(current_gwp_alpha != 0.)
				atom_ptr->gwp_spin = 1;
			else
				atom_ptr->gwp_spin = 0;

			atom_ptr->site_neighbor_id = current_site_neighbor;
			atom_ptr->next = calloc(1, sizeof(atom_t));
			memnullcheck(atom_ptr->next,sizeof(atom_t),129);
			prev_atom_ptr  = atom_ptr;
			atom_ptr       = atom_ptr->next;
		}

		memset(linebuf, 0, MAXLINE);
		n = fgets(linebuf, MAXLINE, fp);
	}

	// terminate the atom list 
	prev_atom_ptr->next = NULL;
	free(atom_ptr);



	// Count the molecules and create an array of pointers, where each pointer
	// points directly to a molecule in the linked list. 
	/////////////////////////////////////////////////////////////////////////////

	// count	
	int molecule_counter = 0;
	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) 
		molecule_counter++;

	// allocate space for array
	system->insertion_molecules_array = malloc( molecule_counter * sizeof(molecule_t*) );
	if(!system->insertion_molecules_array ) {
		error( "INPUT: ERROR - exhausted memory while allocating insertion molecules array.\n" );
		return NULL;
	}

	// point array pointers to their corresponding molecules
	molecule_counter=0;
	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) { 
		system->insertion_molecules_array[ molecule_counter ] = molecule_ptr;
		molecule_counter++;
	}
	system->num_insertion_molecules = molecule_counter;


	  //  ERROR CHECKING OF SOME SORT MAY NEED TO GO HERE
	 //   WHAT'S ALLOWED/DISALLOWED FOR INSERTED MOLECULES?
	/////////////////////////////////////////////////////////

	
	return(molecules);
}



system_t *setup_system(char *input_file) {

	system_t *system;
	char linebuf[MAXLINE];

	/* read in all of the tokens and parameters from the config file */
	system = read_config(input_file);
	if(!system) {
		error("INPUT: error reading config file\n");
		return(NULL);
	} else
		output("INPUT: finished reading config file\n");

	/* validate configuration parameters */
	if(check_system(system) < 0) {
		error("INPUT: invalid config parameters specified\n");
		return(NULL);
	} else
		output("INPUT: config file validated\n");

	/* read in the input pdb and setup the data structures */
	system->molecules = read_molecules(system);
	if(!system->molecules) {
		error("INPUT: error reading in molecules\n");
		return(NULL);
	} else
		output("INPUT: finished reading in molecules\n");

	// Read in the insertion molecules
	if( system->insert_input ) {
		system->insertion_molecules = read_insertion_molecules(system);
		if(!system->insertion_molecules) {
			error("INPUT: error reading in insertion molecules\n");
			return(NULL);
		} else
			output("INPUT: finished reading in insertion molecules\n");
	}

	/* allocate the necessary pairs */
	setup_pairs(system->molecules);
	output("INPUT: finished allocating pair lists\n");

	/* calculate the periodic boundary conditions */
	pbc(system->pbc);
	output("INPUT: finished calculating periodic boundary conditions\n");

	/* get all of the pairwise interactions, exclusions, etc. */
	if(system->cavity_bias) setup_cavity_grid(system);
	pairs(system);

	/* set all pairs to initially have their energies calculated */
	flag_all_pairs(system);
	output("INPUT: finished calculating pairwise interactions\n");

	/* set the ewald gaussian width appropriately */
	if(system->ewald_alpha == EWALD_ALPHA)
		system->ewald_alpha = 3.5/system->pbc->cutoff;
	if(!(system->sg || system->rd_only)) {
		sprintf(linebuf, "INPUT: Ewald gaussian width = %f A\n", system->ewald_alpha);
		output(linebuf);
		sprintf(linebuf, "INPUT: Ewald kmax = %d\n", system->ewald_kmax);
		output(linebuf);
	}

	/* write the frozen lattice configuration */
	if(system->file_pointers.fp_frozen) write_frozen(system->file_pointers.fp_frozen,system);

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
        if(!fp_fit)
        {
            sprintf( errMsg, "INPUT: Unable to open input file: %s\n", filename );
            error( errMsg );
            return (NULL);
        }





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
                while(   (nodePtr->next)   &&   ((nodePtr->next->r) < r)   )
                    // Stop list-traversal if there is no next item, or if next r-value is >= to the current
                    // r-value (consequently maintaining a sorted list).
                    // Recall that for the && operator, the 2nd condition will not be evaluated if the 1st
                    // fails, i.e., we will not check the next r-value if there is no next r-value...
                    nodePtr = nodePtr->next;

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
            sprintf( errMsg, "INPUT: Exhausted memory transferring curve data to array, for curve: %s\n\0", filename );
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
        for( currentCurve=1; currentCurve<nCurves; ++currentCurve )
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






#ifdef DEBUG
void test_list(molecule_t *molecules) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	pair_t *pair_ptr;

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {

		printf("moleculeid = %d\n", molecule_ptr->id);
		printf("moleculetype = %s\n", molecule_ptr->moleculetype);
		printf("molecule_frozen = %d\n", molecule_ptr->frozen);
		printf("molecule_mass = %f\n", molecule_ptr->mass);
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			printf("atomtype = %s x = %f y = %f z = %f\n", atom_ptr->atomtype, atom_ptr->pos[0], atom_ptr->pos[1], atom_ptr->pos[2]);
			printf("atom frozen = %d mass = %f, charge = %f, alpha = %f, eps = %f, sig = %f\n", atom_ptr->frozen, atom_ptr->mass, atom_ptr->charge, atom_ptr->polarizability, atom_ptr->epsilon, atom_ptr->sigma);
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next)
				if(!(pair_ptr->rd_excluded || pair_ptr->es_excluded || pair_ptr->frozen)) printf("pair = 0x%lx eps = %f sig = %f\n", pair_ptr, pair_ptr->epsilon, pair_ptr->sigma);
		}

	}

	fflush(stdout);

}

void test_molecule(molecule_t *molecule) {

	atom_t *atom_ptr;
	pair_t *pair_ptr;

	printf("moleculeid = %d\n", molecule->id);
	printf("moleculetype = %s\n", molecule->moleculetype);
	printf("molecule_frozen = %d\n", molecule->frozen);
	printf("molecule_mass = %f\n", molecule->mass);
	for(atom_ptr = molecule->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
		printf("atomtype = %s x = %f y = %f z = %f\n", atom_ptr->atomtype, atom_ptr->pos[0], atom_ptr->pos[1], atom_ptr->pos[2]);
		printf("atom frozen = %d mass = %f, charge = %f, alpha = %f, eps = %f, sig = %f\n", atom_ptr->frozen, atom_ptr->mass, atom_ptr->charge, atom_ptr->polarizability, atom_ptr->epsilon, atom_ptr->sigma);
		for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
			printf("pair at 0x%lx\n", pair_ptr);fflush(stdout);
		}
	}

printf("...finished\n");fflush(stdout);


}
#endif /* DEBUG */


