/* Includes subroutines used by surface_fit in surface.c */

#include <stdio.h>
#include <stdlib.h>
#include <mc.h>
#include <surf_const.h>

void surface_curve(system_t *system, double r_min, double r_max, double r_inc, double *curve) {
	int i;
	double r;

	for(r = r_min, i = 0; r <= r_max; r += r_inc, i++) {
		surface_dimer_geometry(system, r, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		curve[i] = surface_energy(system, ENERGY_TOTAL);
	}

	return;
}

double min ( double a, double b ) {
	if ( a > b ) return b;
	else return a;
}

// calculates the current error for a surface energy versus the ab-initio curve
double error_calc ( system_t * system, int nCurves, int nPoints, curveData_t * curve, double max_energy ) {
	int i, j;
	double current_error = 0;
	double buffer, weight;
	double a, b;

	// adjust weight constant in the weighting exponential function for calculating sq_error
	double kweight = ((system->surf_weight_constant_on) ? system->surf_weight_constant : WEIGHT_CONSTANT );

	for( j=0; j<nCurves; j++ ) {
			for(i = 0; i < nPoints; i++) {

//changed so all points contribute until both curves (ab-init and fit) exceed max_energy
				a=min(max_energy,curve[j].input[i]);
				b=min(max_energy,curve[j].output[i]);
				weight = exp(kweight*(max_energy-a)/max_energy);
				buffer = a-b;
				buffer *= buffer;
				current_error += curve[j].weight * weight * buffer;

/*					if(curve[j].input[i] < max_energy) {
							weight = exp(kweight*(max_energy-curve[j].input[i])/max_energy);
							buffer = curve[j].input[i] - curve[j].output[i];
							buffer *= buffer;
							current_error += curve[j].weight * weight * buffer;
					}*/
			}
	}

	return current_error;
}

//allocate curves for surface_fit
int alloc_curves ( int nCurves, int nPoints, curveData_t * curve ) {
	int i, memSuccessFlag;

  // This flag will be AND'ed with all memory allocation results
  // so that any 1 failure will set it to false:
	memSuccessFlag=1;

	for ( i=0; i<nCurves; i++ ) {
      // allocate space for the 2d output array
      memSuccessFlag =
      (curve[i].output = calloc(nPoints, sizeof(double)))  &&  memSuccessFlag;
      // allocate space for the 2d global array
      memSuccessFlag =
      (curve[i].global = calloc(nPoints, sizeof(double)))  &&  memSuccessFlag;
	}

	if(!memSuccessFlag) {
		error( "SURFACE: Exhausted memory while allocating input/output/global arrays.\n" );

		// Upon memory allocation failure, return any memory
		// that was actually allocated back to the system.

		for( i=0; i<nCurves; i++ ) {
				if(curve[i].input)
						free(curve[i].input);
				if(curve[i].output)
						free(curve[i].output);
				if(curve[i].global)
						free(curve[i].global);
		}
		free( curve );
		return (-1);
	}

	return 0;
}

//output dimer configurations to pqr for verification
void output_pqrs ( system_t * system, int nCurves, curveData_t * curve ) {

	int i;
  char filename[500];

  for( i=0; i<nCurves; i++ ) {

      // Apply rotation associated with i-th curve and separate COM by 5 Angstroms
      surface_dimer_geometry(system, 5.0, curve[i].alpha1, curve[i].beta1, curve[i].gamma1,
                                          curve[i].alpha2, curve[i].beta2, curve[i].gamma2 );
      sprintf( filename, "%s.pqr\0", curve[i].id );
      write_molecules( system, filename );

      // Restore the molecule to its initial orientation
      surface_dimer_geometry(system, 0.0, -curve[i].gamma1, -curve[i].beta1, -curve[i].alpha1,
                                          -curve[i].gamma2, -curve[i].beta2, -curve[i].alpha2 );

  }

	return;

}

void output_params ( double temperature, double current_error, param_t * params ) {
	param_t * param_ptr;
	printf("temperature = %f, sq_error = %f\n", temperature, current_error);
	for(param_ptr = params; param_ptr; param_ptr = param_ptr->next) {
    printf("\tatomtype %s q = %f eps = %f sig = %f omega = %f dr = %f\n",
    param_ptr->atomtype,
    param_ptr->charge/E2REDUCED,
    param_ptr->epsilon,
    param_ptr->sigma,
    param_ptr->omega,
    param_ptr->dr);
  }
  fflush(stdout);

	return;
}

//record parameters. we will determine which ones need to be adjusted later
param_t * record_params ( system_t * system ) {

	char last_atomtype[MAXLINE];
	last_atomtype[0] = '\0'; //initialize last_atomtype
	atom_t * atom_ptr;

  param_t * rval  = calloc(1, sizeof(param_t)); //return value
	memnullcheck(rval,sizeof(param_t), __LINE__-1, __FILE__);
	param_t * param_ptr = rval; //used to seek through param list
	param_t * prev_param_ptr;

  for(atom_ptr = system->molecules->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

     if(strcasecmp(atom_ptr->atomtype, last_atomtype))
     {
          strcpy(param_ptr->atomtype, atom_ptr->atomtype);

          param_ptr->charge  = atom_ptr->charge;
          param_ptr->epsilon = atom_ptr->epsilon;
          param_ptr->sigma   = atom_ptr->sigma;
          param_ptr->omega   = atom_ptr->omega;

          param_ptr->last_charge  = param_ptr->charge;
          param_ptr->last_epsilon = param_ptr->epsilon;
          param_ptr->last_sigma   = param_ptr->sigma;
          param_ptr->last_omega   = param_ptr->omega;

          param_ptr->next = calloc(1, sizeof(param_t));
					memnullcheck(param_ptr->next,sizeof(param_t), __LINE__-1, __FILE__);
          prev_param_ptr = param_ptr;
          param_ptr = param_ptr->next;

          strcpy(last_atomtype, atom_ptr->atomtype);
		}
	}

	prev_param_ptr->next = NULL;
	free(param_ptr);

	//count ntypes
  for(param_ptr = rval; param_ptr; param_ptr = param_ptr->next) {
      for(atom_ptr = system->molecules->atoms; atom_ptr; atom_ptr = atom_ptr->next)
          if(!strcasecmp(param_ptr->atomtype, atom_ptr->atomtype))
              ++(param_ptr->ntypes);
  }

	return rval;
}

// randomly perturb parameters for surface fitting
void surf_perturb ( system_t * system, double quadrupole, qshiftData_t * qshiftData, param_t * params ) {

	param_t * param_ptr; 

  //override default values if specified in input file
  double scale_q = ((system->surf_scale_q_on) ? system->surf_scale_q : 0.0 );
  double scale_r = ((system->surf_scale_r_on) ? system->surf_scale_r : SCALE_R );
  double scale_epsilon = ((system->surf_scale_epsilon_on) ? system->surf_scale_epsilon : SCALE_EPSILON );
  double scale_sigma = ((system->surf_scale_sigma_on) ? system->surf_scale_sigma : SCALE_SIGMA );
  double scale_omega = ((system->surf_scale_omega_on) ? system->surf_scale_omega : SCALE_OMEGA );

	double delta_q;

	int nchargedsites=0;
	double initcharge;

  //if qshift is on, determine the qshift parameters for this time step
  // these will be applied in the following loop
  if ( system->surf_qshift_on ) qshift_do(system,qshiftData,scale_r,quadrupole);

  //how much to adjust q by, if neccessary
  if ( system->surf_scale_q_on ) {
    delta_q = scale_q * (0.5 - get_rand());
		//determine how many charged sites are present
		for (param_ptr = params; param_ptr; param_ptr = param_ptr->next ) {
			if(!strncasecmp(param_ptr->atomtype,"H2E",3) || !strncasecmp(param_ptr->atomtype,"H2Q",3)){
				nchargedsites+=param_ptr->ntypes;
				initcharge=param_ptr->charge;
			}
		}
	}

  // randomly perturb the parameters
	for(param_ptr = params; param_ptr; param_ptr = param_ptr->next) {

		if ( param_ptr->epsilon > 0.0 )
			param_ptr->epsilon += scale_epsilon*(0.5 - get_rand());
		if(param_ptr->epsilon < 0.0) param_ptr->epsilon = param_ptr->last_epsilon;

		//only compare first three characters, so we don't include H2E, H2EP, H2Q, H2QP, etc.
		if((param_ptr->ntypes > 1) && 
				strncasecmp(param_ptr->atomtype, "H2E", 3) 
				&& strncasecmp(param_ptr->atomtype, "H2Q", 3) )
			param_ptr->dr = scale_r*(0.5 - get_rand());
			else
			param_ptr->dr = 0;

		//if polarvdw is on, also adjust omega
		if  ( system->polarvdw ) {
			if( param_ptr->omega > 0.0 )
				param_ptr->omega += scale_omega*(0.5 - get_rand());
			if( param_ptr->omega < 0.0 ) param_ptr->omega = param_ptr->last_omega;
		}
		//otherwise adjust sigma
		else {
			if ( param_ptr->sigma > 0.0 )
				param_ptr->sigma += scale_sigma*(0.5 - get_rand());
			if( param_ptr->sigma < 0.0 ) param_ptr->sigma = param_ptr->last_sigma;
		}

		//if qshift is on, we will force changes on some parameters
		if( system->surf_qshift_on ) {
			if(!strncasecmp(param_ptr->atomtype,"H2G",3)) 
				param_ptr->charge = qshiftData->qH2G;
			if(!strncasecmp(param_ptr->atomtype,"H2Q",3)) {
				param_ptr->charge = qshiftData->qH2Q;
				param_ptr->dr = qshiftData->drH2Q;
			}
		}

		//adjust charges and charge positions without constraining quadrupole
		if ( system->surf_scale_q_on ) {
			if(!strncasecmp(param_ptr->atomtype,"H2G",3)) {
				param_ptr->charge = -nchargedsites * (initcharge + delta_q);
			}
			if(!strncasecmp(param_ptr->atomtype,"H2Q",3)) {
				if ( param_ptr->charge != 0 ) param_ptr->charge += delta_q;
				param_ptr->dr += scale_r * ( 0.5 - get_rand() );
			}
			if(!strncasecmp(param_ptr->atomtype,"H2E",3)) {
				if ( param_ptr->charge != 0 ) param_ptr->charge += delta_q;
			}
		}
	}
	return ;
}

// output curves for best discovered fit
void output_fit ( int nCurves, int nPoints, curveData_t * curve, double max_energy, double * r_input  ) {

	int i, j;

  printf( "#r-value  " );
  for( i=0; i<nCurves; i++ )
  	printf( "#%s  ", curve[i].id );
  printf( "\n" );

  for(i = 0; i < nPoints; i++)
  {
    printf("%f", r_input[i]);
      for( j=0; j<nCurves; j++ )
      {
//        if(curve[j].global[i] < max_energy)  //why would we want to do this?
          printf(" %f", curve[j].global[i]);
  //      else
    //      printf(" %f", max_energy);
      }
      printf("\n");
      fflush(stdout);
  }

	return;
}

//apply trial parameters
void get_curves ( system_t * system, int nCurves, curveData_t * curve, double r_min, double r_max, double r_inc ) {

	int i;

  for( i=0; i<nCurves; i++ ) {

      // Apply rotation associated with i-th curve
      surface_dimer_geometry(system, r_min, curve[i].alpha1, curve[i].beta1, curve[i].gamma1,
                                            curve[i].alpha2, curve[i].beta2, curve[i].gamma2 );
      // Perform the calculation
      surface_curve(system, r_min, r_max, r_inc, curve[i].output);

      // Restore the molecule to its initial orientation
      surface_dimer_geometry(system, r_min, -curve[i].gamma1, -curve[i].beta1, -curve[i].alpha1,
                                            -curve[i].gamma2, -curve[i].beta2, -curve[i].alpha2 );
  }
  // Reset to origin
  surface_dimer_geometry( system, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );

	return;
}

void revert_parameters ( system_t * system, param_t * params ) {

	param_t * param_ptr;

	for(param_ptr = params; param_ptr; param_ptr = param_ptr->next) {

		param_ptr->charge  = param_ptr->last_charge;
		param_ptr->epsilon = param_ptr->last_epsilon;
		param_ptr->sigma   = param_ptr->last_sigma;
		param_ptr->omega   = param_ptr->last_omega;

		// the dr parameter cannot be reassigned, but rather has to be "undone"
		param_ptr->dr      = -param_ptr->dr;
	}

	// Apply
	surface_dimer_parameters( system, params );

	return;
}

//output new geometry and store optimal curves
void new_global_min ( system_t * system, int nCurves, int nPoints, curveData_t * curve ) {
	int j, i;
	// output the new geometry
	write_molecules(system, "fit_geometry.pqr");
	for( j = 0; j<nCurves; j++ )
		for(i = 0; i < nPoints; i++)
			curve[j].global[i] = curve[j].output[i];

	return;
}

// when accepting a move, store parameters as the "last good values"
void apply_new_parameters ( param_t * params ) {
	param_t * param_ptr;

	for(param_ptr = params; param_ptr; param_ptr = param_ptr->next) {
			param_ptr->last_charge  = param_ptr->charge;
			param_ptr->last_epsilon = param_ptr->epsilon;
			param_ptr->last_sigma   = param_ptr->sigma;
			param_ptr->last_omega = param_ptr->omega;
	}

	return;
}

// free all mem used in surface_fit
void free_all_mem 
( int nCurves, curveData_t * curve, param_t * params, qshiftData_t * qshiftData, double * r_input ) {

	int i;
	param_t * param_ptr, * next;

// free all memory used
  for( i=0; i<nCurves; i++ ) {
      free( curve[i].input  );
      free( curve[i].output );
      free( curve[i].global );
			free( curve[i].filename );
			free( curve[i].id );
  }

  // Free the param linked list
  param_ptr = params;
  while( param_ptr ) {
      next = param_ptr->next;
      free(param_ptr);
      param_ptr = next;
  }

  free(qshiftData);
  free(r_input);
  free(curve);

	return;
}

