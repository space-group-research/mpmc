/*

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

void error(char *msg) {

	if(!rank) fprintf(stderr, "%s", msg);
	fflush(stderr);

}

void output(char *msg) {

	if(!rank) printf("%s", msg);
	fflush(stdout);

}

char * make_filename ( char * basename, int fileno ) {
	size_t len = strlen(basename);
	size_t outlen;
	char * rval = NULL;
	char START[MAXLINE], STOP[MAXLINE];
	int i, j;

	if ( ! strncmp("/dev/null",basename,9) ) {
		rval = malloc(10*sizeof(char));
		memnullcheck(rval,10*sizeof(char),__LINE__-1, __FILE__);
		sprintf(rval,"/dev/null");
		return rval;
	}
	else 
	//check if the file has a three character extension
		if ( len - 4 > 0 )
			if ( basename[len-4] == '.' ) { //only check this if previous check passes so we don't have a memory fault
				for ( i=0; i < len-4; i++ ) //format pre-extension
					START[i]=basename[i];
				START[i]='\0';
				for ( i=len-4, j=0; i<len; i++, j++ ) //format extension
					STOP[j]=basename[i];
				STOP[j]='\0';
				//set string length
				outlen=strlen(START)+strlen(STOP)+7;
				rval = malloc(outlen*sizeof(char));
				memnullcheck(rval,outlen*sizeof(char),__LINE__-1, __FILE__);
				//make filename
				sprintf(rval,"%s-%05d%s", START, fileno, STOP);
			}
	//if rval is still NULL, then it's neither /dev/null nor has a proper 3-character file extension
	//just add the number to the end
	if ( rval == NULL ) {
		outlen = len + 7;
		rval = malloc(outlen*sizeof(char));
		memnullcheck(rval,outlen*sizeof(char),__LINE__-1, __FILE__);
		//make filename
		sprintf(rval,"%s-%05d", basename, fileno);
	}
	
	return rval;	
}

int open_traj_file( system_t * system ) {
	char linebuf[MAXLINE];

	if(system->traj_output) {
#ifdef MPI /*each node will write it's own file*/
		char * filename;
		filename = make_filename(system->traj_output,rank);
		system->file_pointers.fp_traj = fopen(filename, "w");
		filecheck(system->file_pointers.fp_traj,filename, WRITE);
		free(filename);
#else
		system->file_pointers.fp_traj = fopen(system->traj_output, "w");
		filecheck(system->file_pointers.fp_traj,system->traj_output,WRITE);
#endif /*MPI*/
	}

	return 0;
}

int open_files(system_t *system) {

	if(system->energy_output) {
		system->file_pointers.fp_energy = fopen(system->energy_output, "w");
		filecheck(system->file_pointers.fp_energy,system->energy_output,WRITE);
		fprintf(system->file_pointers.fp_energy,
			"#step #energy #coulombic #rd #polar #vdw #kinetic #temp #N #spin_ratio #volume\n");
	}

	if(system->energy_output_csv) {
		system->file_pointers.fp_energy_csv = fopen(system->energy_output_csv, "w");
		filecheck(system->file_pointers.fp_energy_csv,system->energy_output_csv,WRITE);
		fprintf(system->file_pointers.fp_energy_csv,
			"#step,#energy,#coulombic,#rd,#polar,#vdw,#kinetic,#temp,#N,#spin_ratio,#volume\n");
	}

	if(system->dipole_output) {
		system->file_pointers.fp_dipole = fopen(system->dipole_output, "w");
		filecheck(system->file_pointers.fp_dipole,system->dipole_output,WRITE);
	}

	if(system->field_output) {
		system->file_pointers.fp_field = fopen(system->field_output, "w");
		filecheck(system->file_pointers.fp_field,system->field_output,WRITE);
	}

	if(system->histogram_output) {
		system->file_pointers.fp_histogram = fopen(system->histogram_output, "w");
		filecheck(system->file_pointers.fp_histogram,system->histogram_output,WRITE);
	}

	if(system->frozen_output) {
		system->file_pointers.fp_frozen = fopen(system->frozen_output, "w");
		filecheck(system->file_pointers.fp_frozen,system->frozen_output,WRITE);
		//go ahead and write the frozen lattice configuration now

		if(system->file_pointers.fp_frozen) write_frozen(system->file_pointers.fp_frozen,system);
		fclose(system->file_pointers.fp_frozen);
	}

	return(0);

}

void close_files(system_t *system) {

	int i;

	if(system->file_pointers.fp_energy) fclose(system->file_pointers.fp_energy);
	if(system->file_pointers.fp_energy_csv) fclose(system->file_pointers.fp_energy_csv);
	if(system->file_pointers.fp_dipole) fclose(system->file_pointers.fp_dipole);
	if(system->file_pointers.fp_field) fclose(system->file_pointers.fp_field);
	if(system->file_pointers.fp_histogram) fclose(system->file_pointers.fp_histogram);
// no need to keep this file open the whole time (afaik) --kmclaugh
//if(system->file_pointers.fp_frozen) fclose(system->file_pointers.fp_frozen);
	if(system->file_pointers.fp_traj) fclose(system->file_pointers.fp_traj);

}

/* enforce molecule wrapping around the periodic boundaries on output - keep the atoms together */
int wrapall(molecule_t *molecules, pbc_t *pbc) {

	int i, j;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	double d[3], dimg[3];


	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {

		if(!molecule_ptr->frozen) {
			/* get the minimum imaging distance for the com */
			for(i = 0; i < 3; i++) {
				for(j = 0, d[i] = 0; j < 3; j++) {
					d[i] += pbc->reciprocal_basis[i][j]*molecule_ptr->com[j];
				}
				d[i] = rint(d[i]);
			}

			for(i = 0; i < 3; i++) {
				for(j = 0, dimg[i] = 0; j < 3; j++)
					dimg[i] += pbc->basis[i][j]*d[j];

				/* store the wrapped com coordinate */
				molecule_ptr->wrapped_com[i] = dimg[i];
			}

			/* apply the distance to all of the atoms of this molecule */
			for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
				for(i = 0; i < 3; i++)
					atom_ptr->wrapped_pos[i] = atom_ptr->pos[i] - dimg[i];
			}

		} else {

			/* don't wrap frozen */
			for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
				for(i = 0; i < 3; i++)
					atom_ptr->wrapped_pos[i] = atom_ptr->pos[i];
			}

		}

	} /* molecule */


	return(0);

}

/* ensure that the SPECTRE charges are all pulled within the restricted domain */
void spectre_wrapall(system_t *system) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	int p;
	double target[3];
	double d[3], l;

	/* boxlength */
	l = 2.0*system->spectre_max_target;

	/* get the coordinates of the target particle to wrap around */
	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

			if(atom_ptr->target)
				for(p = 0; p < 3; p++)
					target[p] = atom_ptr->pos[p];

		}
	}

	/* wrap SPECTRE charges within the box */
	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

			if(atom_ptr->spectre) {

				for(p = 0; p < 3; p++) {
					d[p] = atom_ptr->pos[p] - target[p];
					atom_ptr->pos[p] -= l*rint(d[p]/l);
				}

			}

		} /* for atom */
	} /* for molecule */


}


/* write out the final system state as a PQR file */
int write_molecules(system_t *system, char *filename) {

	int atom_box, molecule_box, p, q;
	double box_pos[3], box_occupancy[3];
	int l, m, n, box_labels[2][2][2], diff;
	char linebuf[MAXLINE];
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	FILE *fp;
	int i, j, k;
	int ext_output = 0; // By default, PDB compliant coordinates are printed (%8.3f), else extended output is used (%11.6f)

	fp = fopen(filename, "w");
	filecheck(fp,filename,WRITE);

	/* Check if extended coordinate output is needed (CRC) */
	if(system->long_output)
		ext_output = 1;
	else if( (system->pbc->basis[0][0] >= 200.0) || (system->pbc->basis[0][1] >= 200.0) || (system->pbc->basis[0][2] >= 200.0) || 
	         (system->pbc->basis[1][0] >= 200.0) || (system->pbc->basis[1][1] >= 200.0) || (system->pbc->basis[1][2] >= 200.0) || 
	         (system->pbc->basis[2][0] >= 200.0) || (system->pbc->basis[2][1] >= 200.0) || (system->pbc->basis[2][2] >= 200.0) )
		ext_output = 1;

	/* write pqr */
	for(molecule_ptr = system->molecules, i = 1, j = 1; molecule_ptr; molecule_ptr = molecule_ptr->next, j++) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next, i++) {

			fprintf(fp, "ATOM  ");
			fprintf(fp, "%5d", i);		/* give each one a unique id */
			fprintf(fp, " %-4.45s", atom_ptr->atomtype);
			fprintf(fp, " %-3.3s ", molecule_ptr->moleculetype);
			if(atom_ptr->adiabatic)
				fprintf(fp, "%-1.1s", "A");
			else if(atom_ptr->frozen)
				fprintf(fp, "%-1.1s", "F");
			else if(atom_ptr->spectre)
				fprintf(fp, "%-1.1s", "S");
			else if(atom_ptr->target)
				fprintf(fp, "%-1.1s", "T");
			else
				fprintf(fp, "%-1.1s", "M");
			if(system->independent_particle)
				fprintf(fp, " %4d   ", i);
			else
				fprintf(fp, " %4d   ", j);		/* give each molecule a unique id */

			/* Regular (PDB compliant) Coordinate Output */
			if( (system->wrapall) && (ext_output == 0) ) {
				fprintf(fp, "%8.3f", atom_ptr->wrapped_pos[0]);
				fprintf(fp, "%8.3f", atom_ptr->wrapped_pos[1]);
				fprintf(fp, "%8.3f", atom_ptr->wrapped_pos[2]);
			} else if(ext_output == 0){
				fprintf(fp, "%8.3f", atom_ptr->pos[0]);
				fprintf(fp, "%8.3f", atom_ptr->pos[1]);
				fprintf(fp, "%8.3f", atom_ptr->pos[2]);
			}

			/* Extended (PQR) Coordinate Output */
			if( (system->wrapall) && (ext_output == 1) ) {
				fprintf(fp, "%11.6f ", atom_ptr->wrapped_pos[0]);
				fprintf(fp, "%11.6f ", atom_ptr->wrapped_pos[1]);
				fprintf(fp, "%11.6f ", atom_ptr->wrapped_pos[2]);
			} else if (ext_output == 1) {
				fprintf(fp, "%11.6f ", atom_ptr->pos[0]);
				fprintf(fp, "%11.6f ", atom_ptr->pos[1]);
				fprintf(fp, "%11.6f ", atom_ptr->pos[2]);
			}
			fprintf(fp, " %8.5f", atom_ptr->mass);
			fprintf(fp, " %8.5f", atom_ptr->charge/E2REDUCED);	/* convert charge back to real units */
			fprintf(fp, " %8.5f", atom_ptr->polarizability);
			fprintf(fp, " %8.5f", atom_ptr->epsilon);
			fprintf(fp, " %8.5f", atom_ptr->sigma);
			fprintf(fp, " %8.5f", atom_ptr->omega);
			fprintf(fp, " %8.5f", atom_ptr->gwp_alpha);
			fprintf(fp, "\n");

		}
	}

	if(system->wrapall) {

		/* output the box coords as virtual particles for visualization */
		atom_box = i;
		molecule_box = j;
		for(i = 0; i < 2; i++) {
			for(j = 0; j < 2; j++) {
				for(k = 0; k < 2; k++) {

				/* make this frozen */
				fprintf(fp, "ATOM  ");
				fprintf(fp, "%5d", atom_box);
				fprintf(fp, " %-4.45s", "X");
				fprintf(fp, " %-3.3s ", "BOX");
				fprintf(fp, "%-1.1s", "F");
				fprintf(fp, " %4d   ", molecule_box);

				/* box coords */
				box_occupancy[0] = ((double)i) - 0.5;
				box_occupancy[1] = ((double)j) - 0.5;
				box_occupancy[2] = ((double)k) - 0.5;

				for(p = 0; p < 3; p++)
					for(q = 0, box_pos[p] = 0; q < 3; q++)
						box_pos[p] += system->pbc->basis[p][q]*box_occupancy[q];

				for(p = 0; p < 3; p++)
					if(ext_output == 0)
						fprintf(fp, "%8.3f", box_pos[p]);
					else
						fprintf(fp, "%11.6f ", box_pos[p]);

				/* null interactions */
				fprintf(fp, " %8.4f", 0.0);
				fprintf(fp, " %8.4f", 0.0);
				fprintf(fp, " %8.5f", 0.0);
				fprintf(fp, " %8.5f", 0.0);
				fprintf(fp, " %8.5f", 0.0);
				fprintf(fp, "\n");

				box_labels[i][j][k] = atom_box;
				++atom_box;

				}
			}
		}

		for(i = 0; i < 2; i++) {
			for(j = 0; j < 2; j++) {
				for(k = 0; k < 2; k++) {
	
					for(l = 0; l < 2; l++) {
						for(m = 0; m < 2; m++) {
							for(n = 0; n < 2; n++) {

									diff = fabs(i - l) + fabs(j - m) + fabs(k - n);
									if(diff == 1)
										fprintf(fp, "CONECT %4d %4d\n", box_labels[i][j][k], box_labels[l][m][n]);
	
							} /* n */
						} /* m */
					} /* l */


				} /* k */
			} /* j */
		} /* i */

	} /* if wrapall */

	/*write basis to the output file. needed for restarting NPT jobs, or whatever.*/
	fprintf(fp, "REMARK BOX BASIS[0] = %20.14lf %20.14lf %20.14lf\n", 
		system->pbc->basis[0][0], system->pbc->basis[0][1], system->pbc->basis[0][2]);
	fprintf(fp, "REMARK BOX BASIS[1] = %20.14lf %20.14lf %20.14lf\n",
		system->pbc->basis[1][0], system->pbc->basis[1][1], system->pbc->basis[1][2]);
	fprintf(fp, "REMARK BOX BASIS[2] = %20.14lf %20.14lf %20.14lf\n", 
		system->pbc->basis[2][0], system->pbc->basis[2][1], system->pbc->basis[2][2]);

		/* output the connectivity information */
	fprintf(fp, "END\n");
	fflush(fp);

	fclose(fp);
	return(0);

}

void write_states(FILE *fp, system_t * system) {

	molecule_t *molecules = system->molecules;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	char linebuf[MAXLINE];
	double box_pos[3], box_occupancy[3];
	int l, m, n, box_labels[2][2][2], diff;
	int i, j, k;
	int atom_box, molecule_box, p, q;
	int num_frozen_molecules, num_moveable_molecules;
	int num_frozen_atoms, num_moveable_atoms;
	int ext_output = 0; // By default, PDB compliant coordinates are printed (%8.3f), else extended output is used (%11.6f)

	//don't bother if we'd be writing to /dev/null
	if ( ! strncmp("/dev/null",system->traj_output,9) ) return;

	/* count the number of molecules, atoms, etc. */
	num_frozen_molecules = 0, num_moveable_molecules = 0;
	num_frozen_atoms = 0, num_moveable_atoms = 0;
	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {

		if(molecule_ptr->frozen) {
			++num_frozen_molecules;
			for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
				++num_frozen_atoms;
		} else {
			++num_moveable_molecules;
			for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
				++num_moveable_atoms;
		}

	}

	fprintf(fp, "REMARK step=%d\n", system->step);
#ifdef MPI
	fprintf(fp, "REMARK node=%d\n", rank);
#endif

	fprintf(fp, "REMARK total_molecules=%d, total_atoms=%d\n", 
		(num_frozen_molecules + num_moveable_molecules), (num_frozen_atoms + num_moveable_atoms));
	fprintf(fp, "REMARK frozen_molecules=%d, moveable_molecules=%d\n", 
		num_frozen_molecules, num_moveable_molecules);
	fprintf(fp, "REMARK frozen_atoms=%d, moveable_atoms=%d\n", 
		num_frozen_atoms, num_moveable_atoms);

	/* Check if extended coordinate output is needed (CRC) */
	if(system->long_output)
		ext_output = 1;
	else if( (system->pbc->basis[0][0] >= 200.0) || (system->pbc->basis[0][1] >= 200.0) || (system->pbc->basis[0][2] >= 200.0) || 
	       (system->pbc->basis[1][0] >= 200.0) || (system->pbc->basis[1][1] >= 200.0) || (system->pbc->basis[1][2] >= 200.0) || 
	       (system->pbc->basis[2][0] >= 200.0) || (system->pbc->basis[2][1] >= 200.0) || (system->pbc->basis[2][2] >= 200.0) )
		ext_output = 1;

	/* write pqr formatted states */
	for(molecule_ptr = molecules, i = 1, j = 1; molecule_ptr; molecule_ptr = molecule_ptr->next, j++) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next, i++) {

			fprintf(fp, "ATOM  ");
			fprintf(fp, "%5d", i);		/* give each one a unique id */
			fprintf(fp, " %-4.45s", atom_ptr->atomtype);
			fprintf(fp, " %-3.3s ", molecule_ptr->moleculetype);
			if(atom_ptr->adiabatic)
				fprintf(fp, "%-1.1s", "A");
			else if(atom_ptr->frozen)
				fprintf(fp, "%-1.1s", "F");
			else if(atom_ptr->spectre)
				fprintf(fp, "%-1.1s", "S");
			else if(atom_ptr->target)
				fprintf(fp, "%-1.1s", "T");
			else
				fprintf(fp, "%-1.1s", "M");
			fprintf(fp, "%4d    ", j);		/* give each molecule a unique id */

			if(ext_output == 0) {
				/* Regular (PDB compliant) Coordinate Output */
				fprintf(fp, "%8.3f", atom_ptr->wrapped_pos[0]);
				fprintf(fp, "%8.3f", atom_ptr->wrapped_pos[1]);
				fprintf(fp, "%8.3f", atom_ptr->wrapped_pos[2]);
			} else {
				/* Extended (PQR) Coordinate Output */
				fprintf(fp, "%11.6f ", atom_ptr->wrapped_pos[0]);
				fprintf(fp, "%11.6f ", atom_ptr->wrapped_pos[1]);
				fprintf(fp, "%11.6f ", atom_ptr->wrapped_pos[2]);
			}

			fprintf(fp, " %8.4f", atom_ptr->mass);
			fprintf(fp, " %8.4f", atom_ptr->charge/E2REDUCED);	/* convert charge back to real units */
			fprintf(fp, " %8.5f", atom_ptr->polarizability);
			fprintf(fp, " %8.5f", atom_ptr->epsilon);
			fprintf(fp, " %8.5f", atom_ptr->sigma);
			fprintf(fp, " %8.5f", atom_ptr->omega);
			fprintf(fp, " %8.5f", atom_ptr->gwp_alpha);
			fprintf(fp, "\n");

		}
	}

	/*write basis to the output file. needed for restarting NPT jobs, or whatever.*/
	fprintf(fp, "REMARK BOX BASIS[0] = %20.14lf %20.14lf %20.14lf\n", 
		system->pbc->basis[0][0], system->pbc->basis[0][1], system->pbc->basis[0][2]);
	fprintf(fp, "REMARK BOX BASIS[1] = %20.14lf %20.14lf %20.14lf\n",
		system->pbc->basis[1][0], system->pbc->basis[1][1], system->pbc->basis[1][2]);
	fprintf(fp, "REMARK BOX BASIS[2] = %20.14lf %20.14lf %20.14lf\n", 
		system->pbc->basis[2][0], system->pbc->basis[2][1], system->pbc->basis[2][2]);

	fprintf(fp, "ENDMDL\n");
	fflush(fp);

}

int write_performance(int i, system_t *system) {

	char linebuf[MAXLINE];
	time_t current_time;
	double sec_step;
	static time_t last_time;
	static int last_step;

	current_time = time(NULL);
	if(i > system->corrtime) {

		sec_step = difftime(current_time, last_time)/((double)(i - last_step));

		if(system->ensemble == ENSEMBLE_UVT) {
			sprintf(linebuf, "OUTPUT: Grand Canonical Monte Carlo simulation running on %d cores\n", size);
			output(linebuf);
		} else {
			sprintf(linebuf, "OUTPUT: Canonical Monte Carlo simulation running on %d cores\n", size);
			output(linebuf);
		}
		float nsf = system->numsteps; /* numsteps as a float (nsf) */
		sprintf(linebuf, "OUTPUT: Root collecting statistics at %s", ctime(&current_time));
		output(linebuf);
		sprintf(linebuf, "OUTPUT: Completed step %d/%d  (%.3f %%)\n", i, system->numsteps, (i/nsf)*100);
		output(linebuf);
		sprintf(linebuf, "OUTPUT: %f sec/step, ETA = %.3f hrs\n", sec_step, sec_step*(system->numsteps - i)/3600.0);
		output(linebuf);
		if ( system->simulated_annealing) {
			sprintf(linebuf, "OUTPUT: Current temperature = %.5f\n", system->temperature);
			output(linebuf);
		}

	}

	last_step = i;
	last_time = current_time;

	return(0);

}

void write_observables(FILE *fp_energy, system_t * system, observables_t * observables) {

	fprintf(fp_energy, "%d %f %f %f %f %f %f %f %f %f %f", 
		system->step,
		observables->energy, 
		observables->coulombic_energy, 
		observables->rd_energy, 
		observables->polarization_energy, 
		observables->vdw_energy, 
		observables->kinetic_energy, 
		observables->temperature, 
		observables->N, 
		observables->spin_ratio, 
		observables->volume);
	fprintf(fp_energy, "\n");
	fflush(fp_energy);
}

void write_observables_csv(FILE *fp_energy_csv, system_t * system, observables_t * observables) {

	fprintf(fp_energy_csv, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", 
		system->step,
		observables->energy, 
		observables->coulombic_energy, 
		observables->rd_energy, 
		observables->polarization_energy, 
		observables->vdw_energy, 
		observables->kinetic_energy, 
		observables->temperature, 
		observables->N, 
		observables->spin_ratio, 
		observables->volume);
	fprintf(fp_energy_csv, "\n");
	fflush(fp_energy_csv);
}

/* output each molecular dipole (in debye) per line */
void write_dipole(FILE *fp_dipole, molecule_t *molecules) {

	int p;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	double dipole[3];

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(p = 0; p < 3; p++) dipole[p] = 0;
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(p = 0; p < 3; p++)
				dipole[p] += atom_ptr->mu[p];
		}
		if(!molecule_ptr->frozen) fprintf(fp_dipole, "%f %f %f\n", dipole[0]/DEBYE2SKA, dipole[1]/DEBYE2SKA, dipole[2]/DEBYE2SKA);
	}
	fflush(fp_dipole);
	return;
}

/* output the total molecular electrostatic field (in e/A) per line) */
void write_field(FILE *fp_field, molecule_t *molecules) {

	int p;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	double field[3];

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(p = 0; p < 3; p++) field[p] = 0;
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(p = 0; p < 3; p++)
				field[p] += atom_ptr->ef_static[p] + atom_ptr->ef_induced[p];
		}
		if(!molecule_ptr->frozen) fprintf(fp_field, "%f %f %f\n", field[0]/E2REDUCED, field[1]/E2REDUCED, field[2]/E2REDUCED);
	}

	return;
}



int write_averages(system_t *system) {

	avg_observables_t *averages;

	averages = system->avg_observables;

	if(averages->boltzmann_factor > 0.0)
		printf("OUTPUT: BF = %.3f +- %.3f\n", averages->boltzmann_factor, 0.5*averages->boltzmann_factor_error);

	if(averages->acceptance_rate > 0.0) {
		printf("OUTPUT: AR = %.3f (%.3f I/ %.3f R/ %.3f D", 
			averages->acceptance_rate, averages->acceptance_rate_insert, 
			averages->acceptance_rate_remove, averages->acceptance_rate_displace);
		if(averages->acceptance_rate_adiabatic > 0.0) printf("/ %.3f A", averages->acceptance_rate_adiabatic);
		if(averages->acceptance_rate_spinflip > 0.0) printf("/ %.3f S", averages->acceptance_rate_spinflip);
		if(averages->acceptance_rate_volume > 0.0) printf("/ %.3f V", averages->acceptance_rate_volume);
		printf(")\n");
	}

	if(averages->cavity_bias_probability > 0.0)
		printf("OUTPUT: Cavity bias probability = %.5f +- %.5f\n", 
			averages->cavity_bias_probability, 0.5*averages->cavity_bias_probability_error);

	if(system->gwp) 
		printf("OUTPUT: total energy = %.3f +- %.3f eV\n", averages->energy/EV2K, 0.5*averages->energy_error/EV2K);
	else
		printf("OUTPUT: potential energy = %.3f +- %.3f K\n", averages->energy, 0.5*averages->energy_error);

	if(averages->coulombic_energy != 0.0) {
		if(system->gwp) 
			printf("OUTPUT: electrostatic energy = %.3f +- %.3f eV\n",
				averages->coulombic_energy/EV2K, 0.5*averages->coulombic_energy_error/EV2K);
		else
			printf("OUTPUT: electrostatic energy = %.3f +- %.3f K\n",
				averages->coulombic_energy, 0.5*averages->coulombic_energy_error);
	}

	if(averages->rd_energy != 0.0)
		printf("OUTPUT: repulsion/dispersion energy = %.3f +- %.3f K\n", 
			averages->rd_energy, 0.5*averages->rd_energy_error);

	if(averages->polarization_energy != 0.0) {
		printf("OUTPUT: polarization energy = %.5f +- %.5f K", 
			averages->polarization_energy, 0.5*averages->polarization_energy_error);
		if(averages->dipole_rrms_error != 0.0 ) {
			printf(" (iterations = %.1f +- %.1f rrms = %e +- %e)", 
				averages->polarization_iterations, 0.5*averages->polarization_iterations_error, 
				averages->dipole_rrms, 0.5*averages->dipole_rrms_error);
		}
		else if(averages->polarization_iterations != 0.0) {
			printf(" (iterations = %.1f +- %.1f)", 
				averages->polarization_iterations, 0.5*averages->polarization_iterations_error);
		}

		printf("\n");
	}

#ifdef VDW
		printf("OUTPUT: (coupled-dipole) vdw energy = %.5f +- %.5f K\n", 
			averages->vdw_energy, 0.5*averages->vdw_energy_error);
#endif

	if(averages->kinetic_energy > 0.0) {
		if(system->gwp)
			printf("OUTPUT: kinetic energy = %.3f +- %.3f eV\n", 
				averages->kinetic_energy/EV2K, 0.5*averages->kinetic_energy_error/EV2K);
		else
			printf("OUTPUT: kinetic energy = %.3f +- %.3f K\n", 
				averages->kinetic_energy, 0.5*averages->kinetic_energy_error);

		printf("OUTPUT: temperature = %.3f +- %.3f K\n", 
			averages->temperature, 0.5*averages->temperature);
	}

	printf("OUTPUT: N = %.3f +- %.3f molecules\n", averages->N, 0.5*averages->N_error);

	if ( system->sorbateCount == 1 ) { //all based on calculations with assume only one type of sorbate
		printf("OUTPUT: density = %.5f +- %.5f g/cm^3\n", averages->density, 0.5*averages->density_error);
		if(averages->pore_density != 0.0 && system->ensemble != ENSEMBLE_NPT)
			printf("OUTPUT: pore density = %.5f +- %.5f g/cm^3\n", averages->pore_density, 0.5*averages->pore_density_error);
		if(averages->percent_wt > 0.0 ) {
			printf("OUTPUT: wt %% = %.5f +- %.5f %%\n", averages->percent_wt, 0.5*averages->percent_wt_error);
			printf("OUTPUT: wt %% (ME) = %.5f +- %.5f %%\n", averages->percent_wt_me, 0.5*averages->percent_wt_me_error);
		}
		if(averages->excess_ratio > 0.0)
			printf("OUTPUT: excess adsorption ratio = %.5f +- %.5f mg/g\n", averages->excess_ratio, 0.5*averages->excess_ratio_error);
		if((averages->qst > 0.0) && isfinite(averages->qst))
			printf("OUTPUT: qst = %.3f kJ/mol\n", averages->qst);
		if((averages->heat_capacity > 0.0) && (isfinite(averages->heat_capacity)))
			printf("OUTPUT: heat capacity = %.5f +- %.5f kJ/mol K\n", averages->heat_capacity, averages->heat_capacity_error);
		if((averages->compressibility > 0.0) && isfinite(averages->compressibility)) {
			printf("OUTPUT: compressibility = %.6f +- %.6f atm^-1\n", averages->compressibility, averages->compressibility_error);
			printf("OUTPUT: bulk modulus = %.6f +- %.6f GPa\n", ATM2PASCALS*1.0e-9/averages->compressibility, 
				ATM2PASCALS*1.0e-9/averages->compressibility_error);
		}
	}

	if(system->ensemble == ENSEMBLE_NPT)
		printf("OUTPUT: volume = %.5f +- %.5f A^3\n", averages->volume, 0.5*averages->volume_error);

	if(averages->spin_ratio > 0.0) {
		printf("OUTPUT: ortho spin ratio = %.3f +- %.3f %%\n", averages->spin_ratio*100.0, averages->spin_ratio_error*100.0);
	}
	
	if( system->sorbateCount > 1 ){
		sorbateAverages_t *sorbate_ptr;
		for( sorbate_ptr = system->sorbateStats.next; sorbate_ptr; sorbate_ptr = sorbate_ptr->next ) {
			printf( "OUTPUT: Stats for %s\n", sorbate_ptr->id  );
			printf( "             Average_N(%s)= %lf\n",          sorbate_ptr->id, sorbate_ptr->avgN          );
			printf( "             Sorbed_Mass(%s)= %lf g/mol\n",  sorbate_ptr->id, sorbate_ptr->sorbed_mass   );
			printf( "             density(%s)= %E g/cm^3\n",      sorbate_ptr->id, sorbate_ptr->density       );
			printf( "             pore_density(%s)= %E g/cm^3\n", sorbate_ptr->id, sorbate_ptr->pore_density  );
			printf( "             excess_ratio(%s)= %E g/cm^3\n", sorbate_ptr->id, sorbate_ptr->excess_ratio  );
			printf( "             wt_%%(%s)= %lf %%\n",           sorbate_ptr->id, sorbate_ptr->percent_wt    );
			printf( "             wt_%%_(ME)(%s)= %lf %%\n",      sorbate_ptr->id, sorbate_ptr->percent_wt_me );
			printf( "             Selectivity(%s)= %lf\n",        sorbate_ptr->id, sorbate_ptr->selectivity   );
		}
	}
	
	printf("\n");
	fflush(stdout);

	return(0);

}

