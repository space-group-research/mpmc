/*

Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

void error(char *msg) {

	if(!rank) fprintf(stderr, "(ERROR) %s", msg);
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

int open_files(system_t *system) {

	if(system->energy_output) {
		system->file_pointers.fp_energy = fopen(system->energy_output, "w");
		filecheck(system->file_pointers.fp_energy,system->energy_output,WRITE);
		fprintf(system->file_pointers.fp_energy,
			"#step #energy #coulombic #rd #polar #vdw #kinetic #kin_temp #N #spin_ratio #volume #core_temp\n");
	}

	if(system->energy_output_csv) {
		system->file_pointers.fp_energy_csv = fopen(system->energy_output_csv, "w");
		filecheck(system->file_pointers.fp_energy_csv,system->energy_output_csv,WRITE);
		fprintf(system->file_pointers.fp_energy_csv,
			"#step,#energy,#coulombic,#rd,#polar,#vdw,#kinetic,#kin_temp,#N,#spin_ratio,#volume,#core_temp\n");
	}

	// if we're just calculating energy or replaying a trajectory, we need no other output files
	if ( system->ensemble == ENSEMBLE_REPLAY || system->ensemble == ENSEMBLE_TE ) return 0;

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

	if(system->file_pointers.fp_energy) fclose(system->file_pointers.fp_energy);
	if(system->file_pointers.fp_energy_csv) fclose(system->file_pointers.fp_energy_csv);
	if(system->file_pointers.fp_histogram) fclose(system->file_pointers.fp_histogram);
	if(system->file_pointers.fp_traj_replay) fclose(system->file_pointers.fp_traj_replay);
	if(system->file_pointers.fp_surf) fclose(system->file_pointers.fp_surf);
/* open on the fly
	if(system->file_pointers.fp_field) fclose(system->file_pointers.fp_field);
	if(system->file_pointers.fp_dipole) fclose(system->file_pointers.fp_dipole);
	if(system->file_pointers.fp_frozen) fclose(system->file_pointers.fp_frozen);
	if(system->file_pointers.fp_traj) fclose(system->file_pointers.fp_traj);
*/

}

int open_surf_traj_file( system_t * system ) { 

        /* Surface trajectory output */
        if(system->surf_output) {
                system->file_pointers.fp_surf = fopen(system->surf_output, "w");
                filecheck(system->file_pointers.fp_surf,system->surf_output,WRITE);
        }   

        return(0);
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
					d[i] += pbc->reciprocal_basis[j][i]*molecule_ptr->com[j];
				}
				d[i] = rint(d[i]);
			}

			for(i = 0; i < 3; i++) {
				for(j = 0, dimg[i] = 0; j < 3; j++)
					dimg[i] += pbc->basis[j][i]*d[j];

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
	double target[3] = {0,0,0};
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

int write_molecules_wrapper(system_t * system, char * filename) {
	int rval;
	char  filenameold[MAXLINE]; 
	FILE * fp;

#ifdef MPI
	int j;
	char * filenameno;
	//make a new filename with the core/or bath number appended
	if ( system->parallel_tempering )
		filenameno=make_filename(filename,system->ptemp->index[rank]); //append bath index to filename
	else
		filenameno=make_filename(filename,rank); //append core index to filename

	//move most recent state file to file.last
	sprintf(filenameold,"%s.last",filenameno);
	rename(filenameno, filenameold);

	//open the file and free the filename string
	fp = fopen(filenameno, "w");
	filecheck(fp,filenameno,WRITE);
	free(filenameno);

	// we write files one at a time to avoid disk congestion
	for ( j=0; j<size; j++ ) {
		MPI_Barrier(MPI_COMM_WORLD);
		if ( j == rank )
			rval = write_molecules(system,fp);
	}

	//free the file pointer
	fclose(fp);

#else //non-MPI

	//move most recent state file to file.last
	sprintf(filenameold,"%s.last",filename);
	rename(filename, filenameold);

	//open the new file
	fp = fopen(filename, "w");
	filecheck(fp,filename,WRITE);

	//write the file
	rval = write_molecules(system,fp);

	//free the file pointer
	fclose(fp);
#endif

	return rval;
}

/* write out the final system state as a PQR file */
int write_molecules(system_t *system, FILE * fp) {

	int atom_box, molecule_box, p, q;
	double box_pos[3], box_occupancy[3];
	int l, m, n, box_labels[2][2][2], diff;
	// char linebuf[MAXLINE];  (unused variable)
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	int i, j, k;
	int ext_output;
	pbc_t * pbc = system->pbc;

	/* Check if extended coordinate output is needed (CRC) */
	// By default, PDB compliant coordinates are printed (%8.3f), else extended output is used (%11.6f)
	if(system->long_output)
		ext_output = 1;
	else if( (system->pbc->basis[0][0] >= 200.0) || (system->pbc->basis[0][1] >= 200.0) || (system->pbc->basis[0][2] >= 200.0) || 
	         (system->pbc->basis[1][0] >= 200.0) || (system->pbc->basis[1][1] >= 200.0) || (system->pbc->basis[1][2] >= 200.0) || 
	         (system->pbc->basis[2][0] >= 200.0) || (system->pbc->basis[2][1] >= 200.0) || (system->pbc->basis[2][2] >= 200.0) )
		ext_output = 1;
	else
		ext_output = 0;

	/* write PBC data */ //VMD uses a weird convention which essentially reverses alpha <-> beta
	fprintf(fp,"CRYST1");
	fprintf(fp,"%9.3f",sqrt(dddotprod(pbc->basis[0], pbc->basis[0])));
	fprintf(fp,"%9.3f",sqrt(dddotprod(pbc->basis[1], pbc->basis[1])));
	fprintf(fp,"%9.3f",sqrt(dddotprod(pbc->basis[2], pbc->basis[2])));
	fprintf(fp,"%7.2f", 180.0/M_PI*acos( dddotprod(pbc->basis[2],pbc->basis[0]) / sqrt( dddotprod(pbc->basis[0], pbc->basis[0]) * dddotprod(pbc->basis[2], pbc->basis[2]) ) ));
	fprintf(fp,"%7.2f", 180.0/M_PI*acos( dddotprod(pbc->basis[1],pbc->basis[2]) / sqrt( dddotprod(pbc->basis[1], pbc->basis[1]) * dddotprod(pbc->basis[2], pbc->basis[2]) ) ));
	fprintf(fp,"%7.2f", 180.0/M_PI*acos( dddotprod(pbc->basis[0],pbc->basis[1]) / sqrt( dddotprod(pbc->basis[1], pbc->basis[1]) * dddotprod(pbc->basis[0], pbc->basis[0]) ) ));
	fprintf(fp,"\n");
	

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
			fprintf(fp, " %8.5f", atom_ptr->c6);
			fprintf(fp, " %8.5f", atom_ptr->c8);
			fprintf(fp, " %8.5f", atom_ptr->c10);
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
						box_pos[p] += system->pbc->basis[q][p]*box_occupancy[q];

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

	return(0);

}

FILE * open_traj_file( system_t * system ) {
	FILE * fp;
	char * filename;
	static int clobber = 1; //if clobber is set, we will overwrite old files

		//open files for append
	if(system->traj_output) {
#ifdef MPI /*each node will write it's own file*/
		if ( system->parallel_tempering )
			filename = make_filename(system->traj_output,system->ptemp->index[rank]); //append bath index to filename
		else 
			filename = make_filename(system->traj_output,rank); //append core index to filename
#else
		filename = system->traj_output;
#endif /*MPI*/

		if ( clobber == 1 ) {
			fp = fopen(filename, "w");
			filecheck(fp,system->traj_output,WRITE);
			clobber = 0; //don't clobber again
		}
		else {
			fp = fopen(filename, "a");
			filecheck(fp,system->traj_output,APPEND);
		}

#ifdef MPI
		free(filename);
#endif
		return fp;
	}

	return NULL;
}

FILE * open_field_file( system_t * system ) {
	FILE * fp;
	char * filename;
	static int clobber = 1; //if clobber is set, we will overwrite old files

		//open files for append
	if(system->field_output) {
#ifdef MPI /*each node will write it's own file*/
		if ( system->parallel_tempering )
			filename = make_filename(system->field_output,system->ptemp->index[rank]); //append bath index to filename
		else 
			filename = make_filename(system->field_output,rank); //append core index to filename
#else
		filename = system->field_output;
#endif /*MPI*/

		if ( clobber == 1 ) {
			fp = fopen(filename, "w");
			filecheck(fp,system->field_output,WRITE);
			clobber = 0; //don't clobber again
		}
		else {
			fp = fopen(filename, "a");
			filecheck(fp,system->field_output,APPEND);
		}

#ifdef MPI
		free(filename);
#endif
		return fp;
	}

	return NULL;
}

FILE * open_dipole_file( system_t * system ) {
	FILE * fp;
	char * filename;
	static int clobber = 1; //if clobber is set, we will overwrite old files

		//open files for append
	if(system->dipole_output) {
#ifdef MPI /*each  will write it's own file*/
		if ( system->parallel_tempering )
			filename = make_filename(system->dipole_output,system->ptemp->index[rank]); //append bath index to filename
		else 
			filename = make_filename(system->dipole_output,rank); //append core index to filename
#else
		filename = system->dipole_output;
#endif /*MPI*/

		if ( clobber == 1 ) {
			fp = fopen(filename, "w");
			filecheck(fp,system->dipole_output,WRITE);
			clobber = 0; //don't clobber again
		}
		else {
			fp = fopen(filename, "a");
			filecheck(fp,system->dipole_output,APPEND);
		}

#ifdef MPI
		free(filename);
#endif
		return fp;
	}

	return NULL;
}

void write_states(system_t * system) {

	molecule_t *molecules = system->molecules;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	// char linebuf[MAXLINE]; (unused variable)
	// double box_pos[3], box_occupancy[3];  (unused variables)
	// int l, m, n, box_labels[2][2][2], diff;  (unused variables) 
	int i, j; // , k;  (unused variable)
	// int atom_box, molecule_box, p, q;  (unused variables)
	int num_frozen_molecules, num_moveable_molecules;
	int num_frozen_atoms, num_moveable_atoms;
	int ext_output = 0; // By default, PDB compliant coordinates are printed (%8.3f), else extended output is used (%11.6f)
	FILE * fp;
	pbc_t * pbc = system->pbc;

	//don't bother if we'd be writing to /dev/null
	if ( ! strncmp("/dev/null",system->traj_output,9) ) return;
 	else fp = open_traj_file(system);

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
	if ( system->parallel_tempering ) 
		fprintf(fp, "REMARK temperature=%.6lf\n",	system->ptemp->templist[system->ptemp->index[rank]]);
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

	/* write PBC data */
	fprintf(fp,"CRYST1");
	fprintf(fp,"%9.3f",sqrt(dddotprod(pbc->basis[0], pbc->basis[0])));
	fprintf(fp,"%9.3f",sqrt(dddotprod(pbc->basis[1], pbc->basis[1])));
	fprintf(fp,"%9.3f",sqrt(dddotprod(pbc->basis[2], pbc->basis[2])));
	fprintf(fp,"%7.2f", 180.0/M_PI*acos( dddotprod(pbc->basis[1],pbc->basis[2]) / sqrt( dddotprod(pbc->basis[1], pbc->basis[1]) * dddotprod(pbc->basis[2], pbc->basis[2]) ) ));
	fprintf(fp,"%7.2f", 180.0/M_PI*acos( dddotprod(pbc->basis[2],pbc->basis[0]) / sqrt( dddotprod(pbc->basis[0], pbc->basis[0]) * dddotprod(pbc->basis[2], pbc->basis[2]) ) ));
	fprintf(fp,"%7.2f", 180.0/M_PI*acos( dddotprod(pbc->basis[0],pbc->basis[1]) / sqrt( dddotprod(pbc->basis[1], pbc->basis[1]) * dddotprod(pbc->basis[0], pbc->basis[0]) ) ));
	fprintf(fp,"\n");

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
			fprintf(fp, " %8.5f", atom_ptr->c6);
			fprintf(fp, " %8.5f", atom_ptr->c8);
			fprintf(fp, " %8.5f", atom_ptr->c10);
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

	fclose(fp);
}

void write_surface_traj(FILE *fpsurf, system_t * system) {

	molecule_t *molecules = system->molecules;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	// char linebuf[MAXLINE];  (unused variable)
	// double box_pos[3], box_occupancy[3]; (unused variables)
	// int l, m, n, box_labels[2][2][2], diff;  (unused variables)
	int i, j; //, k; (unused variable)
	// int atom_box, molecule_box, p, q; (unused variables) 
	// int num_frozen_molecules, num_moveable_molecules;  (unused variables)
	// int num_frozen_atoms, num_moveable_atoms;  (unused variable)
	// int ext_output = 0; // Don't want to use extended output for surface trajectory file.  (unused variable)

	//don't bother if we'd be writing to /dev/null
	if ( ! strncmp("/dev/null",system->surf_output,9) ) 
		return;

	fprintf(fpsurf, "#;");
#ifdef MPI
	fprintf(fpsurf, "Rn%d;", rank);
#endif

	/* write pqr formatted states */
	for(molecule_ptr = molecules, i = 1, j = 1; molecule_ptr; molecule_ptr = molecule_ptr->next, j++) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next, i++) {
		
			/* REDUCED & CONDENSED OUTPUT */
			fprintf(fpsurf, "@,");
			fprintf(fpsurf, "%d,", i);              /* give each one a unique id */
			fprintf(fpsurf, "%s,", atom_ptr->atomtype);
			fprintf(fpsurf, "%s,", molecule_ptr->moleculetype);
			if(atom_ptr->adiabatic)
				fprintf(fpsurf, "%s,", "A");
			else if(atom_ptr->frozen)
				fprintf(fpsurf, "%s,", "F");
			else if(atom_ptr->spectre)
				fprintf(fpsurf, "%s,", "S");
			else if(atom_ptr->target)
				fprintf(fpsurf, "%s,", "T");
			else
				fprintf(fpsurf, "%s,", "M");
			fprintf(fpsurf, "%d,", j);              /* give each molecule a unique id */
			
			/* Regular (PDB compliant) Coordinate Output */
			if( (atom_ptr->wrapped_pos[0] < 0.0005) && (atom_ptr->wrapped_pos[0] > -0.0005) )
				fprintf(fpsurf, "*,");
			else
				fprintf(fpsurf, "%.3f,", atom_ptr->wrapped_pos[0]);
			
			if( (atom_ptr->wrapped_pos[1] < 0.0005) && (atom_ptr->wrapped_pos[1] > -0.0005) )
				fprintf(fpsurf, "*,");
			else
				fprintf(fpsurf, "%.3f,", atom_ptr->wrapped_pos[1]);
			
			if( (atom_ptr->wrapped_pos[2] < 0.0005) || (atom_ptr->wrapped_pos[2] > -0.0005) )
				fprintf(fpsurf, "*,");
			else
				fprintf(fpsurf, "%.3f", atom_ptr->wrapped_pos[2]);
			
			fprintf(fpsurf, ";");
		
		}
	}

	fprintf(fpsurf, "!;");
	fflush(fpsurf);

}

double calctimediff(struct timeval a, struct timeval b) {
	return a.tv_sec - b.tv_sec + 1.0e-6 * (a.tv_usec - b.tv_usec);
}

int write_performance(int i, system_t *system) {

	static struct timeval current_time, last_time;

	char linebuf[MAXLINE];
	double sec_step;
	static int last_step;

	gettimeofday(&current_time,NULL);
	if(i > system->corrtime) {

		sec_step = calctimediff(current_time, last_time) / ((double)(i - last_step));

		if(system->ensemble == ENSEMBLE_UVT) {
			sprintf(linebuf, "OUTPUT: Grand Canonical Monte Carlo simulation running on %d cores\n", size);
			output(linebuf);
		} else {
			sprintf(linebuf, "OUTPUT: Canonical Monte Carlo simulation running on %d cores\n", size);
			output(linebuf);
		}
		sprintf(linebuf, "OUTPUT: Root collecting statistics at %s", ctime(&(current_time.tv_sec)));
		output(linebuf);
		sprintf(linebuf, "OUTPUT: Completed step %d/%d  (%.3f %%)\n", i, system->numsteps, (i/(double)(system->numsteps))*100);
		output(linebuf);
		sprintf(linebuf, "OUTPUT: %.3lf sec/step, ETA = %.3lf hrs\n", sec_step, sec_step*(system->numsteps - i)/3600.0);
		output(linebuf);

	}	

	last_step = i;
	last_time.tv_sec = current_time.tv_sec;
	last_time.tv_usec = current_time.tv_usec;

	return(0);

}

void write_observables(FILE *fp_energy, system_t * system, observables_t * observables, double core_temp) {

	fprintf(fp_energy, "%d %f %f %f %f %f %f %f %f %f %f %f", 
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
		observables->volume,
		core_temp);
	fprintf(fp_energy, "\n");
	fflush(fp_energy);
}

void write_observables_csv(FILE *fp_energy_csv, system_t * system, observables_t * observables, double core_temp) {

	fprintf(fp_energy_csv, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", 
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
		observables->volume,
		core_temp);
	fprintf(fp_energy_csv, "\n");
	fflush(fp_energy_csv);
}

/* output each molecular dipole (in debye) per line */
void write_dipole(system_t * system) {

	int p;
	FILE * fp;
	molecule_t *molecule_ptr;
	molecule_t *molecules = system->molecules;
	atom_t *atom_ptr;
	double dipole[3];

	//don't bother if we'd be writing to /dev/null
	if ( ! strncmp("/dev/null",system->dipole_output,9) ) return;
 	else fp = open_dipole_file(system);

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(p = 0; p < 3; p++) dipole[p] = 0;
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(p = 0; p < 3; p++)
				dipole[p] += atom_ptr->mu[p];
		}
		if(!molecule_ptr->frozen) fprintf(fp, "%f %f %f\n", dipole[0]/DEBYE2SKA, dipole[1]/DEBYE2SKA, dipole[2]/DEBYE2SKA);
	}
	fflush(fp);
	fclose(fp);
	return;
}

/* output the total molecular electrostatic field (in e/A) per line) */
void write_field(system_t * system) {

	int p;
	FILE * fp;
	molecule_t *molecule_ptr;
	molecule_t *molecules = system->molecules;
	atom_t *atom_ptr;
	double field[3];

	//don't bother if we'd be writing to /dev/null
	if ( ! strncmp("/dev/null",system->field_output,9) ) return;
 	else fp = open_field_file(system);

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(p = 0; p < 3; p++) field[p] = 0;
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(p = 0; p < 3; p++)
				field[p] += atom_ptr->ef_static[p] + atom_ptr->ef_induced[p];
		}
		if(!molecule_ptr->frozen) fprintf(fp, "%f %f %f\n", field[0]/E2REDUCED, field[1]/E2REDUCED, field[2]/E2REDUCED);
	}

	fflush(fp);
	fclose(fp);

	return;
}


int print_observables(system_t *system) {

	observables_t * o = system->observables;
	
	if(system->gwp) printf("OUTPUT: total energy = %.5lf eV\n", o->energy/EV2K);
	else printf("OUTPUT: potential energy = %.5lf K\n", o->energy);

	if(o->coulombic_energy != 0.0) {
		if(system->gwp) printf("OUTPUT: electrostatic energy = %.5lf eV\n", o->coulombic_energy/EV2K);
		else 	printf("OUTPUT: electrostatic energy = %.5lf K\n",	o->coulombic_energy);
	}

	if(o->rd_energy != 0.0)
		printf("OUTPUT: repulsion/dispersion energy = %.5lf K\n", o->rd_energy);

	if(o->polarization_energy != 0.0)
		printf("OUTPUT: polarization energy = %.5f K\n", o->polarization_energy);

#ifdef VDW
		printf("OUTPUT: (coupled-dipole) vdw energy = %.5f K\n", o->vdw_energy);
#endif

	printf("OUTPUT: N = %.5lf molecules\n", o->N);

	printf("OUTPUT: volume = %.5f A^3\n", system->pbc->volume);

	if(o->spin_ratio > 0.0) printf("OUTPUT: ortho spin ratio = %.5lf %%\n", o->spin_ratio*100.0);
	
	printf("\n");
	fflush(stdout);

	return(0);
}

int write_averages(system_t *system) {

	int i;

	avg_observables_t *averages;

	averages = system->avg_observables;

	if(averages->boltzmann_factor > 0.0)
		printf("OUTPUT: BF = %.5lg +- %.5lg\n", averages->boltzmann_factor, averages->boltzmann_factor_error);

	if(averages->acceptance_rate > 0.0) {
		printf("OUTPUT: AR = %.5lf (%.5lf I/ %.5lf R/ %.5lf D", 
			averages->acceptance_rate, averages->acceptance_rate_insert, 
			averages->acceptance_rate_remove, averages->acceptance_rate_displace);
		if(averages->acceptance_rate_adiabatic > 0.0) printf("/ %.5lf A", averages->acceptance_rate_adiabatic);
		if(averages->acceptance_rate_spinflip > 0.0) printf("/ %.5lf S", averages->acceptance_rate_spinflip);
		if(averages->acceptance_rate_volume > 0.0) printf("/ %.5lf V", averages->acceptance_rate_volume);
		if(averages->acceptance_rate_ptemp > 0.0) printf("/ %.5lf PT", averages->acceptance_rate_ptemp);
		printf(")\n");
	}

	//print node's current temperature if doing SA or PT
	if ( system->simulated_annealing )
		printf("OUTPUT: Simulated Annealing Temperature = %.5f K\n", system->temperature);

	if(averages->cavity_bias_probability > 0.0)
		printf("OUTPUT: Cavity bias probability = %.5f +- %.5f\n", 
			averages->cavity_bias_probability, averages->cavity_bias_probability_error);

	if(system->gwp) 
		printf("OUTPUT: total energy = %.5lf +- %.5lf eV\n", averages->energy/EV2K, averages->energy_error/EV2K);
	else
		printf("OUTPUT: potential energy = %.5lf +- %.5lf K\n", averages->energy, averages->energy_error);

	if(averages->coulombic_energy != 0.0) {
		if(system->gwp) 
			printf("OUTPUT: electrostatic energy = %.5lf +- %.5lf eV\n",
				averages->coulombic_energy/EV2K, averages->coulombic_energy_error/EV2K);
		else
			printf("OUTPUT: electrostatic energy = %.5lf +- %.5lf K\n",
				averages->coulombic_energy, averages->coulombic_energy_error);
	}

	if(averages->rd_energy != 0.0)
		printf("OUTPUT: repulsion/dispersion energy = %.5lf +- %.5lf K\n", 
			averages->rd_energy, averages->rd_energy_error);

	if(averages->polarization_energy != 0.0) {
		printf("OUTPUT: polarization energy = %.5f +- %.5f K", 
			averages->polarization_energy, averages->polarization_energy_error);
		if(averages->dipole_rrms_error != 0.0 && system->polar_rrms ) {
			printf(" (iterations = %.1f +- %.1f rrms = %e +- %e)", 
				averages->polarization_iterations, averages->polarization_iterations_error, 
				averages->dipole_rrms, averages->dipole_rrms_error);
		}
		else if(averages->polarization_iterations != 0.0) {
			printf(" (iterations = %.1f +- %.1f)", 
				averages->polarization_iterations, averages->polarization_iterations_error);
		}

		printf("\n");
	}

#ifdef VDW
		printf("OUTPUT: (coupled-dipole) vdw energy = %.5f +- %.5f K\n", 
			averages->vdw_energy, averages->vdw_energy_error);
#endif

	if(averages->kinetic_energy > 0.0) {
		if(system->gwp)
			printf("OUTPUT: kinetic energy = %.5lf +- %.5lf eV\n", 
				averages->kinetic_energy/EV2K, averages->kinetic_energy_error/EV2K);
		else
			printf("OUTPUT: kinetic energy = %.5lf +- %.5lf K\n", 
				averages->kinetic_energy, averages->kinetic_energy_error);

		printf("OUTPUT: kinetic temperature = %.5lf +- %.5lf K\n", 
			averages->temperature, averages->temperature_error);
	}

	printf("OUTPUT: N = %.5lf +- %.5lf molecules\n", averages->N, averages->N_error);

	if ( system->sorbateCount == 1 ) { //all based on calculations with assume only one type of sorbate
		printf("OUTPUT: density = %.5f +- %.5f g/cm^3\n", averages->density, averages->density_error);
		if(averages->pore_density != 0.0 && system->ensemble != ENSEMBLE_NPT)
			printf("OUTPUT: pore density = %.5f +- %.5f g/cm^3\n", averages->pore_density, averages->pore_density_error);
		if(averages->percent_wt > 0.0 ) {
			printf("OUTPUT: wt %% = %.5f +- %.5f %%\n", averages->percent_wt, averages->percent_wt_error);
			printf("OUTPUT: wt %% (ME) = %.5f +- %.5f %%\n", averages->percent_wt_me, averages->percent_wt_me_error);
		}
		if(averages->excess_ratio > 0.0)
			printf("OUTPUT: excess adsorption ratio = %.5f +- %.5f mg/g\n", averages->excess_ratio, averages->excess_ratio_error);
		if((averages->qst > 0.0) && isfinite(averages->qst))
			printf("OUTPUT: qst = %.5lf kJ/mol\n", averages->qst);
		if((averages->compressibility > 0.0) && isfinite(averages->compressibility)) {
			printf("OUTPUT: compressibility = %.6g +- %.6g atm^-1\n", averages->compressibility, averages->compressibility_error);
			printf("OUTPUT: bulk modulus = %.6g +- %.6g GPa\n", ATM2PASCALS*1.0e-9/averages->compressibility, 
				ATM2PASCALS*1.0e-9*averages->compressibility_error/averages->compressibility/averages->compressibility);
		}
	}

	if((averages->heat_capacity > 0.0) && (isfinite(averages->heat_capacity)))
		printf("OUTPUT: heat capacity = %.5g +- %.5g kJ/mol K\n", averages->heat_capacity, averages->heat_capacity_error);

	if(system->ensemble == ENSEMBLE_NPT || system->ensemble == ENSEMBLE_REPLAY)
		printf("OUTPUT: volume = %.5f +- %.5f A^3\n", averages->volume, averages->volume_error);

	if(averages->spin_ratio > 0.0) 
		printf("OUTPUT: ortho spin ratio = %.5lf +- %.5lf %%\n", averages->spin_ratio*100.0, averages->spin_ratio_error*100.0);
	
	if( system->sorbateCount > 1 ){

		for ( i=0; i < system->sorbateCount; i++ ) {

			printf( "OUTPUT: Stats for %s\n", system->sorbateInfo[i].id);
			printf( "             Average_N(%s)= %.5lf +- %.5lf\n", 
				system->sorbateInfo[i].id, system->sorbateGlobal[i].avgN, 
				system->sorbateGlobal[i].avgN_err);
			printf( "             Sorbed_Mass(%s)= %.5lf +- %.5lf g/mol\n",
				system->sorbateInfo[i].id, system->sorbateGlobal[i].avgN*system->sorbateInfo[i].mass, 
				system->sorbateGlobal[i].avgN_err*system->sorbateInfo[i].mass);
			printf( "             density(%s)= %.5le +- %.5le g/cm^3\n",
				system->sorbateInfo[i].id, system->sorbateGlobal[i].density, 
				system->sorbateGlobal[i].density_err);
			if ( system->observables->frozen_mass > 0 ) {
				printf( "             pore_density(%s)= %.5le +- %.5le g/cm^3\n",
					system->sorbateInfo[i].id, system->sorbateGlobal[i].pore_density, 
					system->sorbateGlobal[i].pore_density_err);
				printf( "             excess_ratio(%s)= %.5le +- %.5le g/cm^3\n",
					system->sorbateInfo[i].id, system->sorbateGlobal[i].excess_ratio, 
					system->sorbateGlobal[i].excess_ratio_err);
				printf( "             wt_%%(%s)= %.5lf +- %.5le %%\n",
					system->sorbateInfo[i].id, system->sorbateGlobal[i].percent_wt, 
					system->sorbateGlobal[i].percent_wt_err);
				printf( "             wt_%%(%s)(ME)= %.5lf +- %.5le %%\n",
					system->sorbateInfo[i].id, system->sorbateGlobal[i].percent_wt_me, 
					system->sorbateGlobal[i].percent_wt_me_err);
			}
			printf( "             Selectivity(%s)= %.4lf +- %.4lf\n", 
				system->sorbateInfo[i].id, system->sorbateGlobal[i].selectivity, 
				system->sorbateGlobal[i].selectivity_err);
		}
	}
	
	printf("\n");
	fflush(stdout);

	return(0);

}

void write_virial_output( system_t * system, double tmin, double tmax, double dt ) {
	double t;
	int i;
	FILE * fvirial = fopen(system->virial_output,"w");
	filecheck(fvirial,system->virial_output,WRITE);
		

	printf("### Start Virial Output ###\n");
	printf("#Temperature #B_2\n");
	fprintf(fvirial, "#Temperature #B_2\n");

	for ( i=0, t=tmin; t <= tmax; t+=dt ) {
		printf("%8.3lf %15.10lf\n", t, system->virial_coef[i]);
		fprintf(fvirial,"%8.3lf %15.10lf\n", t, system->virial_coef[i++]);
	}

	printf("### End Virial Output ###\n");
	fclose(fvirial);
	
	return;
}
