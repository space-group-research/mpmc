/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

/* align the molecule with the x-axis for initialization */
void align_molecule(molecule_t *molecule) {

	atom_t *atom_ptr;
	double rotation_matrix[3][3];
	double com[3];
	int i, ii, n;
	double *new_coord_array;
	double theta, phi;

	/* count the number of atoms in the molecule, and allocate new coords array */
	for(atom_ptr = molecule->atoms, n = 0; atom_ptr; atom_ptr = atom_ptr->next) ++n;

	new_coord_array = calloc(n*3, sizeof(double));
	memnullcheck(new_coord_array,n*3*sizeof(double),91);

	/* save the com coordinate */
	com[0] = molecule->com[0];
	com[1] = molecule->com[1];
	com[2] = molecule->com[2];

	/* translate the molecule to the origin */
	for(atom_ptr = molecule->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
		atom_ptr->pos[0] -= com[0];
		atom_ptr->pos[1] -= com[1];
		atom_ptr->pos[2] -= com[2];
	}

	/* make sure that our atom chosen for finding the polar coords doesn't happen to be at the origin */
	for(atom_ptr = molecule->atoms; atom_ptr; atom_ptr = atom_ptr->next)
		if(!((atom_ptr->pos[0] == 0.0) && (atom_ptr->pos[1] == 0.0) && (atom_ptr->pos[2] == 0.0))) break;

	/* get the current polar coordinates */
	if(atom_ptr->pos[2] == 0.0)
		theta = 0.0;
	else
		theta = acos(atom_ptr->pos[2]/sqrt(atom_ptr->pos[0]*atom_ptr->pos[0] + atom_ptr->pos[1]*atom_ptr->pos[1] + atom_ptr->pos[2]*atom_ptr->pos[2]));
	if((atom_ptr->pos[1]/atom_ptr->pos[0]) > MAXDOUBLE)
		phi = atan(MAXDOUBLE);
	else {
		phi = atan(fabs(atom_ptr->pos[1]/atom_ptr->pos[0]));
		if((atom_ptr->pos[1]/atom_ptr->pos[0]) < 0.0) phi *= -1.0;
	}

	/* form the inverse rotation matrix and multiply */
	theta *= -1.0; phi *= -1.0;
	rotation_matrix[0][0] = cos(phi)*cos(theta);
	rotation_matrix[0][1] = -sin(phi);
	rotation_matrix[0][2] = cos(phi)*(-sin(theta));
	rotation_matrix[1][0] = sin(phi)*cos(theta);
	rotation_matrix[1][1] = cos(phi);
	rotation_matrix[1][2] = sin(phi)*(-sin(theta));
	rotation_matrix[2][0] = sin(theta);
	rotation_matrix[2][1] = 0;
	rotation_matrix[2][2] = cos(theta);
	for(atom_ptr = molecule->atoms, i = 0; atom_ptr; atom_ptr = atom_ptr->next, i++) {

		ii = i*3;
		new_coord_array[ii+0] = rotation_matrix[0][0]*atom_ptr->pos[0] + rotation_matrix[0][1]*atom_ptr->pos[1] + rotation_matrix[0][2]*atom_ptr->pos[2];
		new_coord_array[ii+1] = rotation_matrix[1][0]*atom_ptr->pos[0] + rotation_matrix[1][1]*atom_ptr->pos[1] + rotation_matrix[1][2]*atom_ptr->pos[2];
		new_coord_array[ii+2] = rotation_matrix[2][0]*atom_ptr->pos[0] + rotation_matrix[2][1]*atom_ptr->pos[1] + rotation_matrix[2][2]*atom_ptr->pos[2];

	}


	/* set the new coordinates and then translate back from the origin */
	for(atom_ptr = molecule->atoms, i = 0; atom_ptr; atom_ptr = atom_ptr->next, i++) {

		ii = i*3;
		atom_ptr->pos[0] = new_coord_array[ii+0];
		atom_ptr->pos[1] = new_coord_array[ii+1];
		atom_ptr->pos[2] = new_coord_array[ii+2];

		atom_ptr->pos[0] += com[0];
		atom_ptr->pos[1] += com[1];
		atom_ptr->pos[2] += com[2];

	}

	/* free our temporary array */
	free(new_coord_array);

}


/* rotate the selected molecule to an absolute spherical polar position of (theta, phi) */
void rotate_spherical(molecule_t *molecule, double theta, double phi) {

	atom_t *atom_ptr;
	double rotation_matrix[3][3];
	double com[3];
	int i, ii, n;
	double *new_coord_array;

	/* count the number of atoms in the molecule, and allocate new coords array */
	for(atom_ptr = molecule->atoms, n = 0; atom_ptr; atom_ptr = atom_ptr->next) ++n;

	new_coord_array = calloc(n*3, sizeof(double));
	memnullcheck(new_coord_array,n*3*sizeof(double),92);

	/* save the com coordinate */
	com[0] = molecule->com[0];
	com[1] = molecule->com[1];
	com[2] = molecule->com[2];

	/* translate the molecule to the origin */
	for(atom_ptr = molecule->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
		atom_ptr->pos[0] -= com[0];
		atom_ptr->pos[1] -= com[1];
		atom_ptr->pos[2] -= com[2];
	}

	/* form the rotation matrix and multiply */
	rotation_matrix[0][0] = cos(phi)*cos(theta);
	rotation_matrix[0][1] = -sin(phi);
	rotation_matrix[0][2] = cos(phi)*(-sin(theta));
	rotation_matrix[1][0] = sin(phi)*cos(theta);
	rotation_matrix[1][1] = cos(phi);
	rotation_matrix[1][2] = sin(phi)*(-sin(theta));
	rotation_matrix[2][0] = sin(theta);
	rotation_matrix[2][1] = 0;
	rotation_matrix[2][2] = cos(theta);
	for(atom_ptr = molecule->atoms, i = 0; atom_ptr; atom_ptr = atom_ptr->next, i++) {

		ii = i*3;
		new_coord_array[ii+0] = rotation_matrix[0][0]*atom_ptr->pos[0] + rotation_matrix[0][1]*atom_ptr->pos[1] + rotation_matrix[0][2]*atom_ptr->pos[2];
		new_coord_array[ii+1] = rotation_matrix[1][0]*atom_ptr->pos[0] + rotation_matrix[1][1]*atom_ptr->pos[1] + rotation_matrix[1][2]*atom_ptr->pos[2];
		new_coord_array[ii+2] = rotation_matrix[2][0]*atom_ptr->pos[0] + rotation_matrix[2][1]*atom_ptr->pos[1] + rotation_matrix[2][2]*atom_ptr->pos[2];

	}


	/* set the new coordinates and then translate back from the origin */
	for(atom_ptr = molecule->atoms, i = 0; atom_ptr; atom_ptr = atom_ptr->next, i++) {

		ii = i*3;
		atom_ptr->pos[0] = new_coord_array[ii+0];
		atom_ptr->pos[1] = new_coord_array[ii+1];
		atom_ptr->pos[2] = new_coord_array[ii+2];

		atom_ptr->pos[0] += com[0];
		atom_ptr->pos[1] += com[1];
		atom_ptr->pos[2] += com[2];

	}

	/* free our temporary array */
	free(new_coord_array);


}


/* hindered potential of sin(theta)^2 */
/* used for testing, compared with Curl, Hopkins, Pitzer,JCP (1968) 48:4064 */
double hindered_potential(double theta) {

	return(pow(sin(theta), 2.0));

}


/* put the molecule in a (theta,phi) orientation and return the energy */
double rotational_potential(system_t *system, molecule_t *molecule, double theta, double phi) {

	int q;
	molecule_t *molecule_backup;
	atom_t *atom_backup, *atom_ptr;
	double potential;


	if(system->quantum_rotation_hindered) {

		/* use a hindered rotor potential, useful for testing accuracy */
		potential = system->quantum_rotation_B*system->quantum_rotation_hindered_barrier*hindered_potential(theta);

	} else {

		/* backup the molecular coordinates */
		molecule_backup = copy_molecule(system, molecule);

		/* align the molecule with the x-axis initially */
		align_molecule(molecule);

		/* perform the rotation */
		rotate_spherical(molecule, theta, phi);

		/* get the energy for this configuration */
		potential = energy_no_observables(system);

		/* restore the original atomic coordinates */
		for(atom_backup = molecule_backup->atoms, atom_ptr = molecule->atoms; atom_backup; atom_backup = atom_backup->next, atom_ptr = atom_ptr->next)
			for(q = 0; q < 3; q++)
				atom_ptr->pos[q] = atom_backup->pos[q];

		/* restore the c.o.m. coordinates */
		for(q = 0; q < 3; q++) molecule->com[q] = molecule_backup->com[q];


		/* free the backup structure */
		free_molecule(system, molecule_backup);


	}

	return(potential);

}


