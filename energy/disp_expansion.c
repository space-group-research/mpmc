#include <mc.h>

//Copyright 2013 Adam Hogan
//TO DO: long range correction, fix damping, improve speed

double disp_expansion(system_t *system)
{
	double potential = 0.0;

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	pair_t *pair_ptr;

	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if(pair_ptr->recalculate_energy) {

					/* make sure we're not excluded or beyond the cutoff */
					if(!((pair_ptr->rimg > system->pbc->cutoff) || pair_ptr->rd_excluded || pair_ptr->frozen)) {
						double r = pair_ptr->rimg;
						double r2 = r*r;
						double r4 = r2*r2;
						double r6 = r4*r2;
						double r8 = r6*r2;
						double r10 = r8*r2;
						double r12 = r10*r2;

						double c6 = pair_ptr->c6;
						double c8 = pair_ptr->c8;
						double c10 = pair_ptr->c10;
						double c12 = pair_ptr->c12;

						double repulsion = 0.0;

						repulsion = 315.7750382111558307123944638 * exp(-pair_ptr->epsilon*(r-pair_ptr->sigma)); // K = 10^-3 H ~= 316 K

						//pair_ptr->rd_energy = -tang_toennies_damping_function(pair_ptr,6)*c6/r6-tang_toennies_damping_function(pair_ptr,8)*c8/r8-tang_toennies_damping_function(pair_ptr,10)*c10/r10;
						pair_ptr->rd_energy = -c6/r6-c8/r8-c10/r10-c12/r12+repulsion;
					}

				}
				potential += pair_ptr->rd_energy;
			}
		}
	}

	return potential;
}

double disp_expansion_nopbc(system_t *system)
{
	double potential = 0.0;

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	pair_t *pair_ptr;

	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if(pair_ptr->recalculate_energy) {

					/* make sure we're not excluded or beyond the cutoff */
					if(!((pair_ptr->rimg > system->pbc->cutoff) || pair_ptr->rd_excluded || pair_ptr->frozen)) {
						double r = pair_ptr->r;
						double r2 = r*r;
						double r4 = r2*r2;
						double r6 = r4*r2;
						double r8 = r6*r2;
						double r10 = r8*r2;
						double r12 = r10*r2;

						double c6 = pair_ptr->c6;
						double c8 = pair_ptr->c8;
						double c10 = pair_ptr->c10;
						double c12 = pair_ptr->c12;

						double repulsion = 0.0;

						repulsion = 315.7750382111558307123944638 * exp(-pair_ptr->epsilon*(r-pair_ptr->sigma)); // K = 10^-3 H ~= 316 K

						//pair_ptr->rd_energy = -tang_toennies_damping_function(pair_ptr,6)*c6/r6-tang_toennies_damping_function(pair_ptr,8)*c8/r8-tang_toennies_damping_function(pair_ptr,10)*c10/r10;
						pair_ptr->rd_energy = -c6/r6-c8/r8-c10/r10-c12/r12+repulsion;
					}

				}
				potential += pair_ptr->rd_energy;
			}
		}
	}

	return potential;
}

double tang_toennies_factorial(int in)
{
	int i;
	double return_val = (double)in;
	if (in==0) return 1.0;
	for (i=1;i<in;i++)
	{
		return_val *= (double)(in-i);
	}
	return return_val;
}

double tang_toennies_damping_function(pair_t *pair_ptr, int order)
{
	int i;
	double damp = 0.0;
	for (i=0; i<order; i++)
	{
		damp += pow(pair_ptr->rimg/(2.0*pair_ptr->epsilon),(double)i)/tang_toennies_factorial(i);
	}
	return 1.0-exp(-pair_ptr->rimg/(2.0*pair_ptr->epsilon))*damp;
}
