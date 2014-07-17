#include <mc.h>

//Copyright 2013 Adam Hogan
//TO DO: long range correction

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
						double r14 = r12*r2;
						double r16 = r14*r2;
						double r18 = r16*r2;
						double r20 = r18*r2;

						double c6 = pair_ptr->c6;
						double c8 = pair_ptr->c8;
						double c10 = pair_ptr->c10;
						double c12 = pair_ptr->c12;
						double c14 = pair_ptr->c14;
						double c16 = pair_ptr->c16;
						double c18 = pair_ptr->c18;
						double c20 = pair_ptr->c20;

						double repulsion = 0.0;

						if (pair_ptr->epsilon!=0.0&&pair_ptr->sigma!=0.0)
							repulsion = 315.7750382111558307123944638 * exp(-pair_ptr->epsilon*(r-pair_ptr->sigma)); // K = 10^-3 H ~= 316 K

						if (system->damp_dispersion)
							pair_ptr->rd_energy = -tt_damping(6,pair_ptr->epsilon*r)*c6/r6-tt_damping(8,pair_ptr->epsilon*r)*c8/r8-tt_damping(10,pair_ptr->epsilon*r)*c10/r10-tt_damping(12,pair_ptr->epsilon*r)*c12/r12-tt_damping(14,pair_ptr->epsilon*r)*c14/r14-tt_damping(16,pair_ptr->epsilon*r)*c16/r16-tt_damping(18,pair_ptr->epsilon*r)*c18/r18-tt_damping(20,pair_ptr->epsilon*r)*c20/r20+repulsion;
						else
							pair_ptr->rd_energy = -c6/r6-c8/r8-c10/r10+repulsion;
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
						double r = pair_ptr->rimg;
						double r2 = r*r;
						double r4 = r2*r2;
						double r6 = r4*r2;
						double r8 = r6*r2;
						double r10 = r8*r2;
						double r12 = r10*r2;
						double r14 = r12*r2;
						double r16 = r14*r2;
						double r18 = r16*r2;
						double r20 = r18*r2;

						double c6 = pair_ptr->c6;
						double c8 = pair_ptr->c8;
						double c10 = pair_ptr->c10;
						double c12 = pair_ptr->c12;
						double c14 = pair_ptr->c14;
						double c16 = pair_ptr->c16;
						double c18 = pair_ptr->c18;
						double c20 = pair_ptr->c20;

						double repulsion = 0.0;

						if (pair_ptr->epsilon!=0.0&&pair_ptr->sigma!=0.0)
							repulsion = 315.7750382111558307123944638 * exp(-pair_ptr->epsilon*(r-pair_ptr->sigma)); // K = 10^-3 H ~= 316 K

						if (system->damp_dispersion)
							pair_ptr->rd_energy = -tt_damping(6,pair_ptr->epsilon*r)*c6/r6-tt_damping(8,pair_ptr->epsilon*r)*c8/r8-tt_damping(10,pair_ptr->epsilon*r)*c10/r10-tt_damping(12,pair_ptr->epsilon*r)*c12/r12-tt_damping(14,pair_ptr->epsilon*r)*c14/r14-tt_damping(16,pair_ptr->epsilon*r)*c16/r16-tt_damping(18,pair_ptr->epsilon*r)*c18/r18-tt_damping(20,pair_ptr->epsilon*r)*c20/r20+repulsion;
						else
							pair_ptr->rd_energy = -c6/r6-c8/r8-c10/r10+repulsion;
					}

				}
				potential += pair_ptr->rd_energy;
			}
		}
	}

	return potential;
}

long factorial(long n)
{
  if (n == 0)
    return 1;
  else
    return(n * factorial(n-1));
}

double tt_damping(int n, double br)
{
  double sum = 0.0;
  int i;
  for (i=0;i<=n;i++)
  {
    sum += pow(br,i)/(double)(factorial(i));
  }
  return 1.0-exp(-br)*sum;
}
