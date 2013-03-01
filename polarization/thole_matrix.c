/* 

Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

void print_matrix(int N, double **matrix) {

	int i, j;

	printf("\n");
	for(i = 0; i < N; i++) {
		for(j = 0; j < N; j++) {
			printf("%.3f ", matrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");

}

void zero_out_amatrix ( system_t * system, int N ) {
	int i, j;
	/* zero out the matrix */
	for(i = 0; i < 3*N; i++)
		for(j = 0; j < 3*N; j++)
			system->A_matrix[i][j] = 0;
	return;
}
	
/* calculate the dipole field tensor */
void thole_amatrix(system_t *system) {

	int i, j, ii, jj, N, p, q;
	atom_t **atom_array;
	pair_t *pair_ptr;
	double damp1, damp2, wdamp1, wdamp2, v, s;
	double r, r2, ir3, ir5, ir;
	double rcut, rcut2, rcut3, rcut4, rcut5;
	rcut = system->pbc->cutoff;
	rcut2 = rcut*rcut; rcut3 = rcut2*rcut; rcut4 = rcut3*rcut; rcut5 = rcut4*rcut;
	double l, l2, l3, l4;
	l = system->polar_damp;
	l2 = l*l; l3 = l2*l; l4 = l3*l;
	double explr; //exp(-l*r)
	double explrcut = exp(-l*rcut);

	//array of atoms generated in pairs.c
	atom_array = system->atom_array;
	N = system->natoms;

	zero_out_amatrix(system,N);

	/* set the diagonal blocks */
	for(i = 0; i < N; i++) {
		ii = i*3;
		for(p = 0; p < 3; p++) {
			if(atom_array[i]->polarizability != 0.0)
				system->A_matrix[ii+p][ii+p] = 1.0/atom_array[i]->polarizability;
			else
				system->A_matrix[ii+p][ii+p] = MAXVALUE;
		}
	}

	/* calculate each Tij tensor component for each dipole pair */
	for(i = 0; i < (N - 1); i++) {
		ii = i*3;
		for(j = (i + 1), pair_ptr = atom_array[i]->pairs; j < N; j++, pair_ptr = pair_ptr->next) {
			jj = j*3;

			r = pair_ptr->rimg;
			r2 = r*r;

			/* inverse displacements */
			if(pair_ptr->rimg == 0.)
				ir3 = ir5 = MAXVALUE;
			else {
				ir = 1.0/r;
				ir3 = ir*ir*ir;
				ir5 = ir3*ir*ir;
			}

			//evaluate damping factors
			switch (system->damp_type) {
				case DAMPING_OFF:
					if ( pair_ptr->es_excluded )
						damp1 = damp2 = wdamp1 = wdamp2 = 0.0; 
					else 
						damp1 = damp2 = wdamp1 = wdamp2 = 1.0; 
					break;
				case DAMPING_LINEAR:
					s = l * pow(atom_array[i]->polarizability*atom_array[j]->polarizability, 1.0/6.0);
					v = r/s;
					if ( r < s ) {
						damp1 = (4.0 - 3.0*v)*v*v*v;
						damp2 = v*v*v*v;
					} else {
						damp1 = damp2 = 1.0;
					}
					break;
				case DAMPING_EXPONENTIAL:
					explr = exp(-l*r);
					damp1 = 1.0 - explr*(0.5*l2*r2 + l*r + 1.0);
					damp2 = damp1 - explr*(l3*r2*r/6.0);
					if ( system->polar_wolf_full ) { //subtract off damped interaction at r_cutoff
						wdamp1 = 1.0 - explrcut*(0.5*l2*rcut2+l*rcut+1.0);
						wdamp2 = wdamp1 - explrcut*(l3*rcut3/6.0);
					}
					break;
				default:
					error("error: something unexpected happened in thole_matrix.c");
			}

			/* build the tensor */
			for(p = 0; p < 3; p++) {
				for(q = 0; q < 3; q++) {

					system->A_matrix[ii+p][jj+q] = -3.0*pair_ptr->dimg[p]*pair_ptr->dimg[q]*damp2*ir5;
					if (system->polar_wolf_full) 
						system->A_matrix[ii+p][jj+q] -= -3.0*pair_ptr->dimg[p]*pair_ptr->dimg[q]*wdamp2*ir*ir/rcut3;

					/* additional diagonal term */
					if(p == q) {
						system->A_matrix[ii+p][jj+q] += damp1*ir3;
						if ( system->polar_wolf_full ) system->A_matrix[ii+p][jj+q] -= wdamp1/(rcut3);
					}	
				}
			}

			/* set the lower half of the tensor component */
			for(p = 0; p < 3; p++)
				for(q = 0; q < 3; q++)
					system->A_matrix[jj+p][ii+q] = system->A_matrix[ii+p][jj+q];

		} /* end j */
	} /* end i */

	return;
}

/* for uvt runs, resize the A (and B) matrices */
void thole_resize_matrices(system_t *system) {

	int i, N, dN, oldN;

	/* determine how the number of atoms has changed and realloc matrices */
	oldN = 3*system->checkpoint->thole_N_atom; //will be set to zero if first time called
	system->checkpoint->thole_N_atom = countNatoms(system);
	N = 3*system->checkpoint->thole_N_atom;
	dN = N-oldN;

	if(!dN) return;

	// grow A matricies by free/malloc (to prevent fragmentation)
	//free the A matrix
	for (i=0; i < oldN; i++) free(system->A_matrix[i]);
	free(system->A_matrix);

	//if not iterative, free the B matrix
	if (!system->polar_iterative) {
		for (i=0; i < oldN; i++) free(system->B_matrix[i]);
		free(system->B_matrix);
	}

	//(RE)allocate the A matrix
	system->A_matrix=calloc(N,sizeof(double*));
	memnullcheck(system->A_matrix,N*sizeof(double*), __LINE__-1, __FILE__);

	for (i=0; i< N; i++ ) {
		system->A_matrix[i]=malloc(N*sizeof(double));
		memnullcheck(system->A_matrix[i],N*sizeof(double), __LINE__-1, __FILE__);
	}

	//(RE)allocate the B matrix if not iterative
	if (!system->polar_iterative) {
		system->B_matrix=calloc(N,sizeof(double*));
		memnullcheck(system->B_matrix,N*sizeof(double*), __LINE__-1, __FILE__);
		for (i=0; i< N; i++ ) {
			system->B_matrix[i]=malloc(N*sizeof(double));
			memnullcheck(system->B_matrix[i],N*sizeof(double),__LINE__-1, __FILE__);
		}
	}

	return;
}


/* invert the A matrix */
void thole_bmatrix(system_t *system) {

	int N;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;

	/* count the number of atoms */
	N = 0;
	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
			++N;

	invert_matrix(3*N, system->A_matrix, system->B_matrix);

}

/* get the dipoles by vector matrix multiply */
void thole_bmatrix_dipoles(system_t *system) {

	int i, j, ii, p, N;
	atom_t **atom_array;
	double *mu_array, *field_array;

	atom_array = system->atom_array;
	N = system->natoms;

	/* allocate working arrays */
	mu_array = calloc(3*N, sizeof(double));
	memnullcheck(mu_array,3*N*sizeof(double), __LINE__-1, __FILE__);
	field_array = calloc(3*N, sizeof(double));
	memnullcheck(field_array,3*N*sizeof(double), __LINE__-1, __FILE__);
	
	/* copy the field in */
	for(i = 0; i < N; i++) {
		ii = i*3;

		for(p = 0; p < 3; p++)
			field_array[ii+p] = atom_array[i]->ef_static[p] + atom_array[i]->ef_static_self[p];

	}

	/* multiply the supervector with the B matrix */
	for(i = 0; i < 3*N; i++)
		for(j = 0; j < 3*N; j++)
			mu_array[i] += system->B_matrix[i][j]*field_array[j];

	/* copy the dipoles out */
	for(i = 0; i < N; i++) {
		ii = i*3;

		for(p = 0; p < 3; p++)
			atom_array[i]->mu[p] = mu_array[ii+p];

	}

	/* free the working arrays */
	free(mu_array);
	free(field_array);

}


/* numerical recipes routines for inverting a general matrix */

#define TINY	1.0e-20

void LU_decomp(double **a,int n,int *indx,double *d)
{
  int i,j,k,imax=0;
  double big,dum,sum,temp,*vv;

  vv=(double *)malloc(n*sizeof(double));
	memnullcheck(vv,n*sizeof(double), __LINE__-1, __FILE__);
  *d=1.0;
  for (i=0;i<n;i++) {
    big=0.0;
    for (j=0;j<n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
/* note big cannot be zero */
    if(big == 0.0) {
	printf("LU_decomp: matrix to invert can't be singular!\n");
	die(0);
    }
    vv[i]=1.0/big;
  }
  for (j=0;j<n;j++) {
    for (i=0;i<j;i++) {
      sum=a[i][j];
      for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<n;i++) {
      sum=a[i][j];
      for (k=0;k<j;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
        big=dum;
        imax=i;
      }
    }
    if (j != imax) {
      for (k=0;k<n;k++) {
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<n;i++) a[i][j] *= dum;
    }
  }
  free(vv);
}


void LU_bksb(double **a,int n,int *indx,double b[])
{
  int i,j,ii=-1,ip;
  double sum;

  for (i=0;i<n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii != -1) for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }

  for (i=n-1;i>=0;i--) {
    sum=b[i];
    for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}


/* invert a NxN matrix using routines from numerical recipes */
void invert_matrix(int n,double **a,double **ai) {

  int i,j,*indx;
  double *col,d;

  col  = (double *)malloc(n*sizeof(double));
	memnullcheck(col, n*sizeof(double), __LINE__-1, __FILE__);
  indx = (int *)malloc(n*sizeof(int));
	memnullcheck(indx, n*sizeof(int), __LINE__-1, __FILE__);

  LU_decomp(a,n,indx,&d);

  for(i=0;i<n;i++){
    for(j=0;j<n;j++) col[j] = 0.;
    col[i] = 1.;
    LU_bksb(a,n,indx,col);
    for(j=0;j<n;j++) ai[j][i] = col[j];
  }

  free(col);
	free(indx);

}


