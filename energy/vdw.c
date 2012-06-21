/*

Keith McLaughlin, 2011
University of South Florida

Calculation of VDW using coupled dipole. A_matrix is already obtained
via polar.c. Cholesky decomposition is performed on generalized
eigenvalue problem A.u=omega^2 K.u to convert to standard eigenvalue 
problem C.v=omega^2 v. See polar-vdw.pdf for more details.

Frequencies are input in atomic units, and converted to s^(-1) during
energy calculation.


*/

#include <mc.h>
#include <math.h>
#define TWOoverHBAR 2.6184101e11 //K^-1 s^-1
#define cHBAR 7.63822291e-12 //Ks //HBAR is already taken to be in Js
#define halfHBAR 3.81911146e-12 //Ks
#define au2invsec 4.13412763705e16 //s^-1 a.u.^-1
#define FINITE_DIFF 0.01 //too small -> vdw calc noises becomes a problem

//kill MPI before quitting, when neccessary
void die( int code ){
#ifdef MPI
	MPI_Finalize();
#endif
	exit(code);
}

//check that the memory is properly allocated
void checknull ( void * ptr, char * ptr_name, int size ) {
	char errmsg[512];

	if ( size == 0 ) //if size == 0 then whatever
		return;

	if ( ptr == NULL ) {
		sprintf(errmsg, "vdw.c: error: malloc %s for %d bytes failed.\n", ptr_name, size);
		error(errmsg);
		die(-1);
	}
}

// STRUCT MATRIX, ALLOC, AND DELETE ************************************
struct mtx {
	int dim;
	double * val;
};

struct mtx * alloc_mtx ( int dim ) {
	//alloc matrix variable and set dim
	struct mtx * M = NULL;
	M = malloc(sizeof(struct mtx));
	checknull(M,"struct mtx * M", sizeof(struct mtx));
	M->dim=dim;
	//alloc matrix storage space
	M->val=calloc(dim*dim,sizeof(double));
	checknull(M->val,"struct mtx * M->val", dim*dim*sizeof(double));
	return M;
}

void free_mtx ( struct mtx * M ) {
	free(M->val);
	free(M);
	return;
}
// END STRUCT MATRIX STUFF *********************************************

//only run dsyev from LAPACK if compiled with -llapack, otherwise report an error
#ifdef VDW
//prototype for dsyev (LAPACK)
extern void dsyev_(char *, char *, int *, double *, int *, double *, double *, int *, int *);
#else
void dsyev_
(char * a, char *b, int * c, double * d, int * e, double * f, double * g, int * h, int * i) {
	error("ERROR: Not compiled with Linear Algebra VDW.\n");
	die(-1);
}
#endif

//can be used to test shit --- not actually used in any of the code
void print_mtx ( struct mtx * Cm ) {
	int iC, jC;

	printf("\n==============begin===================\n");
	for ( iC=0; iC<Cm->dim; iC++ ) {
		for ( jC=0; jC<Cm->dim; jC++ ) {
			if ( jC>iC ) continue;
			printf("%.1le ", (Cm->val)[iC+jC*(Cm->dim)]);
		}
		printf("\n");
	}
	printf("\n==============end=================\n");

	return;
}

//build C matrix for a given molecule/system, with atom indicies (offset)/3..(offset+dim)/3
struct mtx * build_M ( int dim, int offset, double ** Am, double * sqrtKinv ) {
	int i; //dummy
	int iA, jA; //Am indicies
	int iC, jC; //Cm indicies
	int nonzero; //non-zero col/rows in Am
	struct mtx * Cm; //return matrix Cm

	//count non-zero elements
	nonzero=0;
	for ( i=offset; i<dim+offset; i++ )
		if ( sqrtKinv[i] != 0 ) nonzero++;

	//allocate
	Cm = alloc_mtx(nonzero);

	//build lapack compatible matrix from Am[offset..dim, offset..dim]
	iC=jC=-1; //C index
	for ( iA=offset; iA<dim+offset; iA++ ) {
		if ( sqrtKinv[iA] == 0 ) continue; //suppress rows/cols full of zeros
		iC++; jC=-1;
		for ( jA=offset; jA<=iA; jA++ ) {
			if ( sqrtKinv[jA] == 0 ) continue; //suppress
			jC++;
			(Cm->val)[iC+jC*(Cm->dim)]=
				Am[iA][jA]*sqrtKinv[iA]*sqrtKinv[jA];
		}
	}

	return Cm;
}

void printevects ( struct mtx * M ) {
	int r,c;

	printf("%%vdw === Begin Eigenvectors ===\n");
	for ( r=0; r < (M->dim); r++ ) {
		for ( c=0; c < (M->dim); c++ ) {
			printf("%.2le ", (M->val)[r+c*M->dim]);
		}
		printf("\n");
	}
	printf("%%vdw=== End Eigenvectors ===\n");
}

/* LAPACK using 1D arrays for storing matricies.
	/ 0  3  6 \
	| 1  4  7 |		= 	[ 0 1 2 3 4 5 6 7 8 ]
	\ 2  5  8 /									*/
double * lapack_diag ( struct mtx * M, int jobtype ) {
	char job; //job type
	char uplo='L'; //operate on lower triagle
	double * work; //working space for dsyev
	int lwork; //size of work array
	int rval; //returned from dsyev_
	double * eigvals;

	//eigenvectors or no?
	if ( jobtype = 2 ) job='V';
	else job = 'N';

	if ( M->dim == 0 ) return NULL;

	//allocate eigenvalues array
	eigvals = malloc(M->dim*sizeof(double));
	checknull(eigvals,"double * eigvals",M->dim*sizeof(double));
	//optimize the size of work array
	lwork = -1;
	work = malloc(sizeof(double));
	checknull(work,"double * work",sizeof(double));
	dsyev_(&job, &uplo, &(M->dim), M->val, &(M->dim), eigvals, work, &lwork, &rval);
	//now optimize work array size is stored as work[0]
	lwork=(int)work[0];
	work = realloc(work,lwork*sizeof(double));
	checknull(work,"double * work",lwork*sizeof(double));
	//diagonalize
	dsyev_(&job, &uplo, &(M->dim), M->val, &(M->dim), eigvals, work, &lwork, &rval);

	if ( rval != 0 ) {
		fprintf(stderr,"error: LAPACK: dsyev returned error: %d\n", rval);
		die(-1);
	}

	free(work);

	return eigvals;
}

/* not needed unless T >> 300 
double wtanh ( double w, double T ) {
	if ( w < 1.0E-10 ) TWOoverHBAR*T/au2invsec; //from Taylor expansion
	if ( T == 0 ) return w;
	return w/tanh(halfHBAR*w*au2invsec/T);
}
*/

double eigen2energy ( double * eigvals, int dim, double temperature ) {
	int i;
	double rval=0;
	
	if ( eigvals == NULL ) return 0;

	for ( i=0; i<dim; i++ ) {
		if ( eigvals[i] < 0 ) eigvals[i]=0;
//		rval += wtanh(sqrt(eigvals[i]), temperature);
		rval += sqrt(eigvals[i]);
	}
	return rval;
}

//calculate energies for isolated molecules
//if we don't know it, calculate it and save the value
double calc_e_iso ( system_t * system, double * sqrtKinv, molecule_t * mptr ) {
	int nstart, nsize, curr_dimM;
	double e_iso; //total vdw energy of isolated molecules
	struct mtx * Cm_iso; //matrix Cm_isolated
	double * eigvals; //eigenvalues of Cm_cm
	molecule_t * molecule_ptr;
	atom_t * atom_ptr;

	nstart=nsize=0; //loop through each individual molecule
	for ( molecule_ptr = system->molecules; molecule_ptr; molecule_ptr=molecule_ptr->next ) {
		if ( molecule_ptr != mptr ) {  //count atoms then skip to next molecule
			for ( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) nstart++;
			continue;
		}

		//now that we've found the molecule of interest, count natoms, and calc energy
		for ( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) nsize++;

		//build matrix for calculation of vdw energy of isolated molecule
		Cm_iso = build_M(3*(nsize), 3*nstart, system->A_matrix, sqrtKinv);
		//diagonalize M and extract eigenvales -> calculate energy
		eigvals=lapack_diag(Cm_iso,1); //no eigenvectors
		e_iso=eigen2energy(eigvals,Cm_iso->dim,system->temperature);

		//free memory
		free(eigvals);
		free_mtx(Cm_iso);

		//convert a.u. -> s^-1 -> K
    return e_iso * au2invsec * halfHBAR ;
	}	

	//unmatched molecule
	return NAN; //we should never get here
}

//go through each molecule and determine the VDW energy associated with each isolated molecule
double sum_eiso_vdw ( system_t * system, double * sqrtKinv ) {

	double e_iso = 0;
	molecule_t * mp;
	atom_t * ap;
	vdw_t * vp;
	vdw_t * vpscan;

	//loop through molecules. if not known, calculate, store and count. otherwise just count.
	for ( mp = system->molecules; mp; mp=mp->next ) {
		for ( vp = system->vdw_eiso_info; vp != NULL; vp=vp->next ) { //loop through all vp's
			if ( strncmp(vp->mtype,mp->moleculetype,MAXLINE) == 0 ) {
					e_iso += vp->energy; //count energy
					break; //break out of vp loop. the current molecule is accounted for now. go to the next molecule
			}
			else continue; //not a match, check the next vp
		} //vp loop

		if ( vp == NULL ) { //if the molecule was unmatched, we need to grow the list
			// end of vp list and we haven't matched yet -> grow vdw_eiso_info
			// scan to the last non-NULL element
			if ( system->vdw_eiso_info == NULL ) {
				system->vdw_eiso_info = calloc(1,sizeof(vdw_t)); //allocate space
				checknull(system->vdw_eiso_info,"calloc vdw_t * vdw_eiso_info",sizeof(vdw_t));
				vpscan = system->vdw_eiso_info; //set scan pointer
			} else {
				for ( vpscan = system->vdw_eiso_info; vpscan->next != NULL; vpscan=vpscan->next );
				vpscan->next = calloc(1,sizeof(vdw_t)); //allocate space
				checknull(vpscan->next,"calloc vdw_t * vpscan->next", sizeof(vdw_t));
				vpscan = vpscan->next;
			} //done scanning and malloc'ing
		
			//set values
			strncpy(vpscan->mtype,mp->moleculetype,MAXLINE); //assign moleculetype
			vpscan->energy = calc_e_iso(system,sqrtKinv,mp); //assign energy
			if ( isfinite(vpscan->energy) == 0 ) { //if nan, then calc_e_iso failed
				fprintf(stderr,"VDW: ERROR: Problem in calc_e_iso.\n");
				die(-1);
			}
			//otherwise count the energy and move to the next molecule
			e_iso += vpscan->energy;

		} //vp==NULL
	} //mp loop	

	////all of this logic is actually really bad if we're doing surface fitting, since omega will change... :-(
	//free everything so we can recalc next step
	if ( system->ensemble == ENSEMBLE_SURF_FIT )  {
		free_vdw_eiso(system->vdw_eiso_info);
		system->vdw_eiso_info = NULL;
	}
			
	return e_iso;
}


//build the matrix K^(-1/2) -- see the PDF
double * getsqrtKinv ( system_t * system, int N ) {
	double * sqrtKinv;
	molecule_t * molecule_ptr;
	atom_t * atom_ptr;
	int i = 0;

	//malloc 3*N wastes an insignificant amount of memory, but saves us a lot of index management
	sqrtKinv = malloc(3*N*sizeof(double));
	checknull(sqrtKinv,"double * sqrtKinv",3*N*sizeof(double));

	for ( molecule_ptr = system->molecules; molecule_ptr; molecule_ptr=molecule_ptr->next ) {
		for ( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {
			//Seek through atoms, calculate sqrtKinv*
			sqrtKinv[i] = sqrt(atom_ptr->polarizability)*atom_ptr->omega;
			sqrtKinv[i+1] = sqrt(atom_ptr->polarizability)*atom_ptr->omega;
			sqrtKinv[i+2] = sqrt(atom_ptr->polarizability)*atom_ptr->omega;
			i+=3;
		}
	}
	
	return sqrtKinv;
}

// long-range correction
double lr_vdw_corr ( system_t * system ) {
	molecule_t * molecule_ptr;
	atom_t * atom_ptr;
	pair_t * pair_ptr;
	double w1, w2; //omegas
	double a1, a2; //alphas
	double cC; //leading coefficient to r^-6
	double corr = 0; //correction to the energy

	//skip if PBC isn't set-up
	if ( system->pbc->volume == 0 ) {
		fprintf(stderr,"VDW: PBC not set-up. Did you define your basis? Skipping LRC.\n");
		return 0;
	}

	for ( molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) {
		for ( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {
			for ( pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next ) {
					//skip if frozen
					if ( pair_ptr->frozen ) continue;
					//skip if same molecule  // don't do this... this DOES contribute to LRC
//					if ( molecule_ptr == pair_ptr->molecule ) continue;
					//fetch alphas and omegas
					a1=atom_ptr->polarizability;
					a2=pair_ptr->atom->polarizability;
					w1=atom_ptr->omega;
					w2=pair_ptr->atom->omega;
					if ( w1 == 0 || w2 == 0 || a1 == 0 || a2 == 0 ) continue; //no vdw energy
					// 3/4 hbar/k_B(Ks) omega(s^-1)  Ang^6
					cC=0.75 * cHBAR * sqrt(w1*w2) * au2invsec * a1 * a2;

					// long-range correction
					corr += -4.0/3.0 * M_PI * cC * pow(system->pbc->cutoff,-3) / system->pbc->volume;
			}
		}
	}

	return corr;

}

//calculate T matrix element for a particular separation
double e2body(system_t * system, atom_t * atom, pair_t * pair, double r) {
	double energy;
	double lr = system->polar_damp * r;
	double lr2 = lr*lr;
	double lr3 = lr*lr2;
	double Txx = pow(r,-3)*(-2.0+(0.5*lr3+lr2+2*lr+2)*exp(-lr));
	double Tyy = pow(r,-3)*(1-(0.5*lr2+lr+1)*exp(-lr));
	double * eigvals;
	struct mtx * M = alloc_mtx(6);
	
	//only the sub-diagonals are non-zero
	M->val[1]=M->val[2]=M->val[4]=M->val[5]=M->val[6]=M->val[8]=M->val[9]=M->val[11]=0;
	M->val[12]=M->val[13]=M->val[15]=M->val[16]=M->val[19]=M->val[20]=M->val[22]=M->val[23]=0;
	M->val[24]=M->val[26]=M->val[27]=M->val[29]=M->val[30]=M->val[31]=M->val[33]=M->val[34]=0;

	//true diagonals
	M->val[0]=M->val[7]=M->val[14]=(atom->omega)*(atom->omega);
	M->val[21]=M->val[28]=M->val[35]=(pair->atom->omega)*(pair->atom->omega);

	//sub-diagonals
	M->val[3]=M->val[18]=
		(atom->omega)*(pair->atom->omega)*sqrt(atom->polarizability*pair->atom->polarizability)*Txx;
	M->val[10]=M->val[17]=M->val[25]=M->val[32]=
		(atom->omega)*(pair->atom->omega)*sqrt(atom->polarizability*pair->atom->polarizability)*Tyy;

	eigvals=lapack_diag(M,1);
	energy = eigen2energy(eigvals, 6, system->temperature);

	//subtract energy of atoms at infinity
//	energy -= 3*wtanh(atom->omega, system->temperature);
	energy -= 3*atom->omega;
//	energy -= 3*wtanh(pair->atom->omega, system->temperature);
	energy -= 3*pair->atom->omega;

	free(eigvals);
	free_mtx(M);

  return energy * au2invsec * halfHBAR;
}

// feynman-hibbs correction - molecular pair finite differencing method
double fh_vdw_corr ( system_t * system ) {
	molecule_t * molecule_ptr;
	atom_t * atom_ptr;
	pair_t * pair_ptr;
	double rm; //reduced mass
	double E[5];  //energy at five points, used for finite differencing
	double dv, d2v, d3v, d4v; //derivatives
	double corr = 0; //correction to the energy
	double corr_single; //single vdw interaction energy
	double h = FINITE_DIFF; //small dr used for finite differencing //too small -> vdw calculation noise becomes a problem

	//for each pair
	for ( molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) {
		for ( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {
			for ( pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next ) {
				//skip if frozen
				if ( pair_ptr->frozen ) continue;
				//skip if they belong to the same molecule
				if ( molecule_ptr == pair_ptr->molecule ) continue;
				//skip if distance is greater than cutoff
				if ( pair_ptr->rimg > system->pbc->cutoff ) continue;
				//check if fh is non-zero
				if ( atom_ptr->polarizability == 0 || pair_ptr->atom->polarizability == 0 || 
					atom_ptr->omega == 0 || pair_ptr->atom->omega == 0 ) continue; //no vdw energy

				//calculate two-body energies
				E[0]=e2body(system,atom_ptr,pair_ptr,pair_ptr->rimg-h-h); //smaller r
				E[1]=e2body(system,atom_ptr,pair_ptr,pair_ptr->rimg-h); 
				E[2]=e2body(system,atom_ptr,pair_ptr,pair_ptr->rimg); //current r
				E[3]=e2body(system,atom_ptr,pair_ptr,pair_ptr->rimg+h); //larger r
				E[4]=e2body(system,atom_ptr,pair_ptr,pair_ptr->rimg+h+h);

				//derivatives (Numerical Methods Using Matlab 4E 2004 Mathews/Fink 6.2)
				dv = (E[3]-E[1])/(2.0*h);
				d2v = (E[3]-2.0*E[2]+E[1])/(h*h);
				d3v = (E[4]-2*E[3]+2*E[1]-E[0])/(2*pow(h,3));
				d4v = (E[4]-4*E[3]+6*E[2]-4*E[1]+E[0])/pow(h,4);
				
				// reduced mass
				rm=AMU2KG*(molecule_ptr->mass)*(pair_ptr->molecule->mass)/
					((molecule_ptr->mass)+(pair_ptr->molecule->mass));

				//2nd order correction
				corr_single = pow(METER2ANGSTROM, 2)*(HBAR*HBAR/(24.0*KB*system->temperature*rm))*(d2v + 2.0*dv/pair_ptr->rimg);
				//4th order correction
				if ( system->feynman_hibbs_order >= 4 )
					corr_single += pow(METER2ANGSTROM, 4)*(pow(HBAR, 4) /
						(1152.0*pow(KB*system->temperature*rm, 2))) *
						(15.0*dv/pow(pair_ptr->rimg, 3) + 4.0*d3v/pair_ptr->rimg + d4v);

				corr += corr_single;
			}
		}
	}

	return corr;
}

// feynman-hibbs using 2BE (shitty)
double fh_vdw_corr_2be ( system_t * system ) {
  molecule_t * molecule_ptr;
  atom_t * atom_ptr;
  pair_t * pair_ptr;
  double rm; //reduced mass
  double w1, w2; //omegas
  double a1, a2; //alphas
  double cC; //leading coefficient to r^-6
  double dv, d2v, d3v, d4v; //derivatives
  double corr = 0; //correction to the energy
  double corr_single; //single vdw interaction energy

  //for each pair
  for ( molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) { 
    for ( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) { 
      for ( pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next ) { 
        //skip if frozen
        if ( pair_ptr->frozen ) continue;
        //skip if they belong to the same molecule
        if ( molecule_ptr == pair_ptr->molecule ) continue;
        //skip if distance is greater than cutoff
        if ( pair_ptr->rimg > system->pbc->cutoff ) continue;
        //fetch alphas and omegas
        a1=atom_ptr->polarizability;
        a2=pair_ptr->atom->polarizability;
        w1=atom_ptr->omega;
        w2=pair_ptr->atom->omega;
        if ( w1 == 0 || w2 == 0 || a1 == 0 || a2 == 0 ) continue; //no vdw energy
        // 3/4 hbar/k_B(Ks) omega(s^-1)  Ang^6
        cC=0.75 * cHBAR * sqrt(w1*w2) * au2invsec * a1 * a2;
        // reduced mass
        rm=AMU2KG*(molecule_ptr->mass)*(pair_ptr->molecule->mass)/
          ((molecule_ptr->mass)+(pair_ptr->molecule->mass));

        //derivatives 
        dv = 6.0*cC*pow(pair_ptr->rimg,-7);
        d2v= dv * (-7.0)/pair_ptr->rimg;
        if ( system->feynman_hibbs_order >= 4 ) {
          d3v= d2v* (-8.0)/pair_ptr->rimg;
          d4v= d3v* (-9.0)/pair_ptr->rimg;
        }

        //2nd order correction
        corr_single = pow(METER2ANGSTROM, 2)*(HBAR*HBAR/(24.0*KB*system->temperature*rm))*(d2v + 2.0*dv/pair_ptr->rimg);
        //4th order correction
        if ( system->feynman_hibbs_order >= 4 )
          corr_single += pow(METER2ANGSTROM, 4)*(pow(HBAR, 4) /
            (1152.0*pow(KB*system->temperature*rm, 2))) *
            (15.0*dv/pow(pair_ptr->rimg, 3) + 4.0*d3v/pair_ptr->rimg + d4v);

        corr += corr_single;
      }
    }
  }

  return corr;
}

//with damping
double twobody(system_t * system) {
	molecule_t * molecule_ptr;
	atom_t * atom_ptr;
	pair_t * pair_ptr;
	double energy = 0;

	//for each pair
	for ( molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) {
		for ( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {
			for ( pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next ) {
				//skip if frozen
				if ( pair_ptr->frozen ) continue;
				//skip if they belong to the same molecule
				if ( molecule_ptr == pair_ptr->molecule ) continue;
				//skip if distance is greater than cutoff
				if ( pair_ptr->rimg > system->pbc->cutoff ) continue;
				//check if fh is non-zero
				if ( atom_ptr->polarizability == 0 || pair_ptr->atom->polarizability == 0 || 
					atom_ptr->omega == 0 || pair_ptr->atom->omega == 0 ) continue; //no vdw energy

				//calculate two-body energies
				energy+=e2body(system,atom_ptr,pair_ptr,pair_ptr->rimg);
			}
		}
	}

	return energy;
}
		

//returns interaction VDW energy
double vdw(system_t *system) {

	int N, dimC; //number of atoms, number of non-zero rows in C-Matrix
	double e_total, e_iso; //total energy, isolation energy (atoms @ infinity)
	double * sqrtKinv; //matrix K^(-1/2); cholesky decomposition of K
	double ** Am = system->A_matrix; //A_matrix
	struct mtx * Cm; //C_matrix (we use single pointer due to LAPACK requirements)
	double * eigvals; //eigenvales
	double fh_corr, lr_corr;

	N=getnatoms(system);

	//allocate arrays. sqrtKinv is a diagonal matrix. d,e are used for matrix diag.
	sqrtKinv = getsqrtKinv(system,N);

	//calculate energy vdw of isolated molecules
	e_iso = sum_eiso_vdw ( system, sqrtKinv );

	//Build the C_Matrix
	Cm = build_M (3*N, 0, Am, sqrtKinv);

	//setup and use lapack diagonalization routine dsyev_()
	eigvals = lapack_diag (Cm, system->polarvdw); //eigenvectors if system->polarvdw == 2
	if ( system->polarvdw == 2 )
		printevects(Cm);

	//return energy in inverse time (a.u.) units
	e_total = eigen2energy(eigvals, Cm->dim, system->temperature);
	e_total *= au2invsec * halfHBAR; //convert a.u. -> s^-1 -> K

	//vdw energy comparison
	if ( system->polarvdw == 3 )
		printf("VDW Two-Body | Many Body = %lf | %lf\n", twobody(system),e_total-e_iso);

	if ( system->feynman_hibbs ) {
		if ( system->vdw_fh_2be ) fh_corr = fh_vdw_corr_2be(system); //2be method
		else fh_corr = fh_vdw_corr(system); //mpfd
	}
	else fh_corr=0;

	if ( system->rd_lrc ) lr_corr = lr_vdw_corr(system);
	else lr_corr=0;

//cleanup and return
	free(sqrtKinv);
	free(eigvals);
	free_mtx(Cm);

	return e_total - e_iso + fh_corr + lr_corr;

}
