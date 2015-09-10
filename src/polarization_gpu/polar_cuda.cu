/* 

@2010, Jonathan Belof
University of South Florida

*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <signal.h>
#include <math.h>

#include <structs.h>

#include <cuda.h>

#define BLOCKSIZE	64	/* optimal on Fermi */

//#define MAXATOMS	2048	/* maximum N value */
#define DAMP		2.1304f	/* exponential damping width */
#define MATRIXSIZE	16

/* max float value that we will allow for a squared distance */
/* so that when we calculate r^-5 we will not overflow, this is "infinity" */
#define MAXFVALUE	1.0e14f

extern "C" { void thole_field(system_t *); }


/* XXX TODO: multigpu support */
/* XXX TODO: check with uvt and test for mem leaks */
/* DONE: fixed arrays */
/* DONE: avoid i != j conditional through use of switch */
/* DONE: unroll all matrix op loops by hand to avoid conditional */
/* DONE: make use of constant memory for unit cell info */
/* DONE: add PK correction and preconditioning */


__constant__ float basis[MATRIXSIZE];		/* unit cell basis */
__constant__ float recip_basis[MATRIXSIZE];	/* recip-space basis */

/* kernel that performs a single iteration through the dipole field equations */
__global__ void thole_iterative_cuda(float4 *arg_ptr, int DATA_ARRAY_SIZE ) {

	int bid, tid, bsize, gsize;
	int i, j, k;
	float4 *pos, *estat, *mu_in, *mu_out, *eind_out;
	float4 dr, dri, img;
	float4 posi, posj, estati, mu, eind, mu_out_r;
	float r, r2, r3, r5;
	float damp = DAMP, damp2, damp3, damping_term1, damping_term2;	// exponential damping func 
	float expr;
	float sw; // switch to avoid conditional 
	__shared__ float4 sposj[BLOCKSIZE];
	__shared__ float4 smu[BLOCKSIZE];

	// dipole field tensor 
	// 0 = xx, 1 = xy, 2 = xz
	// 3 = yy, 4 = yz, 5 = zz
	float Tij[8];

	// get block and thead indices
	bid = blockIdx.x;
	tid = threadIdx.x;
	gsize = gridDim.x;
	bsize = blockDim.x;

	// this is the induced dipole that we are summing for 
	i = bid*bsize + tid;

	// set the arrays 
	pos = arg_ptr; arg_ptr += DATA_ARRAY_SIZE;
	estat = arg_ptr; arg_ptr += DATA_ARRAY_SIZE;
	mu_in = arg_ptr; arg_ptr += DATA_ARRAY_SIZE;
	mu_out = arg_ptr; arg_ptr += DATA_ARRAY_SIZE;
	eind_out = arg_ptr;

	// clear the induced field
	eind.x = eind.y = eind.z = eind.w = 0.0f;

	// locate estat fetch close to pos fetching
	estati = estat[i];
	posi = pos[i];

	// for each thread block
	for(j = 0; j < gsize; j++) {

		// fill the share mem dipole array for this block
		smu[tid] = mu_in[j*bsize+tid];
		sposj[tid] = pos[j*bsize+tid];

		// make sure our shared mem update is complete
		__syncthreads();

		// do work with the shared mem array
		for(k = 0; k < bsize; k++) {

			// local registers
			mu = smu[k];
			posj = sposj[k];

			// START MINIMUM IMAGE
			// get the particle displacement
			dr.x = posi.x - posj.x;
			dr.y = posi.y - posj.y;
			dr.z = posi.z - posj.z;

			// this switch will enforce both i != j and the (MAXATOM-N) null pairs 
			sw = (float)(!((int)dr.x) && !((int)dr.y) && !((int)dr.z));

			// matrix multiply with the inverse basis and round
			img.x = recip_basis[0]*dr.x + recip_basis[1]*dr.y + recip_basis[2]*dr.z;
			img.y = recip_basis[3]*dr.x + recip_basis[4]*dr.y + recip_basis[5]*dr.z;
			img.z = recip_basis[6]*dr.x + recip_basis[7]*dr.y + recip_basis[8]*dr.z;
			img.x = rintf(img.x);
			img.y = rintf(img.y);
			img.z = rintf(img.z);

			// matrix multiply to project back into our basis
			dri.x = basis[0]*img.x + basis[1]*img.y + basis[2]*img.z;
			dri.y = basis[3]*img.x + basis[4]*img.y + basis[5]*img.z;
			dri.z = basis[6]*img.x + basis[7]*img.y + basis[8]*img.z;

			// now correct the displacement
			dri.x = dr.x - dri.x;
			dri.y = dr.y - dri.y;
			dri.z = dr.z - dri.z;
			r2 = dri.x*dri.x + dri.y*dri.y + dri.z*dri.z;

			// various powers of r that we need
			r2 += sw*MAXFVALUE;
			r = sqrtf(r2);
			r3 = r2*r;
			r5 = r3*r2;
			r3 = 1.0f/r3;
			r5 = 1.0f/r5;
			// END MINIMUM IMAGE

			// damping terms
			damp2 = damp*damp;
			damp3 = damp2*damp;
			expr =  __expf(-damp*r);
			damping_term1 = 1.0f - expr*(0.5f*damp2*r2 + damp*r + 1.0f);
			damping_term2 = 1.0f - expr*(damp3*r*r2/6.0f + 0.5f*damp2*r2 + damp*r + 1.0f);

			// construct the Tij tensor field, unrolled by hand to avoid conditional on the diagonal terms
			damping_term1 *= r3;
			damping_term2 *= -3.0f*r5;
			// exploit symmetry
			// 0 = xx, 1 = xy, 2 = xz
			// 3 = yy, 4 = yz, 5 = zz
			Tij[0] = dri.x*dri.x*damping_term2 + damping_term1;
			Tij[1] = dri.x*dri.y*damping_term2;
			Tij[2] = dri.x*dri.z*damping_term2;
			Tij[3] = dri.y*dri.y*damping_term2 + damping_term1;
			Tij[4] = dri.y*dri.z*damping_term2;
			Tij[5] = dri.z*dri.z*damping_term2 + damping_term1;

			// contract dipole with the tensor
			eind.x -= Tij[0]*mu.x + Tij[1]*mu.y + Tij[2]*mu.z;
			eind.y -= Tij[1]*mu.x + Tij[3]*mu.y + Tij[4]*mu.z;
			eind.z -= Tij[2]*mu.x + Tij[4]*mu.y + Tij[5]*mu.z;

		} // end k

	} // end j


	// update the ith induced field vector and dipole in global mem
	mu_out_r.x = estati.w*(estati.x + eind.x);
	mu_out_r.y = estati.w*(estati.y + eind.y);
	mu_out_r.z = estati.w*(estati.z + eind.z);
	mu_out[i] = mu_out_r;
	eind_out[i] = eind;


}


extern "C" {

float polar_cuda(system_t *system) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	int current_iteration, i, p, q;	// local counters 
	int N, gridsize, blocksize;
	float *host_basis, *host_recip_basis;
	float4 *host_arg_ptr, *tmp_ptr, *host_pos, *host_estat, *host_mu_in, *host_mu_out, *host_eind_out;
	float4 *arg_ptr;	// on device
	float4 *host_echg_out;	// induced field change for Palmo-Krimm correction
	// polarization energy, to be returned
	float potential;
	cudaError_t error;
	cudaDeviceProp  prop;
	
	error = cudaGetDeviceProperties( &prop, 0 );
	if( error != cudaSuccess ) {
        printf( "POLAR_CUDA: GPU is reporting an error:\n%s\n", cudaGetErrorString( error ) );
		return(-1);
    }
    
    const int PER_BLOCK_MEM_REQUIREMENT = BLOCKSIZE * 5 * sizeof( float4 );
    const int MAXBLOCKS = prop.totalGlobalMem / PER_BLOCK_MEM_REQUIREMENT;
    const int MAXATOMS  = MAXBLOCKS * BLOCKSIZE; 


	// determine N 
	N = system->natoms;
	if(N >= MAXATOMS) {
		fprintf(stderr, "POLAR_CUDA: error, N = %d exceeds MAXATOMS = %d\n", N, MAXATOMS);
		return(-1);
	}

	// grid/block size
	gridsize = N/BLOCKSIZE + (N%BLOCKSIZE == 0?0:1);
	blocksize = BLOCKSIZE;
	const int DATA_ARRAY_SIZE = blocksize * gridsize; 

	// calculate the field vectors
	thole_field(system);

	// allocate temporary host arrays
	host_basis = (float *)calloc(MATRIXSIZE, sizeof(float));
	host_recip_basis = (float *)calloc(MATRIXSIZE, sizeof(float));
	host_arg_ptr = (float4 *)calloc(5*DATA_ARRAY_SIZE, sizeof(float4));
	tmp_ptr = host_arg_ptr;
	host_pos = tmp_ptr; tmp_ptr += DATA_ARRAY_SIZE;
	host_estat = tmp_ptr; tmp_ptr += DATA_ARRAY_SIZE;
	host_mu_in = tmp_ptr; tmp_ptr += DATA_ARRAY_SIZE;
	host_mu_out = tmp_ptr; tmp_ptr += DATA_ARRAY_SIZE;
	host_eind_out = tmp_ptr;
	if(system->polar_palmo) 
		host_echg_out = (float4 *)calloc(DATA_ARRAY_SIZE, sizeof(float4));


	// set some of the above arrays
	for(molecule_ptr = system->molecules, i = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next, i++) {

			host_pos[i].x = (float)atom_ptr->pos[0];
			host_pos[i].y = (float)atom_ptr->pos[1];
			host_pos[i].z = (float)atom_ptr->pos[2];

			host_mu_in[i].x = (float)system->polar_gamma*atom_ptr->polarizability*atom_ptr->ef_static[0];
			host_mu_in[i].y = (float)system->polar_gamma*atom_ptr->polarizability*atom_ptr->ef_static[1];
			host_mu_in[i].z = (float)system->polar_gamma*atom_ptr->polarizability*atom_ptr->ef_static[2];

			host_estat[i].x = (float)atom_ptr->ef_static[0];
			host_estat[i].y = (float)atom_ptr->ef_static[1];
			host_estat[i].z = (float)atom_ptr->ef_static[2];

			// stuff the polarizability in the 4th float of estat
			host_estat[i].w = (float)atom_ptr->polarizability;

		}
	}

	// copy over the basis matrix
	for(p = 0; p < 3; p++) {
		for(q = 0; q < 3; q++) {
			host_basis[p*3+q] = (float)system->pbc->basis[p][q];
			host_recip_basis[p*3+q] = (float)system->pbc->reciprocal_basis[p][q];
		}
	}

	// copy the array elements to constant memory
	cudaMemcpyToSymbol(basis, host_basis, MATRIXSIZE*sizeof(float), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(recip_basis, host_recip_basis, MATRIXSIZE*sizeof(float), 0, cudaMemcpyHostToDevice);


	// alloc aligned memory
	cudaMalloc((void **)&arg_ptr, 5*DATA_ARRAY_SIZE*sizeof(float4));
	cudaMemcpy(arg_ptr, host_arg_ptr, 5*DATA_ARRAY_SIZE*sizeof(float4), cudaMemcpyHostToDevice);

	// prefer 48KB of L1 cache rather than larger shared mem
	cudaFuncSetCacheConfig(thole_iterative_cuda, cudaFuncCachePreferL1);
	//cudaFuncSetCacheConfig(thole_iterative_cuda, cudaFuncCachePreferShared);

	// iterate until we are finished
	for(current_iteration = 0; current_iteration < system->polar_max_iter; current_iteration++) {

		// launch the kernel
		thole_iterative_cuda<<<gridsize, blocksize>>>(arg_ptr, DATA_ARRAY_SIZE);

		// feed dipoles back in for another pass
		cudaMemcpy((arg_ptr+2*DATA_ARRAY_SIZE), (arg_ptr+3*DATA_ARRAY_SIZE), DATA_ARRAY_SIZE*sizeof(float4), cudaMemcpyDeviceToDevice);

	}

	// check for errors after kernel finishes
	cudaThreadSynchronize();
	error = cudaGetLastError();
	if(error != cudaSuccess) {
		fprintf(stderr, "POLAR_CUDA: CUDA error: %s\n", cudaGetErrorString(error));
	}

	// copy back the results
	cudaMemcpy(host_mu_out, (arg_ptr+3*DATA_ARRAY_SIZE), DATA_ARRAY_SIZE*sizeof(float4), cudaMemcpyDeviceToHost);
	cudaMemcpy(host_eind_out, (arg_ptr+4*DATA_ARRAY_SIZE), DATA_ARRAY_SIZE*sizeof(float4), cudaMemcpyDeviceToHost);

	// if PK active, do one more iteration to get the change in induced field
	if(system->polar_palmo) {

		thole_iterative_cuda<<<gridsize, blocksize>>>(arg_ptr, DATA_ARRAY_SIZE );
		cudaMemcpy(host_echg_out, (arg_ptr+4*DATA_ARRAY_SIZE), DATA_ARRAY_SIZE*sizeof(float4), cudaMemcpyDeviceToHost);

	}

	// store the dipole vectors in the linked list on the host
	for(molecule_ptr = system->molecules, i = 0, potential = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next, i++) {

			// store dipole and induced field
			atom_ptr->mu[0] = (double)host_mu_out[i].x;
			atom_ptr->mu[1] = (double)host_mu_out[i].y;
			atom_ptr->mu[2] = (double)host_mu_out[i].z;
			atom_ptr->ef_induced[0] = (double)host_eind_out[i].x;
			atom_ptr->ef_induced[1] = (double)host_eind_out[i].y;
			atom_ptr->ef_induced[2] = (double)host_eind_out[i].z;

			// calculate the polarization energy as 1/2 mu*E
			potential += atom_ptr->mu[0]*atom_ptr->ef_static[0];
			potential += atom_ptr->mu[1]*atom_ptr->ef_static[1];
			potential += atom_ptr->mu[2]*atom_ptr->ef_static[2];

			if(system->polar_palmo) {

				atom_ptr->ef_induced_change[0] = (double)(host_echg_out[i].x - host_eind_out[i].x);
				atom_ptr->ef_induced_change[1] = (double)(host_echg_out[i].y - host_eind_out[i].y);
				atom_ptr->ef_induced_change[2] = (double)(host_echg_out[i].z - host_eind_out[i].z);
				potential += atom_ptr->mu[0]*atom_ptr->ef_induced_change[0];
				potential += atom_ptr->mu[1]*atom_ptr->ef_induced_change[1];
				potential += atom_ptr->mu[2]*atom_ptr->ef_induced_change[2];

			}
		}
	}
	potential *= -0.5;

	// free the host arrays
	free(host_arg_ptr);
	free(host_basis);
	free(host_recip_basis);

	// free the device array
	cudaFree(arg_ptr);

	// return the polarization energy
	return(potential);
}

} // extern "C"

