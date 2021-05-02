/* 

@2012, Edward W. Lowe, Jr.
Vanderbilt University

*/

#include <mc.h>

#include "CL/cl.h"

#define BLOCKSIZE	384	/* optimal on Fermi */

#define MAXATOMS	65536	/* maximum N value */
#define DAMP		2.1304f	/* exponential damping width */
#define MATRIXSIZE	16

#define VERBOSE 0

/* max float value that we will allow for a squared distance */
/* so that when we calculate r^-5 we will not overflow, this is "infinity" */
#define MAXFVALUE	1.0e14f

// Helper function to get error string
void CheckErr( cl_int error, char* CALL)
{
  static const char* s_error_string[] = {
	"CL_SUCCESS",                         //   0
	"CL_DEVICE_NOT_FOUND",                //  -1
	"CL_DEVICE_NOT_AVAILABLE",            //  -2
	"CL_COMPILER_NOT_AVAILABLE",          //  -3
	"CL_MEM_OBJECT_ALLOCATION_FAILURE",   //  -4
	"CL_OUT_OF_RESOURCES",                //  -5
	"CL_OUT_OF_HOST_MEMORY",              //  -6
	"CL_PROFILING_INFO_NOT_AVAILABLE",    //  -7
	"CL_MEM_COPY_OVERLAP",                //  -8
	"CL_IMAGE_FORMAT_MISMATCH",           //  -9
	"CL_IMAGE_FORMAT_NOT_SUPPORTED",      // -10
	"CL_BUILD_PROGRAM_FAILURE",           // -11
	"CL_MAP_FAILURE",                     // -12
	"",                                   // -13
	"",                                   // -14
	"",                                   // -15
	"",                                   // -16
	"",                                   // -17
	"",                                   // -18
	"",                                   // -19
	"",                                   // -20
	"",                                   // -21
	"",                                   // -22
	"",                                   // -23
	"",                                   // -24
	"",                                   // -25
	"",                                   // -26
	"",                                   // -27
	"",                                   // -28
	"",                                   // -29
	"CL_INVALID_VALUE",                   // -30
	"CL_INVALID_DEVICE_TYPE",             // -31
	"CL_INVALID_PLATFORM",                // -32
	"CL_INVALID_DEVICE",                  // -33
	"CL_INVALID_CONTEXT",                 // -34
	"CL_INVALID_QUEUE_PROPERTIES",        // -35
	"CL_INVALID_COMMAND_QUEUE",           // -36
	"CL_INVALID_HOST_PTR",                // -37
	"CL_INVALID_MEM_OBJECT",              // -38
	"CL_INVALID_IMAGE_FORMAT_DESCRIPTOR", // -39
	"CL_INVALID_IMAGE_SIZE",              // -40
	"CL_INVALID_SAMPLER",                 // -41
	"CL_INVALID_BINARY",                  // -42
	"CL_INVALID_BUILD_OPTIONS",           // -43
	"CL_INVALID_PROGRAM",                 // -44
	"CL_INVALID_PROGRAM_EXECUTABLE",      // -45
	"CL_INVALID_KERNEL_NAME",             // -46
	"CL_INVALID_KERNEL_DEFINITION",       // -47
	"CL_INVALID_KERNEL",                  // -48
	"CL_INVALID_ARG_INDEX",               // -49
	"CL_INVALID_ARG_VALUE",               // -50
	"CL_INVALID_ARG_SIZE",                // -51
	"CL_INVALID_KERNEL_ARGS",             // -52
	"CL_INVALID_WORK_DIMENSION",          // -53
	"CL_INVALID_WORK_GROUP_SIZE",         // -54
	"CL_INVALID_WORK_ITEM_SIZE",          // -55
	"CL_INVALID_GLOBAL_OFFSET",           // -56
	"CL_INVALID_EVENT_WAIT_LIST",         // -57
	"CL_INVALID_EVENT",                   // -58
	"CL_INVALID_OPERATION",               // -59
	"CL_INVALID_GL_OBJECT",               // -60
	"CL_INVALID_BUFFER_SIZE",             // -61
	"CL_INVALID_MIP_LEVEL",               // -62
	"CL_INVALID_GLOBAL_WORK_SIZE"         // -63
  };

  int index_err = -( error);

   // valid error
  if( index_err >= 1 || VERBOSE)
  {
	fprintf(stderr, "error for %s: %s\n\n", CALL, s_error_string[ index_err]);
  }

}


float polar_ocl(system_t *system) {
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	int current_iteration, i, p, q;	/* local counters */
	int N, gridsize, blocksize;
	float *host_basis, *host_recip_basis, *host_partial_potential;
	cl_float4 *host_arg_ptr, *tmp_ptr, *host_pos, *host_estat, *host_mu_in, *host_mu_out, *host_eind_out;
	cl_float4 *host_echg_out;	/* induced field change for Palmo-Krimm correction */
	cl_int *host_mol_id, *host_frozen;
	/* polarization energy, to be returned */
	float potential;
	cl_int error = CL_SUCCESS;
	size_t global[ 1];
	size_t local[ 1];

	/* determine N */
	N = system->natoms;
	if(N >= MAXATOMS) {
		fprintf(stderr, "POLAR_OCL: error, N = %d exceeds MAXATOMS = %d\n", N, MAXATOMS);
		return(-1);
	}

	blocksize = BLOCKSIZE;

	int rounded_size = ( ( ( N - 1) / BLOCKSIZE) + 1) * BLOCKSIZE;
	const int nr_blocks = rounded_size / BLOCKSIZE;

	/* calculate the field vectors */
//    thole_field(system);

	/* allocate temporary host arrays */
	host_basis = (float *)calloc(MATRIXSIZE, sizeof(float));
	host_recip_basis = (float *)calloc(MATRIXSIZE, sizeof(float));
	host_arg_ptr = (cl_float4 *)calloc(5*rounded_size, sizeof(cl_float4));
	host_mol_id = (cl_int*)calloc( rounded_size, sizeof( cl_int));
	host_frozen = (cl_int*)calloc( rounded_size, sizeof( cl_int));
	host_partial_potential = (cl_float*)calloc( nr_blocks, sizeof( cl_float));
	tmp_ptr = host_arg_ptr;
	host_pos = tmp_ptr; tmp_ptr += rounded_size;
	host_estat = tmp_ptr; tmp_ptr += rounded_size;
	host_mu_in = tmp_ptr; tmp_ptr += rounded_size;
	host_mu_out = tmp_ptr; tmp_ptr += rounded_size;
	host_eind_out = tmp_ptr;

	if(system->polar_palmo) host_echg_out = (cl_float4 *)calloc(rounded_size, sizeof(cl_float4));

//	fprintf(stdout, "cpu estat\n");

	/* set some of the above arrays */
	for(molecule_ptr = system->molecules, i = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next, i++) {

			host_pos[i].x = (float)atom_ptr->pos[0];
			host_pos[i].y = (float)atom_ptr->pos[1];
			host_pos[i].z = (float)atom_ptr->pos[2];

			/* put charge in pos for estat calc */
			host_pos[ i].w = (float)atom_ptr->charge;

//			fprintf(stdout,"atom %i:  \tx=%f\ty=%f\tz=%f\n", i, (float)atom_ptr->ef_static[0], (float)atom_ptr->ef_static[1], (float)atom_ptr->ef_static[2]); fflush(stdout);

			/* zeroing to calc estat in kernel */
//			host_estat[i].x = 0.0f;//(float)atom_ptr->ef_static[0];
//			host_estat[i].y = 0.0f;//(float)atom_ptr->ef_static[1];
//			host_estat[i].z = 0.0f;//(float)atom_ptr->ef_static[2];
			/* stuff the polarizability in the 4th float of estat */
			host_estat[i].w = (float)atom_ptr->polarizability;
			host_mol_id[i] = molecule_ptr->id;
			/* commenting out since we can also just calc mu_in in kernel */
//			host_mu_in[i].x = (float)system->polar_gamma*atom_ptr->polarizability*atom_ptr->ef_static[0];
//			host_mu_in[i].y = (float)system->polar_gamma*atom_ptr->polarizability*atom_ptr->ef_static[1];
//			host_mu_in[i].z = (float)system->polar_gamma*atom_ptr->polarizability*atom_ptr->ef_static[2];
			/* put frozen into mu_in */
			host_frozen[i] = atom_ptr->frozen;
		}
	}

	/* copy over the basis matrix */
	for(p = 0; p < 3; p++) {
		for(q = 0; q < 3; q++) {
			host_basis[p*3+q] = (float)system->pbc->basis[p][q];
			host_recip_basis[p*3+q] = (float)system->pbc->reciprocal_basis[p][q];
		}
	}
	cl_mem echg;
	/* copy the array elements to constant memory */
	cl_mem basis = clCreateBuffer( system->ocl->context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, MATRIXSIZE*sizeof(float), host_basis, &error);
	CheckErr( error, "clCreateBuffer");

	cl_mem recip_basis = clCreateBuffer( system->ocl->context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, MATRIXSIZE*sizeof(float), host_recip_basis, &error);
	CheckErr( error, "clCreateBuffer");

	cl_mem frozen = clCreateBuffer( system->ocl->context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, rounded_size*sizeof(cl_int), host_frozen, &error);
	CheckErr( error, "clCreateBuffer");

	cl_mem mol_id = clCreateBuffer( system->ocl->context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, rounded_size*sizeof(cl_int), host_mol_id, &error);
	CheckErr( error, "clCreateBuffer");

	/* alloc aligned memory */
	cl_mem arg_ptr = clCreateBuffer( system->ocl->context, CL_MEM_READ_WRITE, 5*rounded_size*sizeof(cl_float4), NULL, &error);
	CheckErr( error, "clCreateBuffer");

	cl_mem partial_potential = clCreateBuffer( system->ocl->context, CL_MEM_WRITE_ONLY, nr_blocks*sizeof(cl_float), NULL, &error);
	CheckErr( error, "clCreateBuffer");

	if(system->polar_palmo)
	{
		echg = clCreateBuffer( system->ocl->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, rounded_size*sizeof(cl_float4), host_echg_out, &error);
		CheckErr( error, "clCreateBuffer");

		error = clSetKernelArg( system->ocl->palmo_echg, 0, sizeof( cl_mem), &arg_ptr);
		CheckErr( error, "clSetKernelArg 0");
		error = clSetKernelArg( system->ocl->palmo_echg, 1, sizeof( cl_mem), &echg);
		CheckErr( error, "clSetKernelArg 1");
		error = clSetKernelArg( system->ocl->palmo_echg, 2, sizeof( cl_mem), &basis);
		CheckErr( error, "clSetKernelArg 2");
		error = clSetKernelArg( system->ocl->palmo_echg, 3, sizeof( cl_mem), &recip_basis);
		CheckErr( error, "clSetKernelArg 3");
		error = clSetKernelArg( system->ocl->palmo_echg, 4, sizeof( cl_uint), &rounded_size);
		CheckErr( error, "clSetKernelArg 4");
	}

	error = clEnqueueWriteBuffer( system->ocl->queue, arg_ptr, CL_FALSE, 0, 5*rounded_size*sizeof(cl_float4), host_arg_ptr, 0, NULL, NULL);
	CheckErr( error, "clEnqueueWriteBuffer");

	global[ 0] = rounded_size;
	local[ 0] = BLOCKSIZE;

	error = clSetKernelArg( system->ocl->kernel, 0, sizeof( cl_mem), &arg_ptr);
	CheckErr( error, "clSetKernelArg 0");
	error = clSetKernelArg( system->ocl->kernel, 1, sizeof( cl_mem), &basis);
	CheckErr( error, "clSetKernelArg 1");
	error = clSetKernelArg( system->ocl->kernel, 2, sizeof( cl_mem), &recip_basis);
	CheckErr( error, "clSetKernelArg 2");
	error = clSetKernelArg( system->ocl->kernel, 3, sizeof( cl_uint), &rounded_size);
	CheckErr( error, "clSetKernelArg 3");

	float gamma = (float)system->polar_gamma;
	float pbc_cut = (float)system->pbc->cutoff;

	/* do estat on gpu */
	error = clSetKernelArg( system->ocl->thole_estat, 0, sizeof( cl_mem), &arg_ptr);
	CheckErr( error, "clSetKernelArg 0");
	error = clSetKernelArg( system->ocl->thole_estat, 1, sizeof( cl_mem), &basis);
	CheckErr( error, "clSetKernelArg 1");
	error = clSetKernelArg( system->ocl->thole_estat, 2, sizeof( cl_mem), &recip_basis);
	CheckErr( error, "clSetKernelArg 2");
	error = clSetKernelArg( system->ocl->thole_estat, 3, sizeof( cl_uint), &rounded_size);
	CheckErr( error, "clSetKernelArg 3");
	error = clSetKernelArg( system->ocl->thole_estat, 4, sizeof( cl_float), &pbc_cut);
	CheckErr( error, "clSetKernelArg 4");
	error = clSetKernelArg( system->ocl->thole_estat, 5, sizeof( cl_float), &gamma);
	CheckErr( error, "clSetKernelArg 5");
	error = clSetKernelArg( system->ocl->thole_estat, 6, sizeof( cl_mem), &frozen);
	CheckErr( error, "clSetKernelArg 6");
	error = clSetKernelArg( system->ocl->thole_estat, 7, sizeof( cl_mem), &mol_id);
	CheckErr( error, "clSetKernelArg 7");

	/* reduction kernel args */
	error = clSetKernelArg( system->ocl->potential_reduction, 0, sizeof( cl_mem), &arg_ptr);
	CheckErr( error, "clSetKernelArg 0");
	error = clSetKernelArg( system->ocl->potential_reduction, 1, sizeof( cl_uint), &rounded_size);
	CheckErr( error, "clSetKernelArg 1");
	error = clSetKernelArg( system->ocl->potential_reduction, 2, sizeof( cl_mem), &partial_potential);
	CheckErr( error, "clSetKernelArg 2");

	error = clEnqueueNDRangeKernel( system->ocl->queue, system->ocl->thole_estat, 1, 0, global, local, 0, NULL, NULL);
	CheckErr( error, "clEnqueueNDRangeKernel");

//	error = clEnqueueReadBuffer( system->ocl->queue, arg_ptr, CL_TRUE, rounded_size*sizeof(cl_float4), rounded_size*sizeof(cl_float4), host_estat, 0, NULL, NULL);
//	CheckErr( error, "clEnqueueReadBuffer");
//
//	fprintf(stdout, "gpu estat\n");
//	for(molecule_ptr = system->molecules, i = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {
//		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next, i++) {
//			fprintf(stdout,"atom %i:  \tx=%f\ty=%f\tz=%f\n", i, host_estat[i].x, host_estat[i].y, host_estat[i].z); fflush(stdout);
//		}
//	}
//	 fflush(stdout);
	/* iterate until we are finished */
	for(current_iteration = 0; current_iteration < system->polar_max_iter; current_iteration++) {

		/* launch the kernel */
		error = clEnqueueNDRangeKernel( system->ocl->queue, system->ocl->kernel, 1, 0, global, local, 0, NULL, NULL);
		CheckErr( error, "clEnqueueNDRangeKernel");
//		error = clFinish( system->ocl->queue);
//		CheckErr( error, "clFinish");

		/* feed dipoles back in for another pass */
		error = clEnqueueCopyBuffer( system->ocl->queue, arg_ptr, arg_ptr, 3*rounded_size*sizeof(cl_float4), 2*rounded_size*sizeof(cl_float4), rounded_size*sizeof(cl_float4), 0, NULL, NULL);
		CheckErr( error, "clEnqueueCopyBuffer");
//		error = clFinish( system->ocl->queue);
//		CheckErr( error, "clFinish");

	}

	/* copy back the results */
	error = clEnqueueReadBuffer( system->ocl->queue, arg_ptr, CL_FALSE, rounded_size*sizeof(cl_float4), rounded_size*sizeof(cl_float4), host_estat, 0, NULL, NULL);
	CheckErr( error, "clEnqueueReadBuffer");

	error = clEnqueueReadBuffer( system->ocl->queue, arg_ptr, CL_FALSE, 3*rounded_size*sizeof(cl_float4), rounded_size*sizeof(cl_float4), host_mu_out, 0, NULL, NULL);
	CheckErr( error, "clEnqueueReadBuffer");

	error = clEnqueueReadBuffer( system->ocl->queue, arg_ptr, CL_TRUE, 4*rounded_size*sizeof(cl_float4), rounded_size*sizeof(cl_float4), host_eind_out, 0, NULL, NULL);
	CheckErr( error, "clEnqueueReadBuffer");

	/* if PK active, do one more iteration to get the change in induced field */
	if(system->polar_palmo) {

		error = clEnqueueNDRangeKernel( system->ocl->queue, system->ocl->palmo_echg, 1, 0, global, local, 0, NULL, NULL);
		CheckErr( error, "clEnqueueNDRangeKernel");

//		clFinish( system->ocl->queue);
//		CheckErr( error, "clFinis()");

		error = clEnqueueReadBuffer( system->ocl->queue, echg, CL_TRUE, 0, rounded_size*sizeof(cl_float4), host_echg_out, 0, NULL, NULL);
		CheckErr( error, "clEnqueueReadBuffer");

//		clFinish( system->ocl->queue);
//		CheckErr( error, "clFinis()");
	}

	error = clEnqueueNDRangeKernel( system->ocl->queue, system->ocl->potential_reduction, 1, 0, global, local, 0, NULL, NULL);
	CheckErr( error, "clEnqueueNDRangeKernel");

	error = clEnqueueReadBuffer( system->ocl->queue, partial_potential, CL_TRUE, 0, nr_blocks*sizeof(cl_float), host_partial_potential, 0, NULL, NULL);
	CheckErr( error, "clEnqueueReadBuffer");

	potential = 0;

	for( i = 0; i < nr_blocks; ++i)
	{
		potential += host_partial_potential[ i];
	}

	/* store the dipole vectors in the linked list on the host */
	for(molecule_ptr = system->molecules, i = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next, i++) {

			/* store dipole and induced field */
			atom_ptr->mu[0] = (double)host_mu_out[i].x;
			atom_ptr->mu[1] = (double)host_mu_out[i].y;
			atom_ptr->mu[2] = (double)host_mu_out[i].z;
			atom_ptr->ef_induced[0] = (double)host_eind_out[i].x;
			atom_ptr->ef_induced[1] = (double)host_eind_out[i].y;
			atom_ptr->ef_induced[2] = (double)host_eind_out[i].z;
			atom_ptr->ef_static[0] = (double)host_estat[i].x;
			atom_ptr->ef_static[1] = (double)host_estat[i].y;
			atom_ptr->ef_static[2] = (double)host_estat[i].z;


			/* calculate the polarization energy as 1/2 mu*E */
//			potential += atom_ptr->mu[0]*atom_ptr->ef_static[0];
//			potential += atom_ptr->mu[1]*atom_ptr->ef_static[1];
//			potential += atom_ptr->mu[2]*atom_ptr->ef_static[2];

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

	/* free the host arrays */
	free(host_arg_ptr);
	free(host_basis);
	free(host_recip_basis);

	/* free the device array */
	if(system->polar_palmo) clReleaseMemObject( echg);
	clReleaseMemObject(arg_ptr);
	clReleaseMemObject( basis);
	clReleaseMemObject( recip_basis);

	/* return the polarization energy */
	return(potential);
}

