#include "cublas.h"
#include <cuda_runtime.h>
#include "structs.h"
#include "function_prototypes.h"
#include <stdio.h>
#include <mc.h>

#define THREADS 1024
#define MAXFVALUE 1.0e13f

__constant__ float basis[9];
__constant__ float recip_basis[9];

/**
 * Method uses exponential polarization regardless of method requested in input
 * script
 */
__global__ static void build_a(int N, float *A, const float damp, float3 *pos, float *pols) {
    int i = blockIdx.x, j;

    if (i >= N)
        return;

    float r, r2, r3, r5;
    float expr, damping_term1, damping_term2;
    float3 dr, dri, img;

    const float one_over_pol_i = 1.0 / pols[i];
    const float3 pos_i = pos[i];
    const float3 recip_basis_0 =
        make_float3(recip_basis[0], recip_basis[1], recip_basis[2]);
    const float3 recip_basis_1 =
        make_float3(recip_basis[3], recip_basis[4], recip_basis[5]);
    const float3 recip_basis_2 =
        make_float3(recip_basis[6], recip_basis[7], recip_basis[8]);
    const float3 basis_0 = make_float3(basis[0], basis[1], basis[2]);
    const float3 basis_1 = make_float3(basis[3], basis[4], basis[5]);
    const float3 basis_2 = make_float3(basis[6], basis[7], basis[8]);

    const float damp2 = damp * damp;
    const float damp3 = damp2 * damp;

    const int N_per_thread = int(N - 0.5) / THREADS + 1;
    const int threadid = threadIdx.x;
    const int threadid_plus_one = threadIdx.x + 1;
    for (j = threadid * N_per_thread; j < threadid_plus_one * N_per_thread && j < N; j++) {
        if (i == j) {
            A[9 * N * j + 3 * i] = one_over_pol_i;
            A[9 * N * j + 3 * i + 3 * N + 1] = one_over_pol_i;
            A[9 * N * j + 3 * i + 6 * N + 2] = one_over_pol_i;
            A[9 * N * j + 3 * i + 1] = 0.0;
            A[9 * N * j + 3 * i + 2] = 0.0;
            A[9 * N * j + 3 * i + 3 * N] = 0.0;
            A[9 * N * j + 3 * i + 3 * N + 2] = 0.0;
            A[9 * N * j + 3 * i + 6 * N] = 0.0;
            A[9 * N * j + 3 * i + 6 * N + 1] = 0.0;
        } else {
            // START MINIMUM IMAGE
            // get the particle displacement
            dr.x = pos_i.x - pos[j].x;
            dr.y = pos_i.y - pos[j].y;
            dr.z = pos_i.z - pos[j].z;

            // matrix multiply with the inverse basis and round
            img.x = recip_basis_0.x * dr.x + recip_basis_0.y * dr.y +
                recip_basis_0.z * dr.z;
            img.y = recip_basis_1.x * dr.x + recip_basis_1.y * dr.y +
                recip_basis_1.z * dr.z;
            img.z = recip_basis_2.x * dr.x + recip_basis_2.y * dr.y +
                recip_basis_2.z * dr.z;
            img.x = rintf(img.x);
            img.y = rintf(img.y);
            img.z = rintf(img.z);

            // matrix multiply to project back into our basis
            dri.x = basis_0.x * img.x + basis_0.y * img.y + basis_0.z * img.z;
            dri.y = basis_1.x * img.x + basis_1.y * img.y + basis_1.z * img.z;
            dri.z = basis_2.x * img.x + basis_2.y * img.y + basis_2.z * img.z;

            // now correct the displacement
            dri.x = dr.x - dri.x;
            dri.y = dr.y - dri.y;
            dri.z = dr.z - dri.z;
            r2 = dri.x * dri.x + dri.y * dri.y + dri.z * dri.z;

            // various powers of r that we need
            r = sqrtf(r2);
            r3 = r2 * r;
            r5 = r3 * r2;
            r3 = 1.0f / r3;
            r5 = 1.0f / r5;
            // END MINIMUM IMAGE

            // damping terms
            expr = __expf(-damp * r);
            damping_term1 = 1.0f - expr * (0.5f * damp2 * r2 + damp * r + 1.0f);
            damping_term2 = 1.0f - expr * (damp3 * r * r2 / 6.0f + 0.5f * damp2 * r2 +
                    damp * r + 1.0f);

            // construct the Tij tensor field, unrolled by hand to avoid conditional
            // on the diagonal terms
            damping_term1 *= r3;
            damping_term2 *= -3.0f * r5;

            // exploit symmetry
            A[9 * N * j + 3 * i] = dri.x * dri.x * damping_term2 + damping_term1;
            const float tmp1 = dri.x * dri.y * damping_term2;
            A[9 * N * j + 3 * i + 1] = tmp1;
            const float tmp2 = dri.x * dri.z * damping_term2;
            A[9 * N * j + 3 * i + 2] = tmp2;
            A[9 * N * j + 3 * i + 3 * N] = tmp1;
            A[9 * N * j + 3 * i + 3 * N + 1] =
                dri.y * dri.y * damping_term2 + damping_term1;
            const float tmp3 = dri.y * dri.z * damping_term2;
            A[9 * N * j + 3 * i + 3 * N + 2] = tmp3;
            A[9 * N * j + 3 * i + 6 * N] = tmp2;
            A[9 * N * j + 3 * i + 6 * N + 1] = tmp3;
            A[9 * N * j + 3 * i + 6 * N + 2] =
                dri.z * dri.z * damping_term2 + damping_term1;
        }
    }
    return;
}

__global__
void build_c(int matrix_size, int dim, float *A, float *pols, float *omegas, float *device_C_matrix) {
    int i = blockIdx.x;
    if (i >= matrix_size) return;

    const int N_per_thread = int(matrix_size - 0.5) / THREADS + 1;
    const int threadid = threadIdx.x;
    const int threadid_plus_one = threadIdx.x + 1;
    for (int j = threadid * N_per_thread; j < threadid_plus_one * N_per_thread && j < matrix_size; j++) {
        int row = j % dim;
        int col = j / dim;
        // Reverse aligned matrix because of fortran bindings
        device_C_matrix[row * 3 * dim / 3 + col] = A[j] * omegas[col / 3] * omegas[row / 3] * 
                                                   sqrt(pols[col / 3] * pols[row / 3]);
    }
}

__global__
void build_kinvsquared(int matrix_size, int dim, float *pols, float *omegas, float *device_invKsquared_matrix) {
    int i = blockIdx.x;
    if (i >= matrix_size) return;

    const int N_per_thread = int(matrix_size - 0.5) / THREADS + 1;
    const int threadid = threadIdx.x;
    const int threadid_plus_one = threadIdx.x + 1;
    for (int j = threadid * N_per_thread; j < threadid_plus_one * N_per_thread && j < matrix_size; j++) {
        int row = j % dim;
        int col = j / dim;
        if (row != col) {
            continue;
        }
        device_invKsquared_matrix[j] = 1 / (pols[row / 3] * omegas[row / 3] * omegas[row / 3]);
    }
}

__global__ static void print_a(int N, float *A) {
    printf("N: %d\n", N);
    for (int i = 0; i < 3 * 3 * N * N; i++) {
        printf("%8.5f ", A[i]);
        if ((i + 1) % (N * 3) == 0 && i != 0) {
            printf("\n");
        }
    }
    printf("\n");
}

extern "C" {
#include <stdlib.h>
#include <mc.h>

    double vdw_cuda(void *systemptr) {
        system_t *system = (system_t *)systemptr;
        int N = system->natoms;
        int matrix_size = 3 * 3 * N * N;

        double *C_matrix = (double *)malloc(matrix_size * sizeof(double));
        float *device_C_matrix;
        cudaMalloc(&device_C_matrix, matrix_size * sizeof(double));
        atom_t **atom_arr = system->atom_array;
        for (int i = 0; i < 3 * N; i++) {
            for (int j = 0; j < 3 * N; j++) {
                /* LAPACK using 1D arrays for storing matricies.
                   / 0  3  6 \
                   | 1  4  7 |		= 	[ 0 1 2 3 4 5 6 7 8 ]
                   \ 2  5  8 /									*/
                // C_matrix is stored in reverse order because fortran uses oppositely aligned arrays from C
                C_matrix[i * 3 * N + j] = system->A_matrix[j][i] * atom_arr[i / 3]->omega * atom_arr[j / 3]->omega *
                    sqrt(atom_arr[i / 3]->polarizability * atom_arr[j / 3]->polarizability);
            }
        }
        cudaMemcpy(device_C_matrix, C_matrix, matrix_size * sizeof(double), cudaMemcpyHostToDevice);


        molecule_t *molecule_ptr;
        atom_t *atom_ptr;
        int dim = 3 * N;
        float *host_pols, *host_basis, *host_recip_basis, *host_omegas;
        float3 *host_pos;
        host_pols = (float *)calloc(N, sizeof(float));
        host_pos = (float3 *)calloc(N, sizeof(float3));
        host_basis = (float *)calloc(9, sizeof(float));
        host_recip_basis = (float *)calloc(9, sizeof(float));
        host_omegas = (float *)calloc(N, sizeof(float));

        float *device_pols, *device_A, *device_omegas, *device_invKsquared_matrix;
        float3 *device_pos;
        cudaMalloc((void **)&device_pols, N * sizeof(float));
        cudaMalloc((void **)&device_pos, N * sizeof(float3));
        cudaMalloc((void **)&device_A, matrix_size * sizeof(float));
        cudaMalloc((void **)&device_omegas, N * sizeof(float));
        cudaMalloc((void **) &device_invKsquared_matrix, matrix_size * sizeof(float));


        // copy over the basis matrix
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                host_basis[i * 3 + j] = (float)system->pbc->basis[j][i];
                host_recip_basis[i * 3 + j] = (float)system->pbc->reciprocal_basis[j][i];
            }
        }
        cudaMemcpy(basis, host_basis, 9 * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(recip_basis, host_basis, 9 * sizeof(float), cudaMemcpyHostToDevice);

        int i;
        for (molecule_ptr = system->molecules, i = 0; molecule_ptr;
                molecule_ptr = molecule_ptr->next) {
            for (atom_ptr = molecule_ptr->atoms; atom_ptr;
                    atom_ptr = atom_ptr->next, i++) {
                host_pos[i].x = (float)atom_ptr->pos[0];
                host_pos[i].y = (float)atom_ptr->pos[1];
                host_pos[i].z = (float)atom_ptr->pos[2];
                host_pols[i] = (atom_ptr->polarizability == 0.0)
                    ? 1.0f / MAXFVALUE
                    : (float)atom_ptr->polarizability;
            }
        }
        cudaMemcpy(device_pos, host_pos, N * sizeof(float3), cudaMemcpyHostToDevice);
        cudaMemcpy(device_pols, host_pols, N * sizeof(float), cudaMemcpyHostToDevice);

        build_a<<<N, THREADS>>>(N, device_A, system->polar_damp, device_pos, device_pols);
        printf("here\n");

        print_a<<<N, THREADS>>>(N, device_A);
        exit(1);

        for (i = 0; i < N; i++) {
            host_omegas[i] = system->atom_array[i]->polarizability;
        }
        cudaMemcpy(device_omegas, host_omegas, N * sizeof(float), cudaMemcpyHostToDevice);
        build_c<<<N, THREADS>>>(matrix_size, dim, device_A, device_pols, device_omegas, device_C_matrix);
        build_kinvsquared<<<N, THREADS>>>(matrix_size, dim, device_pols, device_omegas, device_invKsquared_matrix);


        return 0.;
    }
}

