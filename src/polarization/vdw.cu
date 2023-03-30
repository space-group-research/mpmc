#include "cublas_v2.h"
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include "defines.h"
#include "structs.h"
#include "function_prototypes.h"
#include <stdio.h>
#include <mc.h>

#define THREADS 1024
#define MAXFVALUE 1.0e13f
#define halfHBAR 3.81911146e-12     //Ks
#define cHBAR 7.63822291e-12        //Ks //HBAR is already taken to be in Js
#define FINITE_DIFF 0.01            //too small -> vdw calc noises becomes a problem
#define TWOoverHBAR 2.6184101e11    //K^-1 s^-1

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
void build_c_matrix(int matrix_size, int dim, float *A, float *pols, float *omegas, float *device_C_matrix) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= matrix_size) return;

    int row = i % dim;
    int col = i / dim;
    // Reverse aligned matrix because of fortran bindings
    device_C_matrix[row * 3 * dim / 3 + col] = A[i] * omegas[col / 3] * omegas[row / 3] * 
                                                sqrt(pols[col / 3] * pols[row / 3]);
}

__global__
void print_arr(int arr_size, float *arr) {
    for (int i = 0; i < arr_size; i++) {
        printf("%.3le\n", arr[i]);
    }
}

__global__
void build_kinvsqrt(int matrix_size, int dim, float *pols, float *omegas, float *device_invKsqrt_matrix) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= matrix_size) return;

    int row = i % dim;
    int col = i / dim;
    if (row != col) {
        return;
    }
    device_invKsqrt_matrix[row * 3 * dim / 3 + col] = sqrt((pols[row / 3] * omegas[row / 3] * omegas[row / 3]));
}

__global__ static void print_matrix(int dim, float *A) {
    printf("N: %d\n", dim);
    for (int i = 0; i < dim * dim; i++) {
        printf("%.3le ", A[i]);
        if ((i + 1) % (dim) == 0 && i != 0) {
            printf("\n");
        }
    }
    printf("\n");
}


extern "C" {
#include <stdlib.h>
#include <mc.h>

    
    //only run dsyev from LAPACK if compiled with -llapack, otherwise report an error
#ifdef VDW
    //prototype for dsyev (LAPACK)
    extern void dsyev_(char *, char *, int *, double *, int *, double *, double *, int *, int *);
#else
    void dsyev_(char *a, char *b, int *c, double *d, int *e, double *f, double *g, int *h, int *i) {
        error(
            "ERROR: Not compiled with Linear Algebra VDW.\n");
        die(-1);
    }
#endif

    //calculate energies for isolated molecules
    //if we don't know it, calculate it and save the value
    static double calc_e_iso(system_t *system, molecule_t *mptr, float *device_A_matrix, float *device_pols,
            float *device_omegas) {
        int nstart, nsize;   // , curr_dimM;  (unused variable)
        double e_iso;        //total vdw energy of isolated molecules
        double *eigvals = (double *) malloc(3 * system->natoms);     //eigenvalues of Cm_cm
        molecule_t *molecule_ptr;
        atom_t *atom_ptr;

        nstart = nsize = 0;  //loop through each individual molecule
        for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
            if (molecule_ptr != mptr) {  //count atoms then skip to next molecule
                for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) nstart++;
                continue;
            }

            //now that we've found the molecule of interest, count natoms, and calc energy
            for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) nsize++;

            //build matrix for calculation of vdw energy of isolated molecule
            //Cm_iso = build_M(3 * (nsize), 3 * nstart, system->A_matrix, sqrtKinv);
            float *device_C_matrix;
            int dim = 3 * nsize;
            int matrix_size = dim * dim;
            cudaMalloc((void **) &device_C_matrix, matrix_size);
            int blocks = (matrix_size + THREADS - 1) / THREADS;
            build_c_matrix<<<blocks, THREADS>>>(matrix_size, dim, device_A_matrix, device_pols, device_omegas, device_C_matrix);
            //diagonalize M and extract eigenvales -> calculate energy

            int *devInfo;
            float *d_work;
            float *d_W;
            int lwork = 0;
            cudaMalloc((void **)&d_W, dim * sizeof(float));
            cudaMalloc((void **)&d_work, dim * sizeof(float));
            cudaMalloc((void **)&devInfo, sizeof(int));
            cusolverDnHandle_t cusolverH;
            cusolverDnCreate(&cusolverH);

            cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_NOVECTOR;
            cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;

            // Find optimal workspace size
            cusolverDnSsyevd_bufferSize(cusolverH, jobz, uplo, dim, device_C_matrix, dim, d_W, &lwork);
            cudaMalloc( (void **) &d_work, lwork * sizeof(float));
            // Solve for eigenvalues
            cusolverDnSsyevd(cusolverH, jobz, uplo, dim, device_C_matrix, dim, d_W, d_work, lwork, devInfo);
            cudaDeviceSynchronize();

            float *host_eigenvalues = (float *)malloc(dim * sizeof(float));
            cudaMemcpy(host_eigenvalues, d_W, dim * sizeof(float), cudaMemcpyDeviceToHost);
            for (int i = 0; i < dim; i++) {
                if (host_eigenvalues[i] < 0) host_eigenvalues[i] = 0;
                //printf("eigs[%d]: %e\n", i, host_eigenvalues[i]);
                e_iso += sqrt(host_eigenvalues[i]);
            }

            //free memory
            free(eigvals);
            cudaFree(device_C_matrix);

            //convert a.u. -> s^-1 -> K
            return e_iso * au2invseconds * halfHBAR;
        }

        //unmatched molecule
        return NAN;  //we should never get here
    }


    static double sum_eiso_vdw(system_t *system, float *device_A_matrix, float *device_pols, float *device_omegas) {
        char linebuf[MAXLINE];
        double e_iso = 0;
        molecule_t *mp;
        // atom_t * ap;  (unused variable)
        vdw_t *vp;
        vdw_t *vpscan;

        //loop through molecules. if not known, calculate, store and count. otherwise just count.
        for (mp = system->molecules; mp; mp = mp->next) {
            for (vp = system->vdw_eiso_info; vp != NULL; vp = vp->next) {  //loop through all vp's
                if (strncmp(vp->mtype, mp->moleculetype, MAXLINE) == 0) {
                    e_iso += vp->energy;  //count energy
                    break;                //break out of vp loop. the current molecule is accounted for now. go to the next molecule
                } else
                    continue;  //not a match, check the next vp
            }                  //vp loop

            if (vp == NULL) {  //if the molecule was unmatched, we need to grow the list
                               // end of vp list and we haven't matched yet -> grow vdw_eiso_info
                               // scan to the last non-NULL element
                if (system->vdw_eiso_info == NULL) {
                    system->vdw_eiso_info = (vdw_t *)calloc(1, sizeof(vdw_t));  //allocate space
                    vpscan = system->vdw_eiso_info;  //set scan pointer
                } else {
                    for (vpscan = system->vdw_eiso_info; vpscan->next != NULL; vpscan = vpscan->next)
                        ;
                    vpscan->next = (vdw_t *)calloc(1, sizeof(vdw_t));  //allocate space
                    vpscan = vpscan->next;
                }  //done scanning and malloc'ing

                //set values
                strncpy(vpscan->mtype, mp->moleculetype, MAXLINE);  //assign moleculetype
                vpscan->energy = calc_e_iso(system, mp, device_A_matrix, device_pols, device_omegas);  //assign energy
                if (isfinite(vpscan->energy) == 0) {                //if nan, then calc_e_iso failed
                    sprintf(linebuf, "VDW: Problem in calc_e_iso.\n");
                    exit(1);
                }
                //otherwise count the energy and move to the next molecule
                e_iso += vpscan->energy;

            }  //vp==NULL
        }      //mp loop

        ////all of this logic is actually really bad if we're doing surface fitting, since omega will change... :-(
        //free everything so we can recalc next step
        if (system->ensemble == ENSEMBLE_SURF_FIT) {
            system->vdw_eiso_info = NULL;
        }

        return e_iso;
    }

    //calculate T matrix element for a particular separation
    static double e2body(system_t *system, atom_t *atom, pair_t *pair, double r) {
        double energy = 0;
        double lr = system->polar_damp * r;
        double lr2 = lr * lr;
        double lr3 = lr * lr2;
        double Txx = pow(r, -3) * (-2.0 + (0.5 * lr3 + lr2 + 2 * lr + 2) * exp(-lr));
        double Tyy = pow(r, -3) * (1 - (0.5 * lr2 + lr + 1) * exp(-lr));
        double *eigvals = (double *) malloc(36 * sizeof(double));
        double *T_matrix = (double *) malloc(36 * sizeof(double));

        //only the sub-diagonals are non-zero
        T_matrix[1] = T_matrix[2] = T_matrix[4] = T_matrix[5] = T_matrix[6] = T_matrix[8] = T_matrix[9] = T_matrix[11] = 0;
        T_matrix[12] = T_matrix[13] = T_matrix[15] = T_matrix[16] = T_matrix[19] = T_matrix[20] = T_matrix[22] = T_matrix[23] = 0;
        T_matrix[24] = T_matrix[26] = T_matrix[27] = T_matrix[29] = T_matrix[30] = T_matrix[31] = T_matrix[33] = T_matrix[34] = 0;

        //true diagonals
        T_matrix[0] = T_matrix[7] = T_matrix[14] = (atom->omega) * (atom->omega);
        T_matrix[21] = T_matrix[28] = T_matrix[35] = (pair->atom->omega) * (pair->atom->omega);

        //sub-diagonals
        T_matrix[3] = T_matrix[18] =
            (atom->omega) * (pair->atom->omega) * sqrt(atom->polarizability * pair->atom->polarizability) * Txx;
        T_matrix[10] = T_matrix[17] = T_matrix[25] = T_matrix[32] =
            (atom->omega) * (pair->atom->omega) * sqrt(atom->polarizability * pair->atom->polarizability) * Tyy;

        //eigvals = lapack_diag(M, 1);
        //energy = eigen2energy(eigvals, 6, system->temperature);
        char job = 'N';
        char uplo = 'L';  //operate on lower triagle
        int workSize = -1;
        int rval = 0;
        int dim = 6;
        double *workArr = (double *) malloc(sizeof(double));
        dsyev_(&job, &uplo, &dim, T_matrix, &dim, eigvals, workArr, &workSize, &rval);
        //now optimize work array size is stored as work[0]
        workSize = (int)workArr[0];
        workArr = (double *) realloc(workArr, workSize * sizeof(double));
        //diagonalize
        dsyev_(&job, &uplo, &dim, T_matrix, &dim, eigvals, workArr, &workSize, &rval);

        //subtract energy of atoms at infinity
        //	energy -= 3*wtanh(atom->omega, system->temperature);
        energy -= 3 * atom->omega;
        //	energy -= 3*wtanh(pair->atom->omega, system->temperature);
        energy -= 3 * pair->atom->omega;
        for (int i = 0; i < dim; i++) {
            energy += sqrt(eigvals[i]);
        }

        free(eigvals);
        free(workArr);
        free(T_matrix);

        return energy * au2invseconds * halfHBAR;
    }

    //with damping
    static double twobody(system_t *system) {
        molecule_t *molecule_ptr;
        atom_t *atom_ptr;
        pair_t *pair_ptr;
        double energy = 0;

        //for each pair
        for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
            for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
                for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
                    //skip if frozen
                    if (pair_ptr->frozen) continue;
                    //skip if they belong to the same molecule
                    if (molecule_ptr == pair_ptr->molecule) continue;
                    //skip if distance is greater than cutoff
                    if (pair_ptr->rimg > system->pbc->cutoff) continue;
                    //check if fh is non-zero
                    if (atom_ptr->polarizability == 0 || pair_ptr->atom->polarizability == 0 ||
                        atom_ptr->omega == 0 || pair_ptr->atom->omega == 0) continue;  //no vdw energy

                    //calculate two-body energies
                    energy += e2body(system, atom_ptr, pair_ptr, pair_ptr->rimg);
                }
            }
        }

        return energy;
    }

    // feynman-hibbs using 2BE (shitty)
    static double fh_vdw_corr_2be(system_t *system) {
        molecule_t *molecule_ptr;
        atom_t *atom_ptr;
        pair_t *pair_ptr;
        double rm;                 //reduced mass
        double w1, w2;             //omegas
        double a1, a2;             //alphas
        double cC;                 //leading coefficient to r^-6
        double dv, d2v, d3v, d4v;  //derivatives
        double corr = 0;           //correction to the energy
        double corr_single;        //single vdw interaction energy

        //for each pair
        for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
            for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
                for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
                    //skip if frozen
                    if (pair_ptr->frozen) continue;
                    //skip if they belong to the same molecule
                    if (molecule_ptr == pair_ptr->molecule) continue;
                    //skip if distance is greater than cutoff
                    if (pair_ptr->rimg > system->pbc->cutoff) continue;
                    //fetch alphas and omegas
                    a1 = atom_ptr->polarizability;
                    a2 = pair_ptr->atom->polarizability;
                    w1 = atom_ptr->omega;
                    w2 = pair_ptr->atom->omega;
                    if (w1 == 0 || w2 == 0 || a1 == 0 || a2 == 0) continue;  //no vdw energy
                    // 3/4 hbar/k_B(Ks) omega(s^-1)  Ang^6
                    cC = 1.5 * cHBAR * w1 * w2 / (w1 + w2) * au2invseconds * a1 * a2;
                    // reduced mass
                    rm = AMU2KG * (molecule_ptr->mass) * (pair_ptr->molecule->mass) /
                        ((molecule_ptr->mass) + (pair_ptr->molecule->mass));

                    //derivatives
                    dv = 6.0 * cC * pow(pair_ptr->rimg, -7);
                    d2v = dv * (-7.0) / pair_ptr->rimg;
                    if (system->feynman_hibbs_order >= 4) {
                        d3v = d2v * (-8.0) / pair_ptr->rimg;
                        d4v = d3v * (-9.0) / pair_ptr->rimg;
                    }

                    //2nd order correction
                    corr_single = pow(METER2ANGSTROM, 2) * (HBAR * HBAR / (24.0 * KB * system->temperature * rm)) * (d2v + 2.0 * dv / pair_ptr->rimg);
                    //4th order correction
                    if (system->feynman_hibbs_order >= 4)
                        corr_single += pow(METER2ANGSTROM, 4) * (pow(HBAR, 4) / (1152.0 * pow(KB * system->temperature * rm, 2))) *
                                    (15.0 * dv / pow(pair_ptr->rimg, 3) + 4.0 * d3v / pair_ptr->rimg + d4v);

                    corr += corr_single;
                }
            }
        }

        return corr;
    }

    // feynman-hibbs correction - molecular pair finite differencing method
    static double fh_vdw_corr(system_t *system) {
        molecule_t *molecule_ptr;
        atom_t *atom_ptr;
        pair_t *pair_ptr;
        double rm;                 //reduced mass
        double E[5];               //energy at five points, used for finite differencing
        double dv, d2v, d3v, d4v;  //derivatives
        double corr = 0;           //correction to the energy
        double corr_single;        //single vdw interaction energy
        double h = FINITE_DIFF;    //small dr used for finite differencing //too small -> vdw calculation noise becomes a problem

        //for each pair
        for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
            for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
                for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
                    //skip if frozen
                    if (pair_ptr->frozen) continue;
                    //skip if they belong to the same molecule
                    if (molecule_ptr == pair_ptr->molecule) continue;
                    //skip if distance is greater than cutoff
                    if (pair_ptr->rimg > system->pbc->cutoff) continue;
                    //check if fh is non-zero
                    if (atom_ptr->polarizability == 0 || pair_ptr->atom->polarizability == 0 ||
                            atom_ptr->omega == 0 || pair_ptr->atom->omega == 0) continue;  //no vdw energy

                    //calculate two-body energies
                    E[0] = e2body(system, atom_ptr, pair_ptr, pair_ptr->rimg - h - h);  //smaller r
                    E[1] = e2body(system, atom_ptr, pair_ptr, pair_ptr->rimg - h);
                    E[2] = e2body(system, atom_ptr, pair_ptr, pair_ptr->rimg);      //current r
                    E[3] = e2body(system, atom_ptr, pair_ptr, pair_ptr->rimg + h);  //larger r
                    E[4] = e2body(system, atom_ptr, pair_ptr, pair_ptr->rimg + h + h);

                    //derivatives (Numerical Methods Using Matlab 4E 2004 Mathews/Fink 6.2)
                    dv = (E[3] - E[1]) / (2.0 * h);
                    d2v = (E[3] - 2.0 * E[2] + E[1]) / (h * h);
                    d3v = (E[4] - 2 * E[3] + 2 * E[1] - E[0]) / (2 * pow(h, 3));
                    d4v = (E[4] - 4 * E[3] + 6 * E[2] - 4 * E[1] + E[0]) / pow(h, 4);

                    // reduced mass
                    rm = AMU2KG * (molecule_ptr->mass) * (pair_ptr->molecule->mass) /
                        ((molecule_ptr->mass) + (pair_ptr->molecule->mass));

                    //2nd order correction
                    corr_single = pow(METER2ANGSTROM, 2) * (HBAR * HBAR / (24.0 * KB * system->temperature * rm)) * (d2v + 2.0 * dv / pair_ptr->rimg);
                    //4th order correction
                    if (system->feynman_hibbs_order >= 4)
                        corr_single += pow(METER2ANGSTROM, 4) * (pow(HBAR, 4) / (1152.0 * pow(KB * system->temperature * rm, 2))) *
                            (15.0 * dv / pow(pair_ptr->rimg, 3) + 4.0 * d3v / pair_ptr->rimg + d4v);

                    corr += corr_single;
                }
            }
        }

        return corr;
    }

    // long-range correction
    static double lr_vdw_corr(system_t *system) {
        molecule_t *molecule_ptr;
        atom_t *atom_ptr;
        pair_t *pair_ptr;
        double w1, w2;    //omegas
        double a1, a2;    //alphas
        double cC;        //leading coefficient to r^-6
        double corr = 0;  //correction to the energy

        //skip if PBC isn't set-up
        if (system->pbc->volume == 0) {
            fprintf(stderr, "VDW: PBC not set-up. Did you define your basis? Skipping LRC.\n");
            return 0;
        }

        for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
            for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
                for (pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
                    //skip if frozen
                    if (pair_ptr->frozen) continue;
                    //skip if same molecule  // don't do this... this DOES contribute to LRC
                    //					if ( molecule_ptr == pair_ptr->molecule ) continue;
                    //fetch alphas and omegas
                    a1 = atom_ptr->polarizability;
                    a2 = pair_ptr->atom->polarizability;
                    w1 = atom_ptr->omega;
                    w2 = pair_ptr->atom->omega;
                    if (w1 == 0 || w2 == 0 || a1 == 0 || a2 == 0) continue;  //no vdw energy
                    // 3/4 hbar/k_B(Ks) omega(s^-1)  Ang^6
                    cC = 1.5 * cHBAR * w1 * w2 / (w1 + w2) * au2invseconds * a1 * a2;

                    // long-range correction
                    corr += -4.0 / 3.0 * M_PI * cC * pow(system->pbc->cutoff, -3) / system->pbc->volume;
                }
            }
        }

        return corr;
    }


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

        float *device_pols, *device_A_matrix, *device_omegas, *device_invKsqrt_matrix;
        float3 *device_pos;
        cudaMalloc((void **)&device_pols, N * sizeof(float));
        cudaMalloc((void **)&device_pos, N * sizeof(float3));
        cudaMalloc((void **)&device_A_matrix, matrix_size * sizeof(float));
        cudaMalloc((void **)&device_omegas, N * sizeof(float));
        cudaMalloc((void **) &device_invKsqrt_matrix, matrix_size * sizeof(float));

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

        build_a<<<N, THREADS>>>(N, device_A_matrix, system->polar_damp, device_pos, device_pols);
        cudaDeviceSynchronize();

        for (i = 0; i < N; i++) {
            host_omegas[i] = system->atom_array[i]->omega;
        }
        cudaMemcpy(device_omegas, host_omegas, N * sizeof(float), cudaMemcpyHostToDevice);
        int blocks = (matrix_size + THREADS - 1) / THREADS;
        build_c_matrix<<<blocks, THREADS>>>(matrix_size, dim, device_A_matrix, device_pols, device_omegas, device_C_matrix);
        build_kinvsqrt<<<blocks, THREADS>>>(matrix_size, dim, device_pols, device_omegas, device_invKsqrt_matrix);
        /*
        cudaDeviceSynchronize();
        print_matrix<<<1, 1>>>(dim, device_C_matrix);
        printf("\n");
        cudaDeviceSynchronize();
        print_matrix<<<1, 1>>>(dim, device_invKsqrt_matrix);
        cudaDeviceSynchronize();
        */

        int *devInfo;
        float *d_work;
        float *d_W;
        int lwork = 0;
        cudaMalloc((void **)&d_W, dim * sizeof(float));
        cudaMalloc((void **)&d_work, dim * sizeof(float));
        cudaMalloc((void **)&devInfo, sizeof(int));
        cusolverDnHandle_t cusolverH;
        cusolverDnCreate(&cusolverH);

        cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_NOVECTOR;
        cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;

        // Find optimal workspace size
        cusolverDnSsyevd_bufferSize(cusolverH, jobz, uplo, dim, device_C_matrix, dim, d_W, &lwork);
        cudaMalloc( (void **) &d_work, lwork * sizeof(float));
        // Solve for eigenvalues
        cusolverDnSsyevd(cusolverH, jobz, uplo, dim, device_C_matrix, dim, d_W, d_work, lwork, devInfo);
        cudaDeviceSynchronize();

        float *host_eigenvalues = (float *)malloc(dim * sizeof(float));
        cudaMemcpy(host_eigenvalues, d_W, dim * sizeof(float), cudaMemcpyDeviceToHost);
        float e_total = 0;
        for (int i = 0; i < dim; i++) {
            if (host_eigenvalues[i] < 0) host_eigenvalues[i] = 0;
            //printf("eigs[%d]: %e\n", i, host_eigenvalues[i]);
            e_total += sqrt(host_eigenvalues[i]);
        }
        e_total *= au2invseconds * halfHBAR;

        //double e_iso = sum_eiso_vdw(system, device_A_matrix, device_pols, device_omegas);
        double e_iso = 0;

        //vdw energy comparison
        if (system->polarvdw == 3) {
            printf("VDW Two-Body | Many Body = %lf | %lf\n", twobody(system), e_total - e_iso);
        }

        double fh_corr, lr_corr;
        if (system->feynman_hibbs) {
            if (system->vdw_fh_2be)
                fh_corr = fh_vdw_corr_2be(system);  //2be method
            else
                fh_corr = fh_vdw_corr(system);  //mpfd
        } else
            fh_corr = 0;

        if (system->rd_lrc)
            lr_corr = lr_vdw_corr(system);
        else
            lr_corr = 0;


        double energy = e_total - e_iso + fh_corr + lr_corr;
        printf("etotal: %le\n", e_total);
        printf("e_iso: %le\n", e_iso);
        printf("fh_corr: %le\n", fh_corr);
        printf("lr_corr: %le\n", lr_corr);
        printf("vdw: %e\n", energy);
        return energy;
    }
}

