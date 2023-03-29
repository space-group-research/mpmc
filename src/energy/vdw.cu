#include "cublas.h"
#include <cuda_runtime.h>
#include "structs.h"
#include "function_prototypes.h"
#include <stdio.h>

#define THREADS 1024

extern "C" {
    double vdw_cuda(void *systemptr) {
        system_t *system = (system_t *)systemptr;
        int N = system->natoms;
        float *A;
    }
}
