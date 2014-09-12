#ifndef MCH
#define MCH

extern int rank, size;

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <signal.h>
#include <limits.h>
#include <stdint.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef MPI
#include <mpi.h>
#endif /* MPI */

#include <defines.h>
#include <structs.h>
#include <function_prototypes.h>


#endif /*ifndef MCH*/
