
extern int rank, size;

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <values.h>
#include <sys/time.h>
#include <time.h>
#include <signal.h>

#include <math.h>

#ifdef MPI
#include <mpi.h>
#endif /* MPI */

#include <physical_constants.h>
#include <structs.h>
#include <function_prototypes.h>

#define ENSEMBLE_UVT				1
#define ENSEMBLE_NVT				2
#define ENSEMBLE_SURF				3
#define ENSEMBLE_SURF_FIT			4
#define ENSEMBLE_NVE				5
#define ENSEMBLE_TE				6
#define ENSEMBLE_NPT				7

#define MOVETYPE_INSERT				1
#define MOVETYPE_REMOVE				2
#define MOVETYPE_DISPLACE			3
#define MOVETYPE_ADIABATIC			4
#define MOVETYPE_SPINFLIP			5
#define MOVETYPE_VOLUME			6

