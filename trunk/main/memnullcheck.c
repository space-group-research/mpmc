#include <stdio.h>
#include <stdlib.h>



/* written to hopefully deal with some memory issues that have been plaguing very large
 * runs on PCN61. I believe that we may be running into memory fragmentation.
 * I am adding a function that will add an option defrag the thole_matricies every
 * X steps.
 * Also, as memory becomes scarce, it is neccessary to check the return values of
 * malloc, calloc and realloc to make sure we're not getting NULL's, particularly
 * for larger requests.
*/

// args are the pointer to check, the size of the memory requested and and
// integer that identifies the parent function
int memnullcheck ( void * ptr, int size, int parent ) {

	if ( ptr != NULL ) return 0;

	if ( ptr == NULL ) {
		//if the size is zero, then this is an expected result
		if ( size == 0 ) return 0;

		fprintf(stderr,"ERROR: Failed to allocate %d bytes.\n", size);
		fprintf(stderr,"ERROR: memnullcheck parent == %d.\n", parent);
		fprintf(stderr,"ERROR: the parent flag should help you identify the failed (m,re,c)alloc call.\n");

		exit(-1);
	}

}


/* parent flag list:
 * 0-14 = pairs.c
 * 15-20 = histogram.c
 * 21-33 = main.c
 * 34-36 = cleanup.c
 * 37-39 = cavity.c
 * 40-52 = mc_moves.c
 * 53    = pimc.c
 * 54-55 = surface.c
 * 56-57 = surf_fit.c
 * 58    = thole_field.c
 * 59-79 = thole_matrix.c
 * 80-81 = thole_iterative.c
 * 82-90 = rotational_eigenspectrum.c
 * 91-92 = rotational_potential.c
 * 93-131 = input.c
*/
