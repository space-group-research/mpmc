#include <stdio.h>
#include <stdlib.h>
#include <mc.h>

// include for gcc stack trace
#ifndef __WIN32__
#include <execinfo.h>
#define STDERR 2
#endif

/* written to hopefully deal with some memory issues that have been plaguing very large
 * runs on PCN61. I believe that we may be running into memory fragmentation.
 * I am adding a function that will add an option defrag the thole_matricies every
 * X steps.
 * Also, as memory becomes scarce, it is neccessary to check the return values of
 * malloc, calloc and realloc to make sure we're not getting NULL's, particularly
 * for larger requests.
*/

// args are the pointer to check, the size of the memory requested and
// the last two parameters should be called with the macros __LINE__ and __FILE__
// The first macro is replaced by the line number in source code, the latter
// is replaced by name of the source file.
// Subtracting one from the __LINE__ macro will indicate that the mem
// allocation failed on the line before the memnullcheck() call.
//
// Eg:
// int *ptr = malloc( sizeof(int) * 100 );
// memnullcheck( ptr, sizeof(int)*100, __LINE__-1, __FILE__);

int memnullcheck(void *ptr, int size, int line, char *file) {
    if (ptr != NULL) return 0;

    if (ptr == NULL) {
        //if the size is zero, then this is an expected result
        if (size == 0) return 0;

        fprintf(stderr,
                "ERROR: Failed to allocate %d bytes.\n", size);
        fprintf(stderr,
                "       Check %s:%d\n", file, line);
#ifndef __WIN32__
        // Print a stack trace
        // In gcc, must compile with -rdynamic option in order to see the function names
        fprintf(stderr,
                "STACK_TRACE:\n");
        void *array[20];
        size_t size;
        size = backtrace(array, 20);
        backtrace_symbols_fd(array, size, STDERR);

#endif
        //fprintf(stderr,"ERROR: memnullcheck parent == %d.\n", parent);
        //fprintf(stderr,"ERROR: the parent flag should help you identify the failed (m,re,c)alloc call.\n");

        die(-1);
    }
    return 0;
}

int filecheck(void *ptr, char *filename, int mode) {
    char linebuf[MAXLINE];

    if (ptr != NULL) return 0;
    if (ptr == NULL) {
        switch (mode) {
            case READ:
                sprintf(linebuf,
                        "ERROR: Failed to open file %s for read\n", filename);
                break;
            case WRITE:
                sprintf(linebuf,
                        "ERROR: Failed to open file %s for write\n", filename);
                break;
            case APPEND:
                sprintf(linebuf,
                        "ERROR: Failed to open file %s for append\n", filename);
                break;
            default:
                sprintf(linebuf,
                        "ERROR: Failed to open file %s\n", filename);
        }
        error(linebuf);
        die(-1);
    }

    return 0;
}
