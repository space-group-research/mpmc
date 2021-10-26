#include <stdlib.h>
#include <mc.h>

double dddotprod(double* a, double* b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double didotprod(double* a, int* b) {
    return a[0] * (double)b[0] + a[1] * (double)b[1] + a[2] * (double)b[2];
}

int iidotprod(int* a, int* b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double min(double a, double b) {
    if (a > b)
        return b;
    else
        return a;
}

int is_singular(double **a, int n) {
    /*
     * tells us if a 2D matrix "a" of size n x n is singular.
     * Returns 1 if a is singular; else 0.
     */

    unsigned int j;
    int i, singular;
    double big, temp, *vv;

    singular = 0;

    vv = (double *) malloc(n * sizeof(double));
    memnullcheck(vv, n * sizeof(double), __LINE__ - 1, __FILE__);
    for (i = 0; i < n; i++) {
        big = 0.0;
        for (j = 0; j < n; j++)
            if ((temp = fabs(a[i][j])) > big) big = temp;
        /* note big cannot be zero */
        if (big == 0.0) {
            singular = 1;
            break;
        }
        vv[i] = 1.0 / big;
    }

    return singular;
}