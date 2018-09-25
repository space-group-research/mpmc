#include <stdlib.h>

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
