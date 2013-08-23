#include <stdlib.h>
#include <math.h>

double dddotprod ( double * a, double * b ) {
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

double didotprod ( double * a, int * b ) {
	return a[0]*(double)b[0] + a[1]*(double)b[1] + a[2]*(double)b[2];
}

int iidotprod ( int * a, int * b ) {
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

// get angle between two double vectors
double ddangle ( double * a, double * b ) {
	return acos( dddotprod(a,b) / sqrt( dddotprod(a,a) * dddotprod(b,b) ) );
}
