/****************************************************************/
/* Calculates the structure factor S(k) from a g(r)		*/
/*								*/
/* In latex, the structure factor is:				*/
/* S(k)=1+4\pi\frac{\rho}{k}\int_0^\infty r[g(r) - 1]\sin(kr)dr	*/
/*								*/
/* where rho is the atomic density, k is in distance units	*/
/*								*/
/* In order to eliminate noise form the S(k), the tail end of	*/
/* the g(r) needs to be extrapolated to a cleaner analytical	*/
/* form.  The form used here is:				*/
/* h(r) = \frac{A}{r} e^{-\frac{r}{r_0}} \sin(\frac{r}{r_1})	*/
/*								*/
/* where the parameters A, r0, r1 as well as the extrapolation	*/
/* point, are determined from fitting to the g(r)		*/
/*								*/
/* Fitting procedure:						*/
/* Set the extrapolation point to the third zero value in the	*/
/* h(r).  Set A to the height of the next peak.  Set r0 to a	*/
/* large value so that the exponential damping is minimial -	*/
/* then play with r1 until the phase of the sin wave is a close	*/
/* fit.  Then decrease r0 to increase the exponential damping	*/
/* until the fit is final - A may need minor tweaking at that	*/
/* point.							*/
/*								*/
/* usage: ./sf <-f> [density] [xp] [A] [r0] [r1] [datafile]	*/
/* <-f> - turns on fitting mode, h(r) extrapolated is output	*/
/* side-by-side with the g(r)-based h(r)			*/
/* [density] - density in units of molecules/A^3		*/
/* [datfile] - two column file containing the g(r) and domain	*/
/*								*/
/* compilation: gcc -o sf sf.c -lm				*/
/*								*/
/* @2006							*/
/* Jonathan Belof						*/
/* Space Research Group						*/
/* Department of Chemistry					*/
/* University of South Florida					*/
/****************************************************************/

#define RESOLUTION		0.001
#define MAXLINE			4000000
#define SK_LOWER_BOUND		0
#define SK_UPPER_BOUND		10
#define H_LOWER_BOUND		0
#define H_UPPER_BOUND		100

struct gor_t {
	double domain;
	double range;
};

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

/* generate a curve to fit A, r0 and r1 */
/* first column output is g(r) - 1, the second column is the extrapolated form */
void fit_gor(struct gor_t *gor, int num_gor, double extrapolation_point, double A, double r0, double r1) {

	int i, j;				/* working indices */
	double *h;				/* h(r) = g(r) - 1 */
	double *h_x;				/* h(r) = extrapolated */
	int lower_bound, upper_bound;		/* lower/upper bounds for the domain of h(r) */
	double resolution = RESOLUTION;

	lower_bound = (int)rint(((double)H_LOWER_BOUND) / resolution);
	upper_bound = (int)rint(((double)H_UPPER_BOUND) / resolution);
	h = calloc(upper_bound, sizeof(double));
	h_x = calloc(upper_bound, sizeof(double));

	printf("# extrapolation point = %f A\n", extrapolation_point);
	printf("# A = %f\n", A);
	printf("# r0 = %f\n", r0);
	printf("# r1 = %f\n", r1);

	/* generate h(r) = g(r) - 1 */
	for(i = 0; i < num_gor; i++) {
		for(j = lower_bound; j < upper_bound; j++) {
			if(fabs((gor + i)->domain - ((double)j)*resolution) < resolution)
				h[j] = (gor + i)->range - 1.0;
		}
	}

	/* generate h(r) = extrapolated */
	for(i = lower_bound ; ((double)i)*resolution < extrapolation_point; i++) {	/* use the g(r) data */
		h_x[i] = h[i];
	}
	for(; i < upper_bound; i++) {	/* start extrapolating */
		h_x[i] = (A / (((double)i)*resolution))*exp(-((double)i)*resolution/r0)*sin(((double)i)*resolution/r1);
	}

	/* output the two h(r)'s */
	for(i = lower_bound; i < upper_bound; i++)
		printf("%f %f %f\n", ((double)i)*resolution, h[i], h_x[i]);

}

/* generate the extrapolated h(r) and return a point to it */
double *make_h(struct gor_t *gor, int num_gor, double extrapolation_point, double A, double r0, double r1) {

	int i, j;				/* working indices */
	double *h;				/* h(r) = g(r) - 1 */
	double *h_x;				/* h(r) = extrapolated */
	int lower_bound, upper_bound;		/* lower/upper bounds for the domain of h(r) */
	double resolution = RESOLUTION;

	lower_bound = (int)rint(((double)H_LOWER_BOUND) / resolution);
	upper_bound = (int)rint(((double)H_UPPER_BOUND) / resolution);
	h = calloc(upper_bound, sizeof(double));
	h_x = calloc(upper_bound, sizeof(double));

	printf("# extrapolation point = %f A\n", extrapolation_point);
	printf("# A = %f\n", A);
	printf("# r0 = %f\n", r0);
	printf("# r1 = %f\n", r1);

	/* generate h(r) = g(r) - 1 */
	for(i = 0; i < num_gor; i++) {
		for(j = lower_bound; j < upper_bound; j++) {
			if(fabs((gor + i)->domain - ((double)j)*resolution) < resolution)
				h[j] = (gor + i)->range - 1.0;
		}
	}

	/* generate h(r) = extrapolated */
	for(i = lower_bound ; ((double)i)*resolution < extrapolation_point; i++) {	/* use the g(r) data */
		h_x[i] = h[i];
	}
	for(; i < upper_bound; i++) {	/* start extrapolating */
		h_x[i] = (A / (((double)i)*resolution))*exp(-((double)i)*resolution/r0)*sin(((double)i)*resolution/r1);
	}

	free(h);

	return(h_x);

}

/* this function calculates the structure factor S(k) */
void structure_factor(struct gor_t *gor, int num_gor, double density, double extrapolation_point, double A, double r0, double r1) {

	int i, j;
	double *h;
	int sk_lower_bound, sk_upper_bound;
	int h_lower_bound, h_upper_bound;
	double integral, k, r;
	double resolution = RESOLUTION;

	/* get the integration boundaries */
	/* S(k) boundaries */
	sk_lower_bound = (int)rint(((double)SK_LOWER_BOUND) / resolution);
	/* advance the lower_bound so that we don't divide by zero */
	++sk_lower_bound;
	sk_upper_bound = (int)rint(((double)SK_UPPER_BOUND) / resolution);
	/* h(r) boundaries */
	h_lower_bound = (int)rint(((double)H_LOWER_BOUND) / resolution);
	++h_lower_bound;
	h_upper_bound = (int)rint(((double)H_UPPER_BOUND) / resolution);

	/* get the extrapolated h(r) */
	h = make_h(gor, num_gor, extrapolation_point, A, r0, r1);

	/* get the structure factor S(k) */
	for(i = sk_lower_bound; i < sk_upper_bound; i++) {

		k = ((double)i)*resolution;
		/* integrate h(r)*r*sin(kr) trapezoidally */
		integral = 0.0;
		for(j = h_lower_bound ; j < h_upper_bound; j++) {
			r = ((double)j)*resolution;
			integral += h[j]*r*sin(k*r);
		}
		/* finish the integral */
		integral -= 0.5*(h[h_lower_bound] + h[h_upper_bound - 1]);
		integral *= resolution;

		/* multiply coefficients */
		integral *= (4.0*M_PI*density / k);
		integral += 1.0;

		/* output the current S(k) value */
		printf("%f %f\n", k, integral);
	}

	free(h);
}


void usage(char *progname) {

	fprintf(stderr, "usage: %s <-f> [density] [xp] [A] [r0] [r1] [datafile]\n", progname);
	fprintf(stderr, "<-f> - turns on fitting mode, h(r) extrapolated is output\n");
	fprintf(stderr, "side-by-side with the g(r)-based h(r)\n");
	fprintf(stderr, "[density] - density in units of molecules/A^3\n");
	fprintf(stderr, "[datfile] - two column file containing the g(r) and domain\n");
	exit(1);

}

int main(int argc, char **argv) {

	int i, n, c, fitting = 0;
	FILE *fp;
	char *datfile;
	double density, xp, A, r0, r1;
	struct gor_t *gor;
	extern char *optarg;
	extern int optind, optopt;

	if(argc < 7) usage((char *)argv[0]);

	/* get arguments */
	c = getopt(argc, (char **)argv, "f");
	if(c != -1) {
		switch (c) {
			case 'f':
				fitting = 1;
				break;
			case '?':
			usage((char *)argv[0]);
		}
	}

	/* get the density of particles for S(k) */
	density = atof(argv[optind]);

	/* get the parameters */
	xp = atof(argv[++optind]);
	A = atof(argv[++optind]);
	r0 = atof(argv[++optind]);
	r1 = atof(argv[++optind]);

	/* get the filename containing the g(r) */
	datfile = argv[++optind];
	if(!datfile) {
		fprintf(stderr, "error: invalid datafile\n");
		usage((char *)argv[0]);
	}
	else
		fp = fopen(datfile, "r");

	if(!fp) {
		fprintf(stderr, "error: could not open datafile %s\n", datfile);
		usage((char *)argv[0]);
	}


	/* load in the data */
	gor = (struct gor_t *)calloc(MAXLINE, sizeof(struct gor_t));
	if(!gor) {
		fprintf(stderr, "error: couldn't allocate working memory buffer\n");
		exit(1);
	}

	for(i = 0; n != EOF; i++) {
		n = fscanf(fp, "%lg %lg", &(gor + i)->domain, &(gor + i)->range);
	}
	fclose(fp);
	--i;

	if(!i) {
		fprintf(stderr, "error: datfile %s is empty\n", datfile);
		exit(1);
	}

	/* realloc to save memory */
	gor = realloc(gor, sizeof(struct gor_t)*i);
	if(!gor) {
		fprintf(stderr, "error: couldn't realloc working memory buffer\n");
		exit(1);
	}

	if(fitting)
		fit_gor(gor, i, xp, A, r0, r1);				/* get the fitted h(r) */
	else
		structure_factor(gor, i, density, xp, A, r0, r1);	/* get the structure factor */

	/* free our g(r) structure and exit */
	free(gor);
	exit(0);

}

