/****************************************************************************************/
/* Determines the correlation time interval for a random sequence                       */
/*                                                                                      */
/* The input signal is specified in an input file where each line holds a single        */
/* floating point value.  The signal values are all shifted by the average value in     */
/* order that the signal may fluctuate about the x-axis.  The normalization constant is */
/* then calculated so that the autocorrelation function will be normalized.             */
/* A random signal will see a rapid decay time, whereas a signal with greater           */
/* periodicity would be evident by judging decay time until crossing the x-axis, thus   */
/* determining the correlation interval.  The default method is to integrate the        */
/* autocorrelation function to get the correlation time - an alternative method         */
/* based upon the slope of a semilog scale can also be calculated by the -l option      */
/* on the command line.  Use of both methods yields a reasonable value for the          */
/* correlation time of the data.                                                        */
/*                                                                                      */
/* This implementation utilizes the autocorrelation function as follows. The j-th value */
/* of the autocorrelation function (in latex) is given by:                              */
/*                                                                                      */
/*                      ac(j) = \sum_{i = j}^{N - i} s(j)*s(i+j)                        */
/*                                                                                      */
/* where s(k) is the value of the signal at time k and N is the total number of data    */
/* points in the signal.  The average of this sum is then taken by dividing by          */
/* the number of points processed and a normalization constant, where the               */
/* normalization constant is given by:                                                  */
/*                                                                                      */
/*                               \sum_i^N s(i)^2                                        */
/*                                                                                      */
/* usage:                                                                               */
/*      autocorr <-l> [datfile] > output_file                                           */
/*                                                                                      */
/* compilation:                                                                         */
/*      gcc -o autocorr autocorr.c -lm                                                  */
/*                                                                                      */
/* References:                                                                          */
/*      Ifeachor, Jervis "Digital Signal Processing: A Practical Approach"              */
/*      Mitra, "Digital Signal Processing"                                              */
/*      Frenkel, "Understanding Molecular Simulation"                                   */
/*                                                                                      */
/* @2005 Jonathan Belof                                                                 */
/* Space Research Group                                                                 */
/* Department of Chemistry                                                              */
/* University of South Florida                                                          */
/****************************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#define NUM_LINES		4000000

#define THRESHOLD		0.01		/* if the ac gets close to 0 within this amount, consider a "hit" */
#define DECAY_POINT		2		/* if this many "hits" occurs, consider the ac function decayed */

int autocorr(double *dat, int num, int log_scale) {

	int i, j;		/* working indices */
	double avg = 0;		/* pre-processed average for the signal */
	double norm = 0;	/* normalization constant */
	double ac = 0;		/* current autocorrelation value */
	int hits = 0;		/* number of times the ac has passed the x-axis */
	int decayed = 0;	/* used to flag when ac is decayed */
	double ac_integral;	/* integration value */
	int corr_time = 0;	/* correlation time estimated by integration */
	int log_corr_time = 0;	/* correlation time estimated by semilog extrapolation */

	/* calculate the uncorrelated average */
	for(i = 0; i < num; i++)
		avg += *(dat + i);
	avg /= num;

	/* translate the signal by the average to fluctuate about the x-axis */
	for(i = 0; i < num; i++)
		*(dat + i) -= avg;

	/* calculate the normalization constant */
	for(i = 0; i < num; i++) {
		norm += *(dat + i)*(*(dat + i));
	}
	norm /= num;

	/* autocorrelate the signal */
	for(i = 0; i < (num/2); i++) {

		ac = 0;
		for(j = 0; j < (num - i); j++) {
			ac +=  *(dat + j)*(*(dat + i + j));
		}
		/* normalize the autocorrelation */
		ac /= (num - i) * norm;

		/* determine if it is crossing the x-axis */
		if(fabs(ac) <= THRESHOLD) ++hits;

		/* determine if we have have sufficiently decayed */
		if(hits >= DECAY_POINT) decayed = 1;

		/* generate semilog scale, and determine correlation time from it's inverse slope */
		if(log_scale) {
			printf("%d %f\n", i, exp(ac));
			if(decayed && !log_corr_time) { /* take inverse slope until decay point */
				log_corr_time = (int)(ceil(((double)i)/(exp(1.0) - exp(ac))));
			}
		}
		else {
			printf("%d %f\n", i, ac);
			/* sum this contribution into the integral, until we get to the noisy part */
			if(!decayed) ac_integral += ac;
		}
	}

	/* return the estimated correlation time */
	if(log_scale)
		corr_time = log_corr_time;
	else
		corr_time = (int)ceil(ac_integral);

	return(corr_time);

}

void usage(char *progname) {

	fprintf(stderr, "usage: %s <-l> [datafile]\n", progname);
	fprintf(stderr, "options: -l output logarithmic scale\n");
	exit(1);

}

int main(int argc, char **argv) {

	int i, n = 0, log_scale = 0;
	char *datfile;
	FILE *fp;
	double *dat;
	int c;
	extern char *optarg;
	extern int optind, optopt;

	if(argc < 2) usage((char *)argv[0]);

	while ((c = getopt(argc, (char **)argv, "l")) != -1) {
		switch (c) {
			case 'l':
				log_scale = 1;
				break;
			case '?':
			usage((char *)argv[0]);
		}
	}

	datfile = argv[optind];
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

	/* load in the data  - XXX do without num_lines */
	dat = (double *)calloc(NUM_LINES, sizeof(double));
	if(!dat) {
		fprintf(stderr, "couldn't allocate working memory buffer\n");
		exit(1);
	}

	for(i = 0; n != EOF; i++) {
		n = fscanf(fp, "%lg", (dat + i));
	}
	fclose(fp);

	printf("# correlation time = %d\n", autocorr(dat, (i - 1), log_scale));

	free(dat);
	exit(0);

}

