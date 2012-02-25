/************************************************/
/*						*/
/* Calculate the virial coefficients using "a"	*/
/* coefficients from the CRC			*/
/*						*/
/* Jonathan Belof				*/
/* Department of Chemistry			*/
/* University of South Florida			*/
/*						*/
/* compile with:				*/
/*	gcc -o v ./B_iso.c -lm			*/
/*						*/
/************************************************/


#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#define T0				298.15			/* standard temperature (K) */
#define DT				0.1			/* small interval in K */

void usage(char *progname) {

	fprintf(stderr, "usage: %s <start temp> <stop temp> <points> <a1> [a2] [a3] ...\n", progname);
	fprintf(stderr, "\t<start temp>: starting temperature (K)\n");
	fprintf(stderr, "\t<stop temp>: stopping temperature (K)\n");
	fprintf(stderr, "\t<points>: number of data points to calculate (integer)\n");
	fprintf(stderr, "\t<ai>: the i-th a-coefficient from the CRC\n");
	fprintf(stderr, "\nuseful parameters:\n");
	fprintf(stderr, "\tH2:\ta1=15.4 a2=-9.0 a3=-0.2\n");
	fprintf(stderr, "\tHe:\ta1=12 a2=-1\n");
	fprintf(stderr, "\tH2O:\ta1=-1158 a2=-5157 a3=-10301 a4=-10597 a5=-4415\n");
	fprintf(stderr, "\tUF6:\ta1=-1204 a2=-2690 a3=-2144\n");
	fprintf(stderr, "\tN2:\ta1=-4 a2=-56 a3=-12\n");
	fprintf(stderr, "\tO2:\ta1=-16 a2=-62 a3=-8 a4=-3\n");
	fprintf(stderr, "\tCO:\ta1=-9 a2=-58 a3=-18\n");
	fprintf(stderr, "\tNO:\ta1=-12 a2=-119 a3=89 a4=-73\n");
	fprintf(stderr, "\tCO2:\ta1=-127 a2=-288 a3=-118\n");
	fprintf(stderr, "\tAr:\ta1=-16 a2=-60 a3=-10\n");
	fprintf(stderr, "\tKr:\ta1=-51 a2=-118 a3=-29 a4=-5\n");
	fprintf(stderr, "\tXe:\ta1=-130 a2=-262 a3=-87\n");
	
	exit(1);

}

int main(int argc, char **argv) {

	double dT = DT;
	double temperature, temperature_lb, temperature_ub;	/* lower and upper temperature bounds (K) */
	int points;						/* number of data points to calculate */
	double B;						/* 2nd virial coefficient (cm^3/mol) */
	int i;							/* used for arg parsing */
	double current_a;					/* used in parsing out the "a" coefficients */

	if(argc < 5)
		usage(argv[0]);
	else { ++argv; --argc; };

	/* read in the main arguments */
	temperature_lb = atof(*argv); ++argv; --argc;
	temperature_ub = atof(*argv); ++argv; --argc;
	points = atoi(*argv); ++argv; --argc;

	if((temperature_lb > temperature_ub) || (temperature_lb < 0) || (temperature_ub < 0)) {
		fprintf(stderr, "%s: invalid temperature range specified\n", argv[0]);
		usage(argv[0]);
	}

	if(points < 1) {
		fprintf(stderr, "%s: invalid number of data points\n", argv[0]);
		usage(argv[0]);
	}

	for(temperature = temperature_lb; temperature <= temperature_ub; temperature += dT) {

		for(i = 0, B = 0; i < argc; i++) {
			current_a = atof(*(argv + i));
			B += current_a*(pow((T0/temperature - 1), i));
		}
		printf("%f %f\n", temperature, B);
	}


	exit(0);

}

