/************************************************/
/*						*/
/* Calculate the virial isotherm using the "a"	*/
/* coefficients from the CRC			*/
/*						*/
/* Jonathan Belof				*/
/* Department of Chemistry			*/
/* University of South Florida			*/
/*						*/
/* compile with:				*/
/*	gcc -std=c99 -o v ./virial_iso.c -lm	*/
/*						*/
/************************************************/


#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#define T0				298.15			/* standard temperature (K) */
#define NA				6.02214e23		/* Avogadro's number */
#define KB				1.3806503e-23		/* Boltzmann's constant */

#define ATM2PSC				101325			/* convert from atm to Pascals */
#define DP				100.0			/* small interval in Pascals */

void usage(char *progname) {

	fprintf(stderr, "usage: %s <mw> <temp> <start pressure> <stop pressure> <points> <a1> [a2] [a3] ...\n", progname);
	fprintf(stderr, "\t<mw>: molecular/atomic weight (g/mol)\n");
	fprintf(stderr, "\t<temp>: temperature of isotherm (K)\n");
	fprintf(stderr, "\t<start pressure>: starting pressure of isotherm (atm)\n");
	fprintf(stderr, "\t<stop pressure>: stopping pressure of isotherm (atm)\n");
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

	double mw;						/* molecular weight (g/mol) */
	double temperature;					/* temperature to perform the isotherm at (K) */
	double kT;						/* thermal energy (J) */
	double P;						/* current pressure (Pascals) */
	double dP;						/* pressure interval (Pascals) */
	double pressure_lb, pressure_ub;			/* lower and upper pressure bounds (atm) */
	int points;						/* number of data points to calculate */
	double B;						/* 2nd virial coefficient (cm^3/mol) */
	double density;						/* density, to be solved for (g/cm^3) */
	int i;							/* used for arg parsing */
	double current_a;					/* used in parsing out the "a" coefficients */

	if(argc < 7)
		usage(argv[0]);
	else { ++argv; --argc; };

	/* read in the main arguments */
	mw = atof(*argv); ++argv; --argc;
	temperature = atof(*argv); ++argv; --argc;
	pressure_lb = atof(*argv); ++argv; --argc;
	pressure_ub = atof(*argv); ++argv; --argc;
	points = atoi(*argv); ++argv; --argc;

	/* argument checking */
	if((mw < 1) || (mw > 267)) {	/* synthesize a new isotope lately? */
		fprintf(stderr, "%s: invalid molecular weight\n", argv[0]);
		usage(argv[0]);
	}

	if(temperature < 0) {
		fprintf(stderr, "%s: invalid temperature\n", argv[0]);
		usage(argv[0]);
	}

	if((pressure_lb > pressure_ub) || (pressure_lb < 0) || (pressure_ub < 0)) {
		fprintf(stderr, "%s: invalid pressure range specified\n", argv[0]);
		usage(argv[0]);
	}

	if(points < 1) {
		fprintf(stderr, "%s: invalid number of data points\n", argv[0]);
		usage(argv[0]);
	}

	/* read in the "a" coefficients */
	for(i = 0, B = 0; i < argc; i++) {
		current_a = atof(*(argv + i));
		B += current_a*(pow((T0/temperature - 1), i));
	}

	/* unit conversions to SI */
	pressure_lb *= ATM2PSC;				/* convert from atm to Pascals */
	pressure_ub *= ATM2PSC;				/* ditto */
	dP = (pressure_ub - pressure_lb)/points;	/* interval to use in calculating the function */
	B /= NA*1.0e6;					/* convert from cm^3/mol to m^3/molec*/
	kT = KB*temperature;				/* calculate kT (J) */

	for(P = pressure_lb; P <= pressure_ub; P += dP) {

		/* solve for the density quadratically */
		density = (-1 + sqrt(1 + 4*B*P/kT))/(2*B);
		density *= mw/(NA*1.0e6);	/* convert to g/cm^3 */

		if(!finite(density))		/* if we're in non-physical territory, don't bother with NaN's */
			break;
		else
			printf("%f %f\n", P/ATM2PSC, density);	/* output the pressure in atm, density in g/cm^3 */

	}

	exit(0);

}

