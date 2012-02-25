/* 

Calculate the quantum second virial coefficient to second-order via the Wigner-Kirkwood
semiclassical expansion, including the third-order ideal quantum exchange term


@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define H			6.626068e-34		/* Planck's constant in J s */
#define HBAR			1.054571e-34		/* above divided by 2pi in J s */
#define KB			1.3806503e-23		/* Boltzmann's constant in J/K */
#define NA			6.0221415e23		/* Avogadro's number */
#define H2_MASS			3.348e-27		/* mass of H2 molecule in kg */

#define MAXLINE 256

void usage(char *progname) {

	fprintf(stderr, "usage: %s <min temp> <max temp> <inc temp> <pe file>\n", progname);
	exit(1);

}

int main(int argc, char **argv) {

	int i, N;
	FILE *fp_fit;
	char linebuf[MAXLINE];
	double r, r_min, r_max, r_inc;
	double first_derivative, second_derivative, integrand;
	double *r_input, *fit_input;
	char fit_file[MAXLINE];
	double temperature, temperature_min, temperature_max, temperature_inc;
	double B, B_classical, B_quantum_first, B_quantum_second, B_quantum_ideal;

	if(argc != 5)
		usage(argv[0]);

	temperature_min = atof(argv[1]);
	temperature_max = atof(argv[2]);
	temperature_inc = atof(argv[3]);
	strcpy(fit_file, argv[4]);

	/* read in the number of lines */
	fp_fit = fopen(fit_file, "r");
	for(N = 0; fgets(linebuf, MAXLINE, fp_fit); N++);
	fclose(fp_fit);

	/* allocate space */
	r_input = calloc(N, sizeof(double));
	fit_input = calloc(N, sizeof(double));

	/* read in the isotropic function */
	fp_fit = fopen(fit_file, "r");
	for(i = 0; i < N; i++)
		fscanf(fp_fit, "%lg %lg\n", &r_input[i], &fit_input[i]);
	fclose(fp_fit);

	/* determine independent domain */
	r_min = r_input[0];
	r_max = r_input[N-1];
	r_inc = r_input[1] - r_input[0];

	for(temperature = temperature_min; temperature <= temperature_max; temperature += temperature_inc) {

		/* calculate the classical part */
		B_classical = 0;
		for(r = r_min, i = 0; r < r_max; r += r_inc, i++) {
			integrand = (exp(-fit_input[i]/temperature) - 1.0)*pow(r, 2.0)*r_inc;

			B_classical += integrand;
		}
		B_classical *= -2.0*M_PI;
		B_classical *= 1.0e-24*NA;	/* convert to cm^3/mol */

		/* calculate the first order quantum correction */
		B_quantum_first = 0;
		for(r = r_min, i = 0; r < r_max; r += r_inc, i++) {

			/* take the central derivative */
			first_derivative = KB*(fit_input[i+1] - fit_input[i])/r_inc;
			integrand = exp(-fit_input[i]/temperature)*pow(first_derivative*r, 2.0)*r_inc;

			B_quantum_first += integrand;

		}
		B_quantum_first *= (H*H/H2_MASS)/(24.0*M_PI*pow(KB*temperature, 3.0));
		B_quantum_first *= 1.0e-4*NA;	/* convert to cm^3/mol */


		/* calculate the second order quantum correction */
		B_quantum_second = 0;
		for(r = r_min, i = 0; r < (r_max - r_inc); r += r_inc, i++) {

			first_derivative = KB*(fit_input[i+1] - fit_input[i])/r_inc;
			second_derivative =  KB*((fit_input[i+2] - fit_input[i+1])/r_inc - first_derivative)/r_inc;

			integrand = pow(second_derivative, 2.0) + 2.0*pow((first_derivative/r), 2.0);
			integrand += (10.0/(9.0*KB*temperature))*pow(first_derivative, 3.0)/r;
			integrand -= (5.0/(36.0*pow(KB*temperature, 2.0)))*pow(first_derivative, 4.0)/r;
			integrand *= exp(-fit_input[i]/temperature)*pow(r, 2.0)*r_inc;

			B_quantum_second += integrand;
		}
		B_quantum_second *= pow(H*H/H2_MASS, 2.0)/(960.0*M_PI*M_PI*M_PI*pow(KB*temperature, 4.0));
		B_quantum_second *= -1.0e11*NA;

		/* the exchange term for an ideal Bose-Einstein gas */
		B_quantum_ideal = -NA*pow((H*H/(2.0*M_PI*H2_MASS*KB*temperature)), (3.0/2.0))*pow(2.0, (-5.0/2.0));
		B_quantum_ideal *= 1.0e6;

		/* the total virial coefficient */
		B = B_classical + B_quantum_first + B_quantum_second + B_quantum_ideal;
		printf("%f %f\n", temperature, B);
		/*printf("%f %f %f %f %f %f\n", temperature, B, B_classical, B_quantum_first, B_quantum_second, B_quantum_ideal);*/

	}


	return(0);

}


