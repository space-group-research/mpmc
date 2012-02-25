/* 

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

#include <mc.h>

/* support the first 5 spherical harmonics, they are hard-coded but necessary for gaussian quadrature integration */
/* (doing something slick like a func ptr array [l][m] would be just as ugly) */
double rotational_basis(int type, int l, int m, double theta, double phi) {

	int m_abs;
	double legendre, Ylm;

	/* calculate the prefactor*legendre part */
	m_abs = abs(m);

	if(l == 0) {

		legendre = (1.0/2.0)*sqrt(1.0/M_PI);

	} else if(l == 1) {

		if(m_abs == 0)
			legendre = (1.0/2.0)*sqrt(3.0/M_PI)*cos(theta);
		else if(m_abs == 1)
			legendre = -1.0*(1.0/2.0)*sqrt(3.0/(2.0*M_PI))*sin(theta);

	} else if(l == 2) {

		if(m_abs == 0)
			legendre = (1.0/4.0)*sqrt(5.0/M_PI)*(3.0*pow(cos(theta), 2.0) - 1.0);
		else if(m_abs == 1)
			legendre = -1.0*(1.0/2.0)*sqrt(15.0/(2.0*M_PI))*(sin(theta)*cos(theta));
		else if(m_abs == 2)
			legendre = (1.0/4.0)*sqrt(15.0/(2.0*M_PI))*pow(sin(theta), 2.0);

	} else if(l == 3) {

		if(m_abs == 0)
			legendre = (1.0/4.0)*sqrt(7.0/M_PI)*(5.0*pow(cos(theta), 3.0) -3.0*cos(theta));
		else if(m_abs == 1)
			legendre = -1.0*(1.0/8.0)*sqrt(21.0/M_PI)*sin(theta)*(5.0*pow(cos(theta), 2.0) - 1.0);
		else if(m_abs == 2)
			legendre = (1.0/4.0)*sqrt(105.0/(2.0*M_PI))*pow(sin(theta), 2.0)*cos(theta);
		else if(m_abs == 3)
			legendre = -1.0*(1.0/8.0)*sqrt(35.0/M_PI)*pow(sin(theta), 3.0);

	} else if(l == 4) {

		if(m_abs == 0)
			legendre = (3.0/16.0)*sqrt(1.0/M_PI)*(35.0*pow(cos(theta), 4.0) -30.0*pow(cos(theta), 2.0) +3.0);
		else if(m_abs == 1)
			legendre = -1.0*(3.0/8.0)*sqrt(5.0/M_PI)*sin(theta)*(7.0*pow(cos(theta), 3.0) - 3.0*cos(theta));
		else if(m_abs == 2)
			legendre = (3.0/8.0)*sqrt(5.0/(2.0*M_PI))*pow(sin(theta), 2.0)*(7.0*pow(cos(theta), 2.0) - 1.0);
		else if(m_abs == 3)
			legendre = -1.0*(3.0/8.0)*sqrt(35.0/M_PI)*pow(sin(theta), 3.0)*cos(theta);
		else if(m_abs == 4)
			legendre = (3.0/16.0)*sqrt(35.0/(2.0*M_PI))*pow(sin(theta), 4.0);

	} else if(l == 5) {

		if(m_abs == 0)
			legendre = (1.0/16.0)*sqrt(11.0/M_PI)*(63.0*pow(cos(theta), 5.0) -70*pow(cos(theta), 3.0) +15*cos(theta));
		else if(m_abs == 1)
			legendre = -1.0*(1.0/16.0)*sqrt(165.0/(2.0*M_PI))*sin(theta)*(21.0*pow(cos(theta), 4.0) -14.0*pow(cos(theta), 2.0) + 1.0);
		else if(m_abs == 2)
			legendre = (1.0/8.0)*sqrt(1155.0/(2.0*M_PI))*pow(sin(theta), 2.0)*(3.0*pow(cos(theta), 3.0) - cos(theta));
		else if(m_abs == 3)
			legendre = -1.0*(1.0/32.0)*sqrt(385.0/M_PI)*pow(sin(theta), 3.0)*(9.0*pow(cos(theta), 2.0) - 1.0);
		else if(m_abs == 4)
			legendre = (3.0/16.0)*sqrt(385.0/(2.0*M_PI))*pow(sin(theta), 4.0)*cos(theta);
		else if(m_abs == 5)
			legendre = -1.0*(3.0/32.0)*sqrt(77.0/M_PI)*pow(sin(theta), 5.0);

	} else if(l == 6) {

		if(m_abs == 0)
			legendre = (1.0/32.0)*sqrt(13.0/M_PI)*(231.0*pow(cos(theta), 6.0) - 315.0*pow(cos(theta), 4.0) + 105.0*pow(cos(theta), 2.0) - 5.0);
		else if(m_abs == 1)
			legendre = -1.0*(1.0/16.0)*sqrt(273.0/(2.0*M_PI))*sin(theta)*(33.0*pow(cos(theta), 5.0) - 30.0*pow(cos(theta), 3.0) + 5.0*cos(theta));
		else if(m_abs == 2)
			legendre = (1.0/64.0)*sqrt(1365.0/M_PI)*pow(sin(theta), 2.0)*(33.0*pow(cos(theta), 4.0) - 18.0*pow(cos(theta), 2.0) + 1.0);
		else if(m_abs == 3)
			legendre = -1.0*(1.0/32.0)*sqrt(1365.0/M_PI)*pow(sin(theta), 3.0)*(11.0*pow(cos(theta), 3.0) - 3.0*cos(theta));
		else if(m_abs == 4)
			legendre = (3.0/32.0)*sqrt(91.0/(2.0*M_PI))*pow(sin(theta), 4.0)*(11.0*pow(cos(theta), 2.0) - 1.0);
		else if(m_abs == 5)
			legendre = -1.0*(3.0/32.0)*sqrt(1001.0/M_PI)*pow(sin(theta), 5.0)*cos(theta);
		else if(m_abs == 6)
			legendre = (1.0/64.0)*sqrt(3003.0/M_PI)*pow(sin(theta), 6.0);

	} else if(l == 7) {

		if(m_abs == 0)
			legendre = (1.0/32.0)*sqrt(15.0/M_PI)*(429.0*pow(cos(theta), 7.0) - 693.0*pow(cos(theta), 5.0) + 315.0*pow(cos(theta), 3.0) - 35.0*cos(theta));
		else if(m_abs == 1)
			legendre = -1.0*(1.0/64.0)*sqrt(105.0/(2.0*M_PI))*sin(theta)*(429.0*pow(cos(theta), 6.0) - 495.0*pow(cos(theta), 4.0) + 135.0*pow(cos(theta), 2.0) - 5.0);
		else if(m_abs == 2)
			legendre = (3.0/64.0)*sqrt(35.0/M_PI)*pow(sin(theta), 2.0)*(143.0*pow(cos(theta), 5.0) - 110.0*pow(cos(theta), 3.0) + 15.0*cos(theta));
		else if(m_abs == 3)
			legendre = -1.0*(3.0/64.0)*sqrt(35.0/(2.0*M_PI))*pow(sin(theta), 3.0)*(143.0*pow(cos(theta), 4.0) - 66.0*pow(cos(theta), 2.0) + 3.0);
		else if(m_abs == 4)
			legendre = (3.0/32.0)*sqrt(385.0/(2.0*M_PI))*pow(sin(theta), 4.0)*(13.0*pow(cos(theta), 3.0) - 3.0*cos(theta));
		else if(m_abs == 5)
			legendre = -1.0*(3.0/64.0)*sqrt(385.0/(2.0*M_PI))*pow(sin(theta), 5.0)*(13.0*pow(cos(theta), 2.0) - 1.0);
		else if(m_abs == 6)
			legendre = (3.0/64.0)*sqrt(5005.0/M_PI)*pow(sin(theta), 6.0)*cos(theta);
		else if(m_abs == 7)
			legendre = -1.0*(3.0/64.0)*sqrt(715.0/(2.0*M_PI))*pow(sin(theta), 7.0);

	} else if(l == 8) {

		if(m_abs == 0)
			legendre = (1.0/256.0)*sqrt(17.0/M_PI)*(6435.0*pow(cos(theta), 8.0) - 12012.0*pow(cos(theta), 6.0) + 6930.0*pow(cos(theta), 4.0) - 1260.0*pow(cos(theta), 2.0) + 35.0);
		else if(m_abs == 1)
			legendre = -1.0*(3.0/64.0)*sqrt(17.0/(2.0*M_PI))*sin(theta)*(715.0*pow(cos(theta), 7.0) - 1001.0*pow(cos(theta), 5.0) + 385.0*pow(cos(theta), 3.0) - 35.0*cos(theta));
		else if(m_abs == 2)
			legendre = (3.0/128.0)*sqrt(595.0/M_PI)*pow(sin(theta), 2.0)*(143.0*pow(cos(theta), 6.0) - 143.0*pow(cos(theta), 4.0) + 33.0*pow(cos(theta), 2.0) - 1.0);
		else if(m_abs == 3)
			legendre = -1.0*(1.0/64.0)*sqrt(19635.0/(2.0*M_PI))*pow(sin(theta), 3.0)*(39.0*pow(cos(theta), 5.0) - 26.0*pow(cos(theta), 3.0) + 3.0*cos(theta));
		else if(m_abs == 4)
			legendre = (3.0/128.0)*sqrt(1309.0/(2.0*M_PI))*pow(sin(theta), 4.0)*(65.0*pow(cos(theta), 4.0) - 26.0*pow(cos(theta), 2.0) + 1.0);
		else if(m_abs == 5)
			legendre = -1.0*(3.0/64.0)*sqrt(17017.0/(2.0*M_PI))*pow(sin(theta), 5.0)*(5.0*pow(cos(theta), 3.0) - cos(theta));
		else if(m_abs == 6)
			legendre = (1.0/128.0)*sqrt(7293.0/M_PI)*pow(sin(theta), 6.0)*(15.0*pow(cos(theta), 2.0) - 1.0);
		else if(m_abs == 7)
			legendre = -1.0*(3.0/64.0)*sqrt(12155.0/(2.0*M_PI))*pow(sin(theta), 7.0)*cos(theta);
		else if(m_abs == 8)
			legendre = (3.0/256.0)*sqrt(12155.0/(2.0*M_PI))*pow(sin(theta), 8.0);

	}


	/* tack on the e^{i*m*phi} */
	if(m < 0) {	/* use the fact that Y_l(-m) = (-1)^m*cc(Y_lm) */

		if(type == REAL)
			Ylm = legendre*cos(fabs((double)m)*phi);
		else if(type == IMAGINARY)
			Ylm = -1.0*legendre*sin(fabs((double)m)*phi);

		Ylm *= pow(-1.0, fabs((double)m));

	} else {

		if(type == REAL)
			Ylm = legendre*cos(((double)m)*phi);
		else if(type == IMAGINARY)
			Ylm = legendre*sin(((double)m)*phi);

	}



	return(Ylm);

}

