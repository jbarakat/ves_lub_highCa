/* FUNCTION EVALUATION
 *  Evaluate right-hand side vector f in the system of ODEs dy/dt = f(y,t).
 *
 * REFERENCES
 *  Secomb et al, J Fluid Mech (1986)
 *  
 * PARAMETERS
 *  ny					[input]		number of ODEs
 *  t						[input]		abscissas
 *  y						[input]		solution
 *  f						[input]		function
 *  r   = y[0]  [input]		radius
 *  psi = y[1]  [input]		tilt angle
 *  cs  = y[2]  [input]		meridional curvature
 *  qs  = y[3]  [input]		transverse shear tension
 *  p   = y[4]  [input]		pressure
 *  sig = y[5]  [input]		mean tension
 *  A   = y[6]  [input]		total surface area
 *  V   = y[7]  [input]		total volume
 *  Q2  = y[8]  [input]		leakback flux
 *  S   = y[9]  [input]		total meridional arc length
 */

#ifndef FUNC_H
#define FUNC_H


/* HEADER FILES */
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>

using namespace std;

/* PROTOTYPES */

/* IMPLEMENTATIONS */
void func(int ny, double Ca,
          double t, double *y, double *f){
  int i;
  
  if (ny != 8)
    cout << "Error: ny should equal 8." << endl;
  
  // define variables
  double r   = y[0];
  double psi = y[1];
  double p   = y[2];
  double sig = y[3];
  double A   = y[4];
  double V   = y[5];
  double Q2  = y[6];
  double S   = y[7];
  //double F   = y[10];
	
	// check bounds on radius
	if (r > 0.999999999)
		r = 0.999999999;
	
	if (r <= 0)
		r = 1e-12;
	
	if (r > 0.002){ // check if far from end caps
  	double r2  = r*r;
  //	double a   = 2*M_PI*r;
  	double cos = gsl_sf_cos(psi);
  	double sin = gsl_sf_sin(psi);
  	double log = gsl_sf_log(r  );
  	double cosr = cos/r;
  	double sinr = sin/r;
  	double g, e;

  	// calculate local flow coefficient
  	g  = (8.0/(1.0 - r2));
		g *= (Q2 - 1.0 - (1.0 - r2)/(2.0*log));
  	g /= (1.0 + r2 + (1.0 - r2)/log);

  	// calculate local shear rate
  	// (expand e in Taylor series about r = 1)
  	double Dr  = 1.0 - r; 
  	double Dr2 = Dr*Dr; 
  	double Dr3 = Dr2*Dr;
  	double Dr4 = Dr2*Dr2;

  	e  = -(3.0/Dr2)*Q2;
  	e += (2.0/Dr )*(1.0 - Q2);
  	e += 1.0 - (29.0/20.0)*Q2;
  	e += Dr *(7.0/10.0 - (32.0/20.0)*Q2);
  	e += Dr2*(11.0/20.0 - (2749.0/2800.0)*Q2);
  	e += Dr3*(657.0/1400.0 - (309.0/350.0)*Q2);
  	e += Dr4*(237.0/560.0 - (45967.0/56000.0)*Q2);

  	// calculate function
  	f[0] = -sin;
  	f[1] = -p/sig - cosr;
  	f[2] = g*cos;
  	f[3] = -e;
  	f[4] = 2.0*M_PI*r;
  	f[5] = M_PI*r2*cos;
		f[6] = 0.0;
  	f[7] = 0.0;
	}
	else { // approximate p, sig as constant 
	       // and set dpsi/ds = -cos(psi)/r
  	double cos = gsl_sf_cos(psi);
  	double sin = gsl_sf_sin(psi);
		
		f[0] = -sin;
		f[1] = -p/(2.0*sig);
		f[2] = 0.0;
		f[3] = 0.0;
		f[4] = 2.0*M_PI*r;
		f[5] = M_PI*r*r*cos;
		f[6] = 0.0; 
		f[7] = 0.0;
	}

	for (i = 0; i < ny; i++){
		f[i] *= S;
	}
}



#endif
