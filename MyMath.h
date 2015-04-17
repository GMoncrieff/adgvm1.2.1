#ifndef MyMath_h___
#define MyMath_h___


#include <math.h>
#include <stdlib.h>

#define MODMULT(a, b, c, m, s) q = s/a; s = b*(s-a*q)-c*q; if (s < 0) s += m; 


static long s1 = 1; 
static long s2 = 1; 


// Minimum of two values
double MyMin(double x, double y);

// Maximum of two values
double MyMax(double x, double y);

// Transform degree to rad
double DegToRad(double lat);

// absolute value
double MyAbs( double x );

// get random number
double GetRand(void);

// initialize random number generator
void initLCG(long init_s1, long init_s2); 

// get exponentially distributed random number
double rExp( double mean );

// get normal distributed random number
double rNormal( double x, double mu, double sigma );

// round a number to integer
int MyRound( double x );

// get ratio from humidity in %
double GetRatio( double humidity );


double SlidingAverage( double average, double new_val, double old_val, int num );


double MySignum( double x );

double MyFak( double val );

// -----------------------------------------------------------------------------


double MyMin(double x, double y)
{
	return x<y ? x : y;
}

double MyMax(double x, double y)
{
	return x>y ? x : y;
}


// -----------------------------------------------------------------------------


double DegToRad(double lat)
{
	return lat*M_PI/180.0;
}


// -----------------------------------------------------------------------------


double MyAbs( double x )
{
	if (x<0)
		return -x;
	else
		return x;
}


// -----------------------------------------------------------------------------


double GetRand(void)
{
	long q, z; 
	 
	MODMULT(53668, 40014, 12211, 2147483563L, s1) 
	MODMULT(52774, 40692, 3791, 2147483399L, s2) 
	z = s1 - s2; 
	if (z < 1) 
		z += 2147483562; 
		
	return z * 4.656613e-10; 
} 


// -----------------------------------------------------------------------------


void initLCG(long init_s1, long init_s2) 
{ 
	s1 = init_s1; 
	s2 = init_s2; 
}


// -----------------------------------------------------------------------------


double rExp( double mean )
{
	double r;
	
	r = drand48();
	
	return -log(r)*mean;
}


// -----------------------------------------------------------------------------


double rNormal( double x, double mu, double sigma )
{
	return 1./(sqrt(2.*M_PI)*sigma)*exp(-pow(x-mu,2)/(2.*pow(sigma,2)));
}


// -----------------------------------------------------------------------------


int MyRound( double x )
{
    if ( x-floor(x) >= 0.5 )
        return (int) ceil(x);
    else
        return (int) floor(x);
}


// -----------------------------------------------------------------------------


double GetRatio( double humidity )
{
	return humidity/100.0;
}


// -----------------------------------------------------------------------------


double SlidingAverage( double average, double new_val, double old_val, int num )
{
	return ( (double)num*average+new_val-old_val )/(double)num;
}



double MySignum( double x ) { return x>0 ? 1. : 0.; }


double MyFak( double val )
{
    if (val > 1)
      return MyFak(val-1) * val;
    
    return 1;
}


#endif



















