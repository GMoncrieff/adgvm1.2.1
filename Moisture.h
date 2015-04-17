#ifndef Moisture_h___
#define Moisture_h___

#include "rGamma.h"

using namespace std;


// This file contains functions to calculate rainfall sequneces
// and to update the soil bucket model.


// for each soil layer work out Gtheta; see Arora et al. 2003 for equations
// Output : GetST
void GetSoilTheta( double *diff_fc_wp, double *Theta, double *ThetaWP, double *GetST, int soil_layers );

// Simulate rain sequendce for year
// Input  : days, rdo, pwet, ralpha, rbeta
// Output : RYear
void RainFallYear( int days, arry12 rdo, arry12 pwet, arry12 ralpha, arry12 rbeta, arryYear RYear );

// Read values of daily rainfall from file
// void RainFallYearFromFile( char *filename, arryYear RYear );


// Update water content in soil layers
// Input  : drain, theta, ThetaFC, depth
// Output : BIn
void BucketIn( double drain, double *Theta, double *ThetaFC, double *thickness, int soil_layers );


// Input  : Et, Theta, ThetaWP, depth
// Output : BOut
void BucketOut( double Et, double *Theta, double *ThetaWP, double *thickness, int soil_layers );


clInData ManipulateSeasonality( double factor, arry12 ralpha, arry12 rbeta, arry12 pwet, arry12 rdo, clInData IData );

// ----------------------------------------------------------------------------------------------------------


void GetSoilTheta( double *diff_fc_wp, double *Theta, double *ThetaWP, double *GetST, int soil_layers )
{
	int i;
	double beta;
	
	for ( i=0; i<soil_layers; i++ )
	{
		beta = MyMax(0.,MyMin(1.,(Theta[i]-ThetaWP[i])*diff_fc_wp[i]));
		GetST[i] = 2.*beta-beta*beta;
	}
	
	return;
}



// This procedure generates a vector RYear, which contains the rain/day.
void RainFallYear( int days, arry12 rdo, arry12 pwet, arry12 ralpha, arry12 rbeta, arryYear RYear )
{
	int month;
	int lastmonth;
	int d;
	double monthlyrain;
	double eventsize = 0;
	double rain_sum;
	
	lastmonth = -1;
	rain_sum = 0;
	// for all days in the year
	for ( d=1; d<=days; d++ )
	{
		// compute the month
		month = (int) floor( ((double) d)/30.42 );
		
		// If a new month begins, we must compute the total rain for it.
		// Rain is assumed to be gamma-distributed over a month
		// with the parameters ralpha and rbeta from climate files.
		if ( month != lastmonth )
		{
			monthlyrain = MyGamma(ralpha[month])*rbeta[month];
			
			lastmonth = month;
			// if the number of rainy days in the current month is 0, we
			// will have no rain in this month.
			if ( rdo[month]<=0 )
				eventsize = 0; 
			else
				// else we compute the average amount of rain on a
				// rainy day.
				eventsize = monthlyrain/rdo[month];
		}
		
		if ( drand48()  <= pwet[month] )
			RYear[d-1] = rExp(eventsize); 
		else 
			RYear[d-1] = 0;
		
		rain_sum += RYear[d-1];
	}
	
	return;
}



void BucketIn( double drain, double *Theta, double *ThetaFC, double *thickness, int soil_layers )
{
	int i;
	double soilin;			   	//rain input (m/day)
	double sumSD;
	double *ThetaTmp = new double[soil_layers];
	double *SD       = new double[soil_layers];				//this is how much water can take up by soil SoilDry (m)
	
	sumSD    = 0;
	soilin   = drain*1e-3;
	
	for ( i=0; i<soil_layers; i++ )
		ThetaTmp[i] = Theta[i];
	
	
	for ( i=0; i<soil_layers; i++ )
	{
		SD[i] = (ThetaFC[i]-ThetaTmp[i])*thickness[i];
		sumSD = sumSD+SD[i];			// how much water can be taken up by all soil layers
	}
	
	if ( soilin > sumSD )				// more rain than can be taken up
	{
		for ( i=0; i<soil_layers; i++ )
			ThetaTmp[i] = ThetaFC[i];	// set water contents to fc
			
			soilin = 0;
	}
	
	if ( soilin > 0 )
	{
		for ( i=0; i<soil_layers; i++ )
		{
			if ( soilin>SD[i] )
			{
				ThetaTmp[i] = ThetaFC[i];
				soilin -= SD[i];
				soilin  = MyMax( 0, soilin );
			}
			else
			{
				ThetaTmp[i] = ThetaTmp[i]+soilin/thickness[i];
				soilin = 0;
			}
		}
	}
	
	for ( i=0; i<soil_layers; i++ )
		Theta[i] = ThetaTmp[i];
	
	delete[] ThetaTmp;
	delete[] SD;
	
	return;
}


void BucketOut( double Et, double *Theta, double *ThetaWP, double *thickness, int soil_layers )
{
	int i;
	double soilout;
	double sumSL;
	double *SL         = new double[soil_layers];
	double *ThetaTmp   = new double[soil_layers];
	double *ThetaWPTmp = new double[soil_layers];
	
	for ( i=0; i<soil_layers; i++ ) ThetaWPTmp[i] = ThetaWP[i];
	
	
	for ( i=0; i<soil_layers; i++ )
		ThetaTmp[i] = Theta[i];
	
	
	sumSL   = 0;
	soilout = Et*1e-3; 							//convert mm/day into m/day
	
	for ( i=0; i<soil_layers; i++ )
	{
		SL[i] = (ThetaTmp[i]-ThetaWPTmp[i])*thickness[i]; 			//how much water can be lost from soil (m)
		sumSL += SL[i];
	}
	
	if ( soilout > sumSL )
	{
		for ( i=0; i<soil_layers; i++ )
			ThetaTmp[i] = ThetaWPTmp[i];
		
		soilout = 0;
	}
	
	if ( soilout>0 )
	{
		for ( i=0; i<soil_layers; i++ )
		{
			if ( soilout>SL[i] )
			{
				ThetaTmp[i] = ThetaWPTmp[i];
				soilout -= SL[i];
				soilout  = MyMax( 0, soilout );
			}
			else
			{
				ThetaTmp[i] = ThetaTmp[i]-soilout/thickness[i];
				soilout = 0;
			}
		}
	}
	
	for ( i=0; i<soil_layers; i++ )
		Theta[i] = ThetaTmp[i];
	
	delete[] SL;
	delete[] ThetaTmp;
	delete[] ThetaWPTmp;
	
	return;
}


clInData ManipulateSeasonality( double factor, arry12 ralpha, arry12 rbeta, arry12 pwet, arry12 rdo, clInData IData )
{
	
	// manipulate parameters for rainfall generator
	if ( factor == -5 )  // same precipitation for all month
	{
		double alpha_mean = 0.;
		double beta_mean  = 0.;
		double pwet_mean  = 0.;
		double rdo_mean   = 0.;
		
		for ( int i=0; i<12; i++ )
		{
			alpha_mean += IData.ralpha_[i];
			beta_mean  += IData.rbeta_[i];
			pwet_mean  += IData.pwet_[i];
			rdo_mean   += IData.rdo_[i];
		}
		alpha_mean /= 12.;
		beta_mean  /= 12.;
		pwet_mean  /= 12.;
		rdo_mean   /= 12.;
		
		
		double sumold = 0.;
		for ( int i=0; i<12; i++ ) sumold += IData.rbeta_[i]/IData.ralpha_[i];
		
		for ( int i=0; i<12; i++ )
		{
			IData.ralpha_[i] = alpha_mean;
			IData.rbeta_[i]  = beta_mean;
			IData.pwet_[i]   = pwet_mean;
			IData.rdo_[i]    = rdo_mean;
		}
		
		double sumnew = 0.;
		for ( int i=0; i<12; i++ ) sumnew += IData.rbeta_[i]/IData.ralpha_[i];
		
		for ( int i=0; i<12; i++ )
			IData.rbeta_[i] = IData.rbeta_[i]*sumold/sumnew;
		
	}
	else
	{
		double alpha_new[12];
		double beta_new[12];
		double pwet_new[12];
		
		int ialpha[12];
		int ibeta[12];
		int ipwet[12];
		
		for ( int i=0; i<12; i++ ) ialpha[i] = ibeta[i] = ipwet[i] = 0;
		
		int ialpha1;
		int ialpha2;
		int ibeta1;
		int ibeta2;
		int ipwet1;
		int ipwet2;
		
		double salpha;
		double sbeta;
		double spwet;
		
		for ( int i=0; i<12; i++  )
		{
			for ( int j=0; j<12; j++ )
			{
				if ( ralpha[j] < ralpha[i] ) ialpha[i]++;
				if ( rbeta[j]  < rbeta[i]  ) ibeta[i]++;
				if ( pwet[j]   < pwet[i]   ) ipwet[i]++;
			}
		}
		
		// The following manipulates the parameters
		for ( int i=0; i<6; i++ )
		{
			for ( int j=0; j<12; j++ )
			{
				if ( ialpha[j]==i )      ialpha1 = j;
				if ( ialpha[j]==(11-i) ) ialpha2 = j;
				if ( ibeta[j] ==i )      ibeta1  = j;
				if ( ibeta[j] ==(11-i) ) ibeta2  = j;
				if ( ipwet[j] ==i )      ipwet1  = j;
				if ( ipwet[j] ==(11-i) ) ipwet2  = j;
			}
			
			salpha =      factor*ralpha[ialpha2];
			alpha_new[ialpha1] = ralpha[ialpha1]+salpha;
			alpha_new[ialpha2] = ralpha[ialpha2]-salpha;
			
			sbeta =      factor*rbeta[ibeta2];
			beta_new[ibeta1]  = rbeta[ibeta1]+sbeta;
			beta_new[ibeta2]  = rbeta[ibeta2]-sbeta;
			
			spwet =      factor*pwet[ipwet2];
			pwet_new[ipwet1]  = pwet[ipwet1]+spwet;
			pwet_new[ipwet2]  = pwet[ipwet2]-spwet;
		}
		
		double sumold=0.;
		double sumnew=0.;
		
		for ( int i=0; i<12; i++ )
		{
			sumold += rbeta[i]/ralpha[i];
			sumnew += beta_new[i]/alpha_new[i];
		}
		// 		cout << sumold << endl;
		
		for ( int i=0; i<12; i++ )
			beta_new[i] = beta_new[i]*sumold/sumnew;
		
		
		for ( int i=0; i<12; i++ )
		{
			IData.ralpha_[i] = alpha_new[i];
			IData.rbeta_[i]  = beta_new[i];
			IData.pwet_[i]   = pwet_new[i];
			IData.rdo_[i]    = pwet_new[i]*30.42;
		}
		
	}
	
	
	return IData;
}


#endif















