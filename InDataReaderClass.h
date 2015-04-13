#ifndef InDataReaderClass_h___
#define InDataReaderClass_h___

#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>

#include "Globals.h"
#include "GlobalTypes.h"
#include "InDataClass.h"

using namespace std;

// This file contains routines to read input data from global databases or
// from the shortlists.


class clInDataReader
{
	public:
		clInDataReader()  {};
		~clInDataReader() {};
		clInData getInData( double longitude, double latitude );
	
	private:
		double latitude_;
		double longitude_;
		clInData MyInData_;
		
		void ReadShortList();
		void ReadFiles();
		double ReadSOIL(const char *filename);
};




clInData clInDataReader::getInData( double longitude,  double latitude )
{
	latitude_  = latitude;
	longitude_ = longitude;
	string short_list_name;
	short_list_name = IN_DATA_HOME + "shortlists/shortlist_" + doubleToString(longitude_) + "_"
					+ doubleToString(latitude_) + "_0.dat";
	
	string shortlists_file_name = IN_DATA_HOME + "shortlists.slt";
	ifstream short_lists_file(shortlists_file_name.c_str());
	
	string   read_line;
	int      read_flag = 0;
	
	while ( short_lists_file >> read_line )
		if ( read_line == short_list_name ) read_flag = 1;
		
	short_lists_file.close();
	
	if ( read_flag == 0 )
	{
		cerr << "INI I must read in " << short_list_name << ", this may take some time... " << endl;
		
		// Read data from the files and assign them to variables.
		// Create file with importand information
		ReadFiles();
		
		short_lists_file.close();
		
		ofstream short_lists_file(shortlists_file_name.c_str(), ios::app);
		short_lists_file << short_list_name << endl;
		short_lists_file.close();
		
		cerr << "INI     ... data successfully read in. Now the simulation starts." << endl << endl;
	}
	
	ReadShortList();
	
	
	
	for ( int i=0; i<12; i++)
	{
		MyInData_.rdo_[i] = MyInData_.pwet_[i]*30.42;
		MyInData_.tmp_day_[i] = MyInData_.tmp_[i] + (MyInData_.tmp_max_[i]-MyInData_.tmp_[i])/2.;
	}
	
	// Model for soil layers not jet integrated, each layer has the same values
// 	for ( int i=0; i<MyInData_.soil_layers_; i++ ) MyInData_.theta_wp_[i]  = MyInData_.theta_wp_[0];
// 	for ( int i=0; i<MyInData_.soil_layers_; i++ ) MyInData_.theta_fc_[i]  = MyInData_.theta_fc_[0];
// 	for ( int i=0; i<MyInData_.soil_layers_; i++ ) MyInData_.bulk_dens_[i] = MyInData_.bulk_dens_[0];
	for ( int i=0; i<MyInData_.soil_layers_; i++ )
		MyInData_.theta_[i] = MyInData_.theta_wp_[i]+0.25*(MyInData_.theta_fc_[i]-MyInData_.theta_wp_[i]);
	
	
	MyInData_.atm_press_ = 101.325*pow((293.0-0.0065*MyInData_.elev_)/293.0,5.26)*1000.;
	MyInData_.latitude_  = latitude_;
	MyInData_.longitude_ = longitude_;
	
	return MyInData_;
}




void clInDataReader::ReadShortList()
{
	string file_name = IN_DATA_HOME + "shortlists/shortlist_" + doubleToString(longitude_) + "_" 
						+ doubleToString(latitude_) + "_0.dat";
	
	#ifdef S_GCM_180
	DEBUG( DEBUG_INPUT_DATA, "INI  USE PALEO CLIMATE FROM GCM, 180 PPM")
	file_name = IN_DATA_HOME + "shortlists/shortlist_" + doubleToString(longitude_) + "_" 
						+ doubleToString(latitude_) + "_180.dat";
	#endif
	
	#ifdef S_GCM_280
	DEBUG( DEBUG_INPUT_DATA, "INI  USE PALEO CLIMATE FROM GCM, 280 PPM")
	file_name = IN_DATA_HOME + "shortlists/shortlist_" + doubleToString(longitude_) + "_" 
						+ doubleToString(latitude_) + "_280.dat";
	#endif
	
	#ifdef S_GCM_400
	DEBUG( DEBUG_INPUT_DATA, "INI  USE PALEO CLIMATE FROM GCM, 400 PPM")
	file_name = IN_DATA_HOME + "shortlists/shortlist_" + doubleToString(longitude_) + "_" 
						+ doubleToString(latitude_) + "_400.dat";
	#endif
	
	DEBUG( DEBUG_INPUT_DATA,
				   "INI  Read shortlist                      " << file_name << endl <<
				   "INI -------------------------------------------------------------------------------------------------------------------------------" )
    
	ifstream Quelle( file_name.c_str() );
	
	if(!Quelle)
	{
		cerr << "INI ERROR Can't open file: " << file_name << endl;
		cerr << "INI ERROR Maybe the filename is listed in shortlists.slt but the file does not exist." << endl;
		exit(1);
	}
	
	for ( int i=0; i<12; i++ ) Quelle >> MyInData_.tmp_min_[i];
	for ( int i=0; i<12; i++ ) Quelle >> MyInData_.tmp_max_[i];
	for ( int i=0; i<12; i++ ) Quelle >> MyInData_.tmp_[i];
	for ( int i=0; i<12; i++ ) Quelle >> MyInData_.ralpha_[i];
	for ( int i=0; i<12; i++ ) Quelle >> MyInData_.rbeta_[i];
	for ( int i=0; i<12; i++ ) Quelle >> MyInData_.reh_[i];
	for ( int i=0; i<12; i++ ) Quelle >> MyInData_.pwet_[i];
	for ( int i=0; i<12; i++ ) Quelle >> MyInData_.sun_[i];
	for ( int i=0; i<12; i++ ) Quelle >> MyInData_.wnd_[i];
	for ( int i=0; i<12; i++ ) Quelle >> MyInData_.frost_[i];
	for ( int i=0; i<12; i++ ) MyInData_.frost_[i] = MyInData_.frost_[i]/31.;
	
// #   ifdef ELEPHANTS
//     for ( int i=0; i<12; i++ ) MyInData_.tmp_day_[i] += 0.45;
//     for ( int i=0; i<12; i++ ) MyInData_.tmp_min_[i] += 0.45;
//     for ( int i=0; i<12; i++ ) MyInData_.tmp_max_[i] += 0.45;
//     for ( int i=0; i<12; i++ ) MyInData_.tmp_[i]     += 0.45;
// #   endif

// #   ifdef ELEPHANTS_FUT   // kruger
//     for ( int i=0; i<12; i++ ) MyInData_.tmp_day_[i] += 4.272614;
//     for ( int i=0; i<12; i++ ) MyInData_.tmp_min_[i] += 4.272614;
//     for ( int i=0; i<12; i++ ) MyInData_.tmp_max_[i] += 4.272614;
//     for ( int i=0; i<12; i++ ) MyInData_.tmp_[i]     += 4.272614;
// #   endif
// #   ifdef ELEPHANTS_FUT   // chobe
// 	for ( int i=0; i<12; i++ ) MyInData_.tmp_day_[i] += 6.039497;
// 	for ( int i=0; i<12; i++ ) MyInData_.tmp_min_[i] += 6.039497;
// 	for ( int i=0; i<12; i++ ) MyInData_.tmp_max_[i] += 6.039497;
// 	for ( int i=0; i<12; i++ ) MyInData_.tmp_[i]     += 6.039497;
// #   endif
#   ifdef ELEPHANTS_FUT   // tsavo
for ( int i=0; i<12; i++ ) MyInData_.tmp_day_[i] += 4.396683;
for ( int i=0; i<12; i++ ) MyInData_.tmp_min_[i] += 4.396683;
for ( int i=0; i<12; i++ ) MyInData_.tmp_max_[i] += 4.396683;
for ( int i=0; i<12; i++ ) MyInData_.tmp_[i]     += 4.396683;
#   endif

// #   ifndef CLIM_SCEN
//     for ( int i=0; i<12; i++ ) MyInData_.tmp_day_[i] += 0.47;
//     for ( int i=0; i<12; i++ ) MyInData_.tmp_min_[i] += 0.47;
//     for ( int i=0; i<12; i++ ) MyInData_.tmp_max_[i] += 0.47;
//     for ( int i=0; i<12; i++ ) MyInData_.tmp_[i]     += 0.47;
// #   endif

// #   ifdef PLUS_TEMP
//     for ( int i=0; i<12; i++ ) MyInData_.tmp_day_[i] += 3.803;
//     for ( int i=0; i<12; i++ ) MyInData_.tmp_min_[i] += 3.803;
//     for ( int i=0; i<12; i++ ) MyInData_.tmp_max_[i] += 3.803;
//     for ( int i=0; i<12; i++ ) MyInData_.tmp_[i]     += 3.803;
// #   endif
	
	Quelle >> MyInData_.elev_;
	Quelle >> MyInData_.soil_layers_;
	//  	cout << MyInData_.elev_ << "  " <<  MyInData_.soil_layers_ << endl;
	
	MyInData_.soil_N_    = new double[MyInData_.soil_layers_];
	MyInData_.soil_C_    = new double[MyInData_.soil_layers_];
	MyInData_.theta_     = new double[MyInData_.soil_layers_];
	MyInData_.theta_wp_  = new double[MyInData_.soil_layers_];
	MyInData_.theta_fc_  = new double[MyInData_.soil_layers_];
	MyInData_.bulk_dens_ = new double[MyInData_.soil_layers_];
	MyInData_.g_theta_   = new double[MyInData_.soil_layers_];
	MyInData_.thickness_ = new double[MyInData_.soil_layers_];
	MyInData_.depth_     = new double[MyInData_.soil_layers_];
	
	for ( int i=0; i<MyInData_.soil_layers_; i++ ) Quelle >> MyInData_.theta_wp_[i];
	for ( int i=0; i<MyInData_.soil_layers_; i++ ) Quelle >> MyInData_.theta_fc_[i];
	
	for ( int i=0; i<MyInData_.soil_layers_; i++ ) Quelle >> MyInData_.bulk_dens_[i];
	for ( int i=0; i<MyInData_.soil_layers_; i++ ) Quelle >> MyInData_.soil_N_[i];
	for ( int i=0; i<MyInData_.soil_layers_; i++ ) Quelle >> MyInData_.soil_C_[i];
// 	for ( int i=0; i<MyInData_.soil_layers_; i++ ) MyInData_.soil_N_[i] = 775.;  // mean value for Africa
// 	for ( int i=0; i<MyInData_.soil_layers_; i++ ) MyInData_.soil_C_[i] = 7221.; // mean value for Africa
	for ( int i=0; i<MyInData_.soil_layers_; i++ ) Quelle >> MyInData_.thickness_[i];
	
	for ( int i=0; i<MyInData_.soil_layers_; i++ )
	{
		if ( MyInData_.soil_N_[i] <= 0.000001 ) MyInData_.soil_N_[i] = 0.01;
	}
	
	for ( int i=0; i<MyInData_.soil_layers_; i++ )
	{
		if ( MyInData_.soil_C_[i] <= 0.000001 ) MyInData_.soil_C_[i] = 0.01;
	}
	
	for ( int i=0; i<MyInData_.soil_layers_; i++ )
	{
		if ( MyInData_.soil_C_[i] >  30000. ) MyInData_.soil_C_[i] = 30000.;
	}
	
	// calculate depth of different soil layers
	MyInData_.depth_[0] = MyInData_.thickness_[0];
	for ( int i=1; i<MyInData_.soil_layers_; i++ ) MyInData_.depth_[i] = MyInData_.depth_[i-1]+MyInData_.thickness_[i];
	
//  	for ( int i=0; i<MyInData_.soil_layers_; i++ )
//   		cout << setw(14) << MyInData_.theta_wp_[i] << setw(14) << MyInData_.depth_[i] << setw(14) << MyInData_.soil_N_[i] << setw(14) << MyInData_.thickness_[i] << endl;
	
	if ( MyInData_.soil_N_[0]<=-1 )
	{
		for ( int i=0; i<12; i++ ) MyInData_.tmp_min_[i] = -10;
		for ( int i=0; i<12; i++ ) MyInData_.tmp_max_[i] = -10;
		for ( int i=0; i<12; i++ ) MyInData_.tmp_[i]     = -10;
		for ( int i=0; i<12; i++ ) MyInData_.ralpha_[i]  = -10;
		for ( int i=0; i<12; i++ ) MyInData_.rbeta_[i]   = -10;
		for ( int i=0; i<12; i++ ) MyInData_.reh_[i]     = -10;
		for ( int i=0; i<12; i++ ) MyInData_.pwet_[i]    = -10;
		for ( int i=0; i<12; i++ ) MyInData_.sun_[i]     = -10;
		for ( int i=0; i<12; i++ ) MyInData_.wnd_[i]     = -10;
		for ( int i=0; i<12; i++ ) MyInData_.frost_[i]   = -10;
		
		MyInData_.atm_press_    = -10;
		MyInData_.elev_         = -10;
		
		for ( int i=0; i<MyInData_.soil_layers_; i++ ) MyInData_.theta_wp_[i]  = -10;
		for ( int i=0; i<MyInData_.soil_layers_; i++ ) MyInData_.theta_fc_[i]  = -10;
		for ( int i=0; i<MyInData_.soil_layers_; i++ ) MyInData_.bulk_dens_[i] = -10;
		for ( int i=0; i<MyInData_.soil_layers_; i++ ) MyInData_.soil_N_[i]    = -10;
		for ( int i=0; i<MyInData_.soil_layers_; i++ ) MyInData_.soil_C_[i]    = -10;
	}
	
	
// 	for ( int i=0; i<MyInData_.soil_layers_; i++ )
// 	{
// 		MyInData_.theta_wp_[i]  = 0.159;
// 		MyInData_.theta_fc_[i]  = 0.396;
// 		MyInData_.soil_N_[i]    = 919.;
// 		MyInData_.soil_C_[i]    = 8712.;
// 	}
	
	Quelle.close();
	
	return;
}



void clInDataReader::ReadFiles()
{
	string	file_name = IN_DATA_HOME + "filenames.txt";
	char	climate_file_name_tmp[190];
	string  climate_file_name;
	char	line[500];
	double	min;
	double	lat;
	double	lon;
	double	err;
	int		row = 0;
	double	temp;
	arry12  rm;
	arry12  rcv;
	arry12  dtr;
	
	ifstream Quelle(file_name.c_str());
	if(!Quelle)
	{
		cerr << "INI ERROR Can't open file: " << file_name << endl;
		exit(1);
	}
	
			
	// Read the 8 (j=0..7) climate data sets first, grid10min_....
	for ( int j=0; j<=8; j++ )
	{
		// read name of next climate data set from filenames.txt
		Quelle.getline(climate_file_name_tmp, 190);
		climate_file_name = IN_DATA_HOME + (string)climate_file_name_tmp;
		
		// Open this file
		ifstream ClimateFile(climate_file_name.c_str());
	
		if(!ClimateFile)
		{
			cerr << "INI ERROR Can't open file: " << climate_file_name << endl;
			exit(1);
		}
		
		min = 500.;
		
		// find site in climate data which is next to study site given by lat and long
		for ( int i=1; i<=566260; i++ )
		{
			// read lat and long and ignore the rest of the line
			ClimateFile >> lat >> lon;
			ClimateFile.getline(line,500);
			
			err = pow(MyAbs(lon-longitude_),2)+pow(MyAbs(lat-latitude_),2);
			if ( err < min )
			{
				min = err;
				row = i;
			}
		}
		
		
		// set cursor to the start of the file
		ClimateFile.seekg(0);
		
		// go to the in the climate data
		for ( int k=1; k<row; k++ )
			ClimateFile.getline(line, 500);
		
		// read lat and long
		ClimateFile >> lat >> lon;
		
		// Read data of the 12 month
		if      (j==0) { for ( int i=0; i<12; i++ ) { ClimateFile >> dtr[i]; } } 
		else if (j==1) { for ( int i=0; i<12; i++ ) { ClimateFile >> MyInData_.tmp_[i]; } }
		else if (j==2) { for ( int i=0; i<12; i++ ) { ClimateFile >> MyInData_.wnd_[i]; } }
		else if (j==3) { for ( int i=0; i<12; i++ ) { ClimateFile >> MyInData_.sun_[i]; } }
		else if (j==4) { for ( int i=0; i<12; i++ ) { ClimateFile >> MyInData_.rdo_[i]; } }
		else if (j==5) { for ( int i=0; i<12; i++ ) { ClimateFile >> MyInData_.reh_[i]; } }
		else if (j==6) {                              ClimateFile >> MyInData_.elev_; }
		else if (j==7) { for ( int i=0; i<12; i++ ) { ClimateFile >> rm[i];  } 
						 for ( int i=0; i<12; i++ ) { ClimateFile >> rcv[i]; } }
		else if (j==8) { for ( int i=0; i<12; i++ ) { ClimateFile >> MyInData_.frost_[i]; } }
		
		ClimateFile.close();
		
	}
	
	
	MyInData_.soil_layers_ = 12;       // NOTE now number of soil layers is 12, need to read in soil files, when we have them
	
	MyInData_.soil_N_    = new double[MyInData_.soil_layers_];
	MyInData_.soil_C_    = new double[MyInData_.soil_layers_];
	MyInData_.theta_     = new double[MyInData_.soil_layers_];
	MyInData_.theta_wp_  = new double[MyInData_.soil_layers_];
	MyInData_.theta_fc_  = new double[MyInData_.soil_layers_];
	MyInData_.bulk_dens_ = new double[MyInData_.soil_layers_];
	MyInData_.g_theta_   = new double[MyInData_.soil_layers_];
	
	
	// All grid10**** files are read, read soil data now
	Quelle >> climate_file_name_tmp;
	climate_file_name = IN_DATA_HOME + (string)climate_file_name_tmp;
	temp = ReadSOIL( climate_file_name.c_str() );
	for ( int i=0; i<MyInData_.soil_layers_; i++) MyInData_.theta_wp_[i] = temp/1000.;
	
	Quelle >> climate_file_name_tmp;
	climate_file_name = IN_DATA_HOME + (string)climate_file_name_tmp;
	temp = ReadSOIL( climate_file_name.c_str());
	for ( int i=0; i<MyInData_.soil_layers_; i++) MyInData_.theta_fc_[i] = temp/1000.;
	
	Quelle >> climate_file_name_tmp;
	climate_file_name = IN_DATA_HOME + (string)climate_file_name_tmp;
	temp = ReadSOIL( climate_file_name.c_str() );
	for ( int i=0; i<MyInData_.soil_layers_; i++) MyInData_.bulk_dens_[i] = temp;
	
	Quelle >> climate_file_name_tmp;
	climate_file_name = IN_DATA_HOME + (string)climate_file_name_tmp;
	temp = ReadSOIL( climate_file_name.c_str() );
	for ( int i=0; i<MyInData_.soil_layers_; i++) MyInData_.soil_N_[i] = temp;
	
	Quelle >> climate_file_name_tmp;
	climate_file_name = IN_DATA_HOME + (string)climate_file_name_tmp;
	temp = ReadSOIL( climate_file_name.c_str() );
	for ( int i=0; i<MyInData_.soil_layers_; i++) MyInData_.soil_C_[i] = 1000.*temp;
	
	// Bring input data to the format which is required by the model
	for ( int i=0; i<12; i++ )
	{
		MyInData_.tmp_min_[i] = MyInData_.tmp_[i] - 0.5*dtr[i];
		MyInData_.tmp_max_[i] = MyInData_.tmp_[i] + 0.5*dtr[i];
		
		if ( rm[i] == 0 )
			MyInData_.ralpha_[i] = 1.;
		else
// 			MyInData_.ralpha_[i] = 1./pow(rcv[i]/100.,2)/1000;
			MyInData_.ralpha_[i] = 1./pow(rcv[i],2)/100.;
		
		if ( rm[i] == 0  ) 
			MyInData_.rbeta_[i] = 0;
		else
// 			MyInData_.rbeta_[i] = rm[i]/pow(rcv[i]/100.,2)/1000;
			MyInData_.rbeta_[i] = rm[i]/pow(rcv[i],2)/100.;
		
// 		cout << setw(14) << MyInData_.ralpha_[i] << setw(14) << MyInData_.rbeta_[i] << endl;
		
		MyInData_.pwet_[i] = MyInData_.rdo_[i]/30.42;
		MyInData_.sun_[i]  = MyInData_.sun_[i]/100.;
	}
	
	
	MyInData_.elev_ = MyInData_.elev_*1000.;
	
	// Write shortlist for the study site
	string output_file = IN_DATA_HOME + "shortlists/shortlist_" + doubleToString(longitude_) + "_"
						+ doubleToString(latitude_) + "_0.dat"; 
	
	ofstream Ausg(output_file.c_str(), ios::trunc);
	
	for ( int i=0; i<12; i++ )   Ausg << MyInData_.tmp_min_[i] << "   ";  // 1
	Ausg << endl;
	for ( int i=0; i<12; i++ )   Ausg << MyInData_.tmp_max_[i] << "   ";  // 2
	Ausg << endl;
	for ( int i=0; i<12; i++ )   Ausg << MyInData_.tmp_[i] << "   ";  // 3
	Ausg << endl;
	for ( int i=0; i<12; i++ )   Ausg << MyInData_.ralpha_[i] << "   ";  // 4
	Ausg << endl;
	for ( int i=0; i<12; i++ )   Ausg << MyInData_.rbeta_[i] << "   ";  // 5
	Ausg << endl;
	for ( int i=0; i<12; i++ )   Ausg << MyInData_.reh_[i] << "   ";  // 6
	Ausg << endl;
	for ( int i=0; i<12; i++ )   Ausg << MyInData_.pwet_[i] << "   ";  // 7
	Ausg << endl;
	for ( int i=0; i<12; i++ )   Ausg << MyInData_.sun_[i] << "   ";  // 8
	Ausg << endl;
	for ( int i=0; i<12; i++ )   Ausg << MyInData_.wnd_[i] << "   ";  // 9
	Ausg << endl;	
	for ( int i=0; i<12; i++ )   Ausg << MyInData_.frost_[i] << "   ";  // 10
	Ausg << endl;	
	
	Ausg << MyInData_.elev_ << endl;  // 11
	Ausg << MyInData_.soil_layers_ << endl;  // 12
	
	for ( int i=0; i<MyInData_.soil_layers_; i++ ) Ausg << MyInData_.theta_wp_[i] << "  ";  // 13
	Ausg << endl;
	for ( int i=0; i<MyInData_.soil_layers_; i++ ) Ausg << MyInData_.theta_fc_[i] << "  ";  // 14
	Ausg << endl;
	for ( int i=0; i<MyInData_.soil_layers_; i++ ) Ausg << MyInData_.bulk_dens_[i] << "  ";  // 15
	Ausg << endl;
	for ( int i=0; i<MyInData_.soil_layers_; i++ ) Ausg << MyInData_.soil_N_[i] << "  ";  // 16
	Ausg << endl;
	for ( int i=0; i<MyInData_.soil_layers_; i++ ) Ausg << MyInData_.soil_C_[i] << "  ";  // 17
	Ausg << endl;
	
// 	Ausg << "0.05 0.05 0.05 0.05 0.05 0.05 0.1 0.1 0.1  0.1  0.15 0.15 0.25 0.25 0.25 0.25" << endl;
	Ausg << "0.05 0.05 0.1  0.1  0.1  0.2  0.2 0.2 0.25 0.25 0.25 0.25" << endl;   // NOTE "standard" profile, need to get that from data
	
	Ausg.close();
	Quelle.close();

}


double clInDataReader::ReadSOIL(const char *filename)
{
	int		col;
	int		row;
	double	val;
	char	line[34560];
	
	ifstream Quelle(filename);
	
	if(!Quelle)
	{
		cerr << "INI ERROR Can't open file: " << filename << endl;
		exit(1);
	}
	
	row = (int) floor(-12.*(latitude_-84.)+1.);
	col = (int) floor(12.*(180.+longitude_)+1.);
	
	// add 6 to skip the first 6 lines of the file
	row += 6;
	
	// read the first (row-1) lines to move curser to desired line
	for (int i=1; i<row; i++) 
		Quelle.getline(line,34560);
	
	// read in value
	for (int i=1; i<=col; i++)
		Quelle >> val;
	
	
	Quelle.close();
	
	return val;
}




#endif














