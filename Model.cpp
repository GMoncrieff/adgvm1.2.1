#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

// for optimization, the target parameters are defined as global variables
#ifdef OPTIM_GLOBALS
double CLD_TREE;
double K_CAN_EXT_TREE[2];
double K_CAN_EXT_TREE_INVERSE[2];
double BETA_STEM_TREE;
double BETA_ROOT_TREE;
double UPSILON_STEM_TREE;
double UPSILON_ROOT_TREE;
double SIGMA_GROW_RESP_TREE[2];
double SIGMA_GROW_RESP_GRASS[2];
double IGNITION_PROB;
double IGNITION_PAR_2;
double SEED_GERM_PROB[2];
double PROB_ROOT_SUCKER[2];
double DEATH_PROB_FROST[2];
double DEATH_PROB_CARBON[2];
double DEATH_PROB_COMP[2];
double LIGHT_COMP_SAVSAV;
double LIGHT_COMP_SAVC4;
double LIGHT_COMP_SAVC3;
double LIGHT_COMP_FORSAV;
double TOP_KILL_CONST[2];
double TOP_KILL_H[2];
double TOP_KILL_I[2];

double L_SAV_SAV;
double L_SAV_FOR;
double L_SAV_C4G;
double L_SAV_C3G;

double L_FOR_SAV;
double L_FOR_FOR;
double L_FOR_C4G;
double L_FOR_C3G;

double L_C4G_SAV;
double L_C4G_FOR;
double L_C4G_C4G;
double L_C4G_C3G;

double L_C3G_SAV;
double L_C3G_FOR;
double L_C3G_C4G;
double L_C3G_C3G;

double LICMP_1[8][8];
double LICMP_2[8][8];


#endif

int GLOB_YEAR;				// GLOBAL VARIABLE
double tmp_min_month = 0;	// GLOBAL VARIABLE!!!! Needed for Aindex
double phen_counter  = 0;	// GLOBAL VARIABLE!!!! Counts swiches between metabolic and dormant state of tree 0
double A0_max_C3_day = 0;	// GLOBAL VARIABLE!!!! Index of day with maximum photosynthesis, needed for reproduction



#ifdef S_MANIP_FIRE
double glob_fire_param = 0.25;   // parameter used to manipulate fire probability
#endif


#ifdef GRAZING
double grazing_rate=1.;
#endif

#pragma warning(disable:981)

#include "MyNrutil.h"
#include "BigModelDebug.h"
#include "MyMath.h"

#ifdef S_AUSTRALIA
  #include "GlobalsAustralia.h"
#elif defined TRANSPLANT
  #include "GlobalsTransplant.h"
#else
  #include "Globals.h"
#endif


#ifdef S_ELEPHANTS
#include "ElephantClass.h"
#endif

#include "GlobalTypes.h"
#include "PenmanMonteith.h"
#include "Radiation.h"
#include "Leaf.h"
#include "InDataClass.h"
#include "InDataReaderClass.h"
#include "Fxn.h"
#include "Moisture.h"
#include "Fire.h"
#include "GrassPopClass.h"
#include "TreePopClass.h"
#include "YearlyDataWriterClass.h"
#include "SoilClass.h"


int main( int argc, char **argv )
{
	
#if DEBUG_ON
	__active_debug_flags__ = DEBUG_INPUT_DATA | DEBUG_YEARLY_DATA | DEBUG_FIRE_DATA;
#endif
	
	
//	abort if number of arguments is wrong
	if ( argc < 9 )
	{
		cerr << endl << "INI ERROR: wrong number of arguments." << endl;
		cerr << "INI ERROR: Needed arguments: ./Model     trees years lon lat fire rseed scenario output_dir" << endl;
		cerr << "INI ERROR: Abort simulation..." << endl << endl;
		return 1;
	}
	
	// initialize variables for simulation
	int		count_years			= 0;	// counter for years
	int		month_in_year		= 0;	// counter for month
	int 	day_in_year			= 0;	// counter for days
	int		fire_num 			= 0;
	int		frost				= 0;
	int     fire_flag			= 0;
	
// 	double radnetsum = 0;
	
	
	double	sunhrs_day;					// sunhours per day
	double	radnet_day;					// net radiation per day
	double	mmsTOkgd;					// micromol to mol;mol to g; growth resp loss, carbon to dry mass; seconds in 12 hours
	double	EtSite				= 0;	// evapotranspiration, total
	double	EtSiteGrass			= 0;	// evapotranspiration, grass
	double	EtSiteTrees			= 0;	// evapotranspiration, trees
	double	EtSiteGround		= 0;	// evapotranspiration, ground
	
	double	EtSiteRef;					// reference evapotranspiration
	double	evapo_sum_year		= 0;
	double	rain_sum_year		= 0;
	double	rain_sum			= 0;
	double	evapo_sum			= 0;
	double	evapo_ref_sum		= 0.;
	double	evapo_grass_sum		= 0;
	double	evapo_soil_sum		= 0;
	
	double	total_fire_intensity= 0.;
	double	total_basal_area	= 0.;
	double	total_npp			= 0.;
	double	total_nee			= 0.;
	
	double	T_fac				= 0;  // needed for f(T) function of MaintResp
	
	double		*diff_fc_wp;				// 1/(theta_fc-theta_wp), needed for GetSoilTheta
	arryYear	Rain;					// array for rain distributon over the year
	arryYear	Ignitions;				// stores days where ignition is possible
	for ( int m=0; m<365; m++ ) Ignitions[m] = 0;
	
	clYearlyDataWriter myYearlyData(YEARLY_DATA_LENGTH);
	double yearly_output_data[YEARLY_DATA_LENGTH];
	
	
	// initialize variables with input arguments
	int			num_of_trees	 = atoi( argv[1] );
	int			years_to_run	 = atoi( argv[2] );
	int			with_fire		 = atoi(argv[5]);
	double		latitude 		 = (double) atof( argv[4] ); // y-coordinate
	double		longitude		 = (double) atof( argv[3] ); // x-coordinate
	double		d_seed			 = (double) atof( argv[6] );
	int			i_scen			 = (int) atoi( argv[7] );
	string		s_seed			 = (string) argv[6];
	string		s_scen			 = (string) argv[7];
	string		output_directory = OUT_DATA_HOME+(string) argv[8]+"/";

    double      d_params[INPUT_PARAMS];
    string      s_params[INPUT_PARAMS];
    for ( int i=0; i<INPUT_PARAMS; i++ ) d_params[i]=0;
    for ( int i=0; i<INPUT_PARAMS; i++ ) s_params[i]="0";


#   ifdef S_ELEPHANTS
    if ( argc < 20 )
    {
        cerr << endl << "INI WARNING: Wrong number of arguments. No elephant impact simulated." << endl;
    }
    else
    {
        for ( int i=0; i<11; i++ ) s_params[i] = (string) argv[9+i];
        for ( int i=0; i<11; i++ ) d_params[i] = atof(argv[9+i]);
    }

    clElephants myElephants( d_params[0], d_params[1], d_params[2], d_params[3],
                             d_params[4], d_params[5], d_params[6], d_params[7],
                             d_params[8], d_params[9], d_params[10] );
#   endif


	#ifdef S_CLIM_SENSITIVITY  // allows to set CO2, temperature and rainfall
	if ( argc != 13 )
	{
		cerr << endl << "INI WARNING: Wrong number of arguments. Standard climate used." << endl;
		// NOTE define "standard climate "
		s_params[0] = "387.";  // CO2
		d_params[0] =  387.;
		s_params[1] = "0.";    // temperature
		d_params[1] =  0.;
		s_params[2] = "1.";    // mean prec
		d_params[2] =  1.;
		s_params[3] = "0.";    // seasonality
		d_params[3] =  0.;
	}
	else
	{
		for ( int i=0; i<4; i++ ) s_params[i] = (string) argv[9+i];
		for ( int i=0; i<4; i++ ) d_params[i] = atof(argv[9+i]);
	}
	#endif  // S_CLIM_SENSITIVITY


#   ifdef OPTIM_GLOBALS
//     cout << "---- OPTIM_GLOBALS -----" << endl;
	if ( argc < 38 )
	{
		cerr << endl << "INI WARNING: wrong number of arguments for optimization. Standard values used." << endl;
		for ( int i=0; i<INPUT_PARAMS; i++ ) s_params[i] = "1";
		for ( int i=0; i<INPUT_PARAMS; i++ ) d_params[i] = 1.;
		
		CLD_TREE                 = 0.02;
		K_CAN_EXT_TREE[0]        = 0.5;
		K_CAN_EXT_TREE[1]        = 0.4;
		K_CAN_EXT_TREE_INVERSE[0]   = 1./K_CAN_EXT_TREE[0];
		K_CAN_EXT_TREE_INVERSE[1]   = 1./K_CAN_EXT_TREE[1];
		IGNITION_PROB            = 0.01;
		IGNITION_PAR_2           = 0.1;
		UPSILON_STEM_TREE        = 150.;
		UPSILON_ROOT_TREE        = UPSILON_STEM_TREE*6./15.;
		SIGMA_GROW_RESP_TREE[0]     = 0.35;
		SIGMA_GROW_RESP_TREE[1]     = 0.35;
		SIGMA_GROW_RESP_GRASS[0]    = 0.35;
		SIGMA_GROW_RESP_GRASS[1]    = 0.35;
		PROB_ROOT_SUCKER[0]         = 0.;
		PROB_ROOT_SUCKER[1]         = 0.;
		BETA_STEM_TREE           = 0.025;
		BETA_ROOT_TREE           = BETA_STEM_TREE;
		SEED_GERM_PROB[0]           = 0.25;
		SEED_GERM_PROB[1]           = 0.25;
		DEATH_PROB_FROST[0]         = 0.001;
		DEATH_PROB_FROST[1]         = 0.001;
		DEATH_PROB_CARBON[0]        = 0.001;
		DEATH_PROB_CARBON[1]        = 0.001;
		DEATH_PROB_COMP[0]          = 0.001;
		DEATH_PROB_COMP[1]          = 0.0005;
		LIGHT_COMP_SAVSAV           = 0.5;
		LIGHT_COMP_SAVC4            = 0.5;
		LIGHT_COMP_SAVC3            = 0.15;
		LIGHT_COMP_FORSAV           = 0.5;
		TOP_KILL_CONST[0]           = 4.3;
		TOP_KILL_CONST[1]           = 6.3;
		TOP_KILL_H[0]               = 5.003;
		TOP_KILL_H[1]               = 3.003;
		TOP_KILL_I[0]               = 0.004408;
		TOP_KILL_I[1]               = 0.006408;
		
		L_SAV_SAV = 0.5;
		L_SAV_FOR = 0.15;
		L_SAV_C4G = 0.5;
		L_SAV_C3G = 0.15;

		L_FOR_SAV = 0.5;
		L_FOR_FOR = 0.15;
		L_FOR_C4G = 0.5;
		L_FOR_C3G = 0.15;

		L_C4G_SAV = 0.5;
		L_C4G_FOR = 0.15;
		L_C4G_C4G = 0.5;
		L_C4G_C3G = 0.15;

		L_C3G_SAV = 0.5;
		L_C3G_FOR = 0.15;
		L_C3G_C4G = 0.5;
		L_C3G_C3G = 0.15;

		LICMP_1[0][0] = L_SAV_SAV;
		LICMP_1[1][0] = L_FOR_SAV;
		LICMP_1[2][0] = L_C4G_SAV;
		LICMP_1[3][0] = L_C4G_SAV;
		LICMP_1[4][0] = L_C4G_SAV;
		LICMP_1[5][0] = L_C3G_SAV;
		LICMP_1[6][0] = L_C3G_SAV;
		LICMP_1[7][0] = L_C3G_SAV;
		LICMP_1[0][1] = L_SAV_FOR;
		LICMP_1[1][1] = L_FOR_FOR;
		LICMP_1[2][1] = L_C4G_FOR;
		LICMP_1[3][1] = L_C4G_FOR;
		LICMP_1[4][1] = L_C4G_FOR;
		LICMP_1[5][1] = L_C3G_FOR;
		LICMP_1[6][1] = L_C3G_FOR;
		LICMP_1[7][1] = L_C3G_FOR;
		LICMP_1[0][2] = L_SAV_C4G;
		LICMP_1[1][2] = L_FOR_C4G;
		LICMP_1[2][2] = L_C4G_C4G;
		LICMP_1[3][2] = L_C4G_C4G;
		LICMP_1[4][2] = L_C4G_C4G;
		LICMP_1[5][2] = L_C3G_C4G;
		LICMP_1[6][2] = L_C3G_C4G;
		LICMP_1[7][2] = L_C3G_C4G;
		LICMP_1[0][3] = L_SAV_C4G;
		LICMP_1[1][3] = L_FOR_C4G;
		LICMP_1[2][3] = L_C4G_C4G;
		LICMP_1[3][3] = L_C4G_C4G;
		LICMP_1[4][3] = L_C4G_C4G;
		LICMP_1[5][3] = L_C3G_C4G;
		LICMP_1[6][3] = L_C3G_C4G;
		LICMP_1[7][3] = L_C3G_C4G;
		LICMP_1[0][4] = L_SAV_C4G;
		LICMP_1[1][4] = L_FOR_C4G;
		LICMP_1[2][4] = L_C4G_C4G;
		LICMP_1[3][4] = L_C4G_C4G;
		LICMP_1[4][4] = L_C4G_C4G;
		LICMP_1[5][4] = L_C3G_C4G;
		LICMP_1[6][4] = L_C3G_C4G;
		LICMP_1[7][4] = L_C3G_C4G;
		LICMP_1[0][5] = L_SAV_C3G;
		LICMP_1[1][5] = L_FOR_C3G;
		LICMP_1[2][5] = L_C4G_C3G;
		LICMP_1[3][5] = L_C4G_C3G;
		LICMP_1[4][5] = L_C4G_C3G;
		LICMP_1[5][5] = L_C3G_C3G;
		LICMP_1[6][5] = L_C3G_C3G;
		LICMP_1[7][5] = L_C3G_C3G;
		LICMP_1[0][6] = L_SAV_C3G;
		LICMP_1[1][6] = L_FOR_C3G;
		LICMP_1[2][6] = L_C4G_C3G;
		LICMP_1[3][6] = L_C4G_C3G;
		LICMP_1[4][6] = L_C4G_C3G;
		LICMP_1[5][6] = L_C3G_C3G;
		LICMP_1[6][6] = L_C3G_C3G;
		LICMP_1[7][6] = L_C3G_C3G;
		LICMP_1[0][7] = L_SAV_C3G;
		LICMP_1[1][7] = L_FOR_C3G;
		LICMP_1[2][7] = L_C4G_C3G;
		LICMP_1[3][7] = L_C4G_C3G;
		LICMP_1[4][7] = L_C4G_C3G;
		LICMP_1[5][7] = L_C3G_C3G;
		LICMP_1[6][7] = L_C3G_C3G;
		LICMP_1[7][7] = L_C3G_C3G;
		
		LICMP_2[0][0] = 1.- L_SAV_SAV;
		LICMP_2[1][0] = 1.- L_FOR_SAV;
		LICMP_2[2][0] = 1.- L_C4G_SAV;
		LICMP_2[3][0] = 1.- L_C4G_SAV;
		LICMP_2[4][0] = 1.- L_C4G_SAV;
		LICMP_2[5][0] = 1.- L_C3G_SAV;
		LICMP_2[6][0] = 1.- L_C3G_SAV;
		LICMP_2[7][0] = 1.- L_C3G_SAV;
		LICMP_2[0][1] = 1.- L_SAV_FOR;
		LICMP_2[1][1] = 1.- L_FOR_FOR;
		LICMP_2[2][1] = 1.- L_C4G_FOR;
		LICMP_2[3][1] = 1.- L_C4G_FOR;
		LICMP_2[4][1] = 1.- L_C4G_FOR;
		LICMP_2[5][1] = 1.- L_C3G_FOR;
		LICMP_2[6][1] = 1.- L_C3G_FOR;
		LICMP_2[7][1] = 1.- L_C3G_FOR;
		LICMP_2[0][2] = 1.- L_SAV_C4G;
		LICMP_2[1][2] = 1.- L_FOR_C4G;
		LICMP_2[2][2] = 1.- L_C4G_C4G;
		LICMP_2[3][2] = 1.- L_C4G_C4G;
		LICMP_2[4][2] = 1.- L_C4G_C4G;
		LICMP_2[5][2] = 1.- L_C3G_C4G;
		LICMP_2[6][2] = 1.- L_C3G_C4G;
		LICMP_2[7][2] = 1.- L_C3G_C4G;
		LICMP_2[0][3] = 1.- L_SAV_C4G;
		LICMP_2[1][3] = 1.- L_FOR_C4G;
		LICMP_2[2][3] = 1.- L_C4G_C4G;
		LICMP_2[3][3] = 1.- L_C4G_C4G;
		LICMP_2[4][3] = 1.- L_C4G_C4G;
		LICMP_2[5][3] = 1.- L_C3G_C4G;
		LICMP_2[6][3] = 1.- L_C3G_C4G;
		LICMP_2[7][3] = 1.- L_C3G_C4G;
		LICMP_2[0][4] = 1.- L_SAV_C4G;
		LICMP_2[1][4] = 1.- L_FOR_C4G;
		LICMP_2[2][4] = 1.- L_C4G_C4G;
		LICMP_2[3][4] = 1.- L_C4G_C4G;
		LICMP_2[4][4] = 1.- L_C4G_C4G;
		LICMP_2[5][4] = 1.- L_C3G_C4G;
		LICMP_2[6][4] = 1.- L_C3G_C4G;
		LICMP_2[7][4] = 1.- L_C3G_C4G;
		LICMP_2[0][5] = 1.- L_SAV_C3G;
		LICMP_2[1][5] = 1.- L_FOR_C3G;
		LICMP_2[2][5] = 1.- L_C4G_C3G;
		LICMP_2[3][5] = 1.- L_C4G_C3G;
		LICMP_2[4][5] = 1.- L_C4G_C3G;
		LICMP_2[5][5] = 1.- L_C3G_C3G;
		LICMP_2[6][5] = 1.- L_C3G_C3G;
		LICMP_2[7][5] = 1.- L_C3G_C3G;
		LICMP_2[0][6] = 1.- L_SAV_C3G;
		LICMP_2[1][6] = 1.- L_FOR_C3G;
		LICMP_2[2][6] = 1.- L_C4G_C3G;
		LICMP_2[3][6] = 1.- L_C4G_C3G;
		LICMP_2[4][6] = 1.- L_C4G_C3G;
		LICMP_2[5][6] = 1.- L_C3G_C3G;
		LICMP_2[6][6] = 1.- L_C3G_C3G;
		LICMP_2[7][6] = 1.- L_C3G_C3G;
		LICMP_2[0][7] = 1.- L_SAV_C3G;
		LICMP_2[1][7] = 1.- L_FOR_C3G;
		LICMP_2[2][7] = 1.- L_C4G_C3G;
		LICMP_2[3][7] = 1.- L_C4G_C3G;
		LICMP_2[4][7] = 1.- L_C4G_C3G;
		LICMP_2[5][7] = 1.- L_C3G_C3G;
		LICMP_2[6][7] = 1.- L_C3G_C3G;
		LICMP_2[7][7] = 1.- L_C3G_C3G;
		

	}
	else
	{
		for ( int i=0; i<INPUT_PARAMS; i++ ) s_params[i] = (string) argv[9+i];
		for ( int i=0; i<INPUT_PARAMS; i++ ) d_params[i] = atof(argv[9+i]);
		
		CLD_TREE                 = d_params[0]*0.1/100.;
		K_CAN_EXT_TREE[0]        = d_params[1]*1./100.;
		K_CAN_EXT_TREE[1]        = d_params[2]*1./100.;
		K_CAN_EXT_TREE_INVERSE[0] = 1./K_CAN_EXT_TREE[0];
		K_CAN_EXT_TREE_INVERSE[1] = 1./K_CAN_EXT_TREE[1];
		IGNITION_PROB            = d_params[3]*0.1/100.;
		IGNITION_PAR_2           = d_params[4]*1./100.;
		UPSILON_STEM_TREE        = d_params[5]*1000./100.;
		UPSILON_ROOT_TREE        = UPSILON_STEM_TREE*6./15.;
		SIGMA_GROW_RESP_TREE[0]  = d_params[6]*1./100.;
		SIGMA_GROW_RESP_TREE[1]  = d_params[7]*1./100.;
		SIGMA_GROW_RESP_GRASS[0]    = SIGMA_GROW_RESP_TREE[0];
		SIGMA_GROW_RESP_GRASS[1]    = SIGMA_GROW_RESP_TREE[1];
		PROB_ROOT_SUCKER[0]         = d_params[8]*1./100.;
		PROB_ROOT_SUCKER[1]         = d_params[9]*1./100.;
		BETA_STEM_TREE           = d_params[10]*0.1/100.;
		BETA_ROOT_TREE           = BETA_STEM_TREE;
		SEED_GERM_PROB[0]           = d_params[11]*1./100.;
		SEED_GERM_PROB[1]           = d_params[12]*1./100.;
		DEATH_PROB_FROST[0]         = d_params[13]*0.01/100.;
		DEATH_PROB_FROST[1]         = d_params[14]*0.01/100.;
		DEATH_PROB_CARBON[0]        = d_params[15]*0.01/100.;
		DEATH_PROB_CARBON[1]        = d_params[16]*0.01/100.;
		DEATH_PROB_COMP[0]          = d_params[17]*0.01/100.;
		DEATH_PROB_COMP[1]          = d_params[18]*0.001/100.;
		LIGHT_COMP_SAVSAV   = d_params[19]*1./100.;
		LIGHT_COMP_SAVC4    = d_params[20]*1./100.;
		LIGHT_COMP_SAVC3    = d_params[21]*1./100.;
		LIGHT_COMP_FORSAV   = d_params[22]*1./100.;
		TOP_KILL_CONST[0]           = d_params[23]*1./10; 
		TOP_KILL_CONST[1]           = d_params[24]*1./10;  
		TOP_KILL_H[0]               = d_params[25]*1./1000;  
		TOP_KILL_H[1]               = d_params[26]*1./1000;  
		TOP_KILL_I[0]               = d_params[27]*1./1000000;   
		TOP_KILL_I[1]               = d_params[28]*1./1000000;  
		
		
		L_SAV_SAV = d_params[19]*1./100.;
		L_SAV_FOR = 0.15;
		L_SAV_C4G = d_params[20]*1./100.;
		L_SAV_C3G = d_params[21]*1./100.;

		L_FOR_SAV = d_params[22]*1./100.;
		L_FOR_FOR = 0.15;
		L_FOR_C4G = 0.5;
		L_FOR_C3G = 0.15;

		L_C4G_SAV = 0.5;
		L_C4G_FOR = 0.15;
		L_C4G_C4G = 0.5;
		L_C4G_C3G = 0.15;

		L_C3G_SAV = 0.5;
		L_C3G_FOR = 0.15;
		L_C3G_C4G = 0.5;
		L_C3G_C3G = 0.15;

		LICMP_1[0][0] = L_SAV_SAV;
		LICMP_1[1][0] = L_FOR_SAV;
		LICMP_1[2][0] = L_C4G_SAV;
		LICMP_1[3][0] = L_C4G_SAV;
		LICMP_1[4][0] = L_C4G_SAV;
		LICMP_1[5][0] = L_C3G_SAV;
		LICMP_1[6][0] = L_C3G_SAV;
		LICMP_1[7][0] = L_C3G_SAV;
		LICMP_1[0][1] = L_SAV_FOR;
		LICMP_1[1][1] = L_FOR_FOR;
		LICMP_1[2][1] = L_C4G_FOR;
		LICMP_1[3][1] = L_C4G_FOR;
		LICMP_1[4][1] = L_C4G_FOR;
		LICMP_1[5][1] = L_C3G_FOR;
		LICMP_1[6][1] = L_C3G_FOR;
		LICMP_1[7][1] = L_C3G_FOR;
		LICMP_1[0][2] = L_SAV_C4G;
		LICMP_1[1][2] = L_FOR_C4G;
		LICMP_1[2][2] = L_C4G_C4G;
		LICMP_1[3][2] = L_C4G_C4G;
		LICMP_1[4][2] = L_C4G_C4G;
		LICMP_1[5][2] = L_C3G_C4G;
		LICMP_1[6][2] = L_C3G_C4G;
		LICMP_1[7][2] = L_C3G_C4G;
		LICMP_1[0][3] = L_SAV_C4G;
		LICMP_1[1][3] = L_FOR_C4G;
		LICMP_1[2][3] = L_C4G_C4G;
		LICMP_1[3][3] = L_C4G_C4G;
		LICMP_1[4][3] = L_C4G_C4G;
		LICMP_1[5][3] = L_C3G_C4G;
		LICMP_1[6][3] = L_C3G_C4G;
		LICMP_1[7][3] = L_C3G_C4G;
		LICMP_1[0][4] = L_SAV_C4G;
		LICMP_1[1][4] = L_FOR_C4G;
		LICMP_1[2][4] = L_C4G_C4G;
		LICMP_1[3][4] = L_C4G_C4G;
		LICMP_1[4][4] = L_C4G_C4G;
		LICMP_1[5][4] = L_C3G_C4G;
		LICMP_1[6][4] = L_C3G_C4G;
		LICMP_1[7][4] = L_C3G_C4G;
		LICMP_1[0][5] = L_SAV_C3G;
		LICMP_1[1][5] = L_FOR_C3G;
		LICMP_1[2][5] = L_C4G_C3G;
		LICMP_1[3][5] = L_C4G_C3G;
		LICMP_1[4][5] = L_C4G_C3G;
		LICMP_1[5][5] = L_C3G_C3G;
		LICMP_1[6][5] = L_C3G_C3G;
		LICMP_1[7][5] = L_C3G_C3G;
		LICMP_1[0][6] = L_SAV_C3G;
		LICMP_1[1][6] = L_FOR_C3G;
		LICMP_1[2][6] = L_C4G_C3G;
		LICMP_1[3][6] = L_C4G_C3G;
		LICMP_1[4][6] = L_C4G_C3G;
		LICMP_1[5][6] = L_C3G_C3G;
		LICMP_1[6][6] = L_C3G_C3G;
		LICMP_1[7][6] = L_C3G_C3G;
		LICMP_1[0][7] = L_SAV_C3G;
		LICMP_1[1][7] = L_FOR_C3G;
		LICMP_1[2][7] = L_C4G_C3G;
		LICMP_1[3][7] = L_C4G_C3G;
		LICMP_1[4][7] = L_C4G_C3G;
		LICMP_1[5][7] = L_C3G_C3G;
		LICMP_1[6][7] = L_C3G_C3G;
		LICMP_1[7][7] = L_C3G_C3G;
		
		LICMP_2[0][0] = 1.- L_SAV_SAV;
		LICMP_2[1][0] = 1.- L_FOR_SAV;
		LICMP_2[2][0] = 1.- L_C4G_SAV;
		LICMP_2[3][0] = 1.- L_C4G_SAV;
		LICMP_2[4][0] = 1.- L_C4G_SAV;
		LICMP_2[5][0] = 1.- L_C3G_SAV;
		LICMP_2[6][0] = 1.- L_C3G_SAV;
		LICMP_2[7][0] = 1.- L_C3G_SAV;
		LICMP_2[0][1] = 1.- L_SAV_FOR;
		LICMP_2[1][1] = 1.- L_FOR_FOR;
		LICMP_2[2][1] = 1.- L_C4G_FOR;
		LICMP_2[3][1] = 1.- L_C4G_FOR;
		LICMP_2[4][1] = 1.- L_C4G_FOR;
		LICMP_2[5][1] = 1.- L_C3G_FOR;
		LICMP_2[6][1] = 1.- L_C3G_FOR;
		LICMP_2[7][1] = 1.- L_C3G_FOR;
		LICMP_2[0][2] = 1.- L_SAV_C4G;
		LICMP_2[1][2] = 1.- L_FOR_C4G;
		LICMP_2[2][2] = 1.- L_C4G_C4G;
		LICMP_2[3][2] = 1.- L_C4G_C4G;
		LICMP_2[4][2] = 1.- L_C4G_C4G;
		LICMP_2[5][2] = 1.- L_C3G_C4G;
		LICMP_2[6][2] = 1.- L_C3G_C4G;
		LICMP_2[7][2] = 1.- L_C3G_C4G;
		LICMP_2[0][3] = 1.- L_SAV_C4G;
		LICMP_2[1][3] = 1.- L_FOR_C4G;
		LICMP_2[2][3] = 1.- L_C4G_C4G;
		LICMP_2[3][3] = 1.- L_C4G_C4G;
		LICMP_2[4][3] = 1.- L_C4G_C4G;
		LICMP_2[5][3] = 1.- L_C3G_C4G;
		LICMP_2[6][3] = 1.- L_C3G_C4G;
		LICMP_2[7][3] = 1.- L_C3G_C4G;
		LICMP_2[0][4] = 1.- L_SAV_C4G;
		LICMP_2[1][4] = 1.- L_FOR_C4G;
		LICMP_2[2][4] = 1.- L_C4G_C4G;
		LICMP_2[3][4] = 1.- L_C4G_C4G;
		LICMP_2[4][4] = 1.- L_C4G_C4G;
		LICMP_2[5][4] = 1.- L_C3G_C4G;
		LICMP_2[6][4] = 1.- L_C3G_C4G;
		LICMP_2[7][4] = 1.- L_C3G_C4G;
		LICMP_2[0][5] = 1.- L_SAV_C3G;
		LICMP_2[1][5] = 1.- L_FOR_C3G;
		LICMP_2[2][5] = 1.- L_C4G_C3G;
		LICMP_2[3][5] = 1.- L_C4G_C3G;
		LICMP_2[4][5] = 1.- L_C4G_C3G;
		LICMP_2[5][5] = 1.- L_C3G_C3G;
		LICMP_2[6][5] = 1.- L_C3G_C3G;
		LICMP_2[7][5] = 1.- L_C3G_C3G;
		LICMP_2[0][6] = 1.- L_SAV_C3G;
		LICMP_2[1][6] = 1.- L_FOR_C3G;
		LICMP_2[2][6] = 1.- L_C4G_C3G;
		LICMP_2[3][6] = 1.- L_C4G_C3G;
		LICMP_2[4][6] = 1.- L_C4G_C3G;
		LICMP_2[5][6] = 1.- L_C3G_C3G;
		LICMP_2[6][6] = 1.- L_C3G_C3G;
		LICMP_2[7][6] = 1.- L_C3G_C3G;
		LICMP_2[0][7] = 1.- L_SAV_C3G;
		LICMP_2[1][7] = 1.- L_FOR_C3G;
		LICMP_2[2][7] = 1.- L_C4G_C3G;
		LICMP_2[3][7] = 1.- L_C4G_C3G;
		LICMP_2[4][7] = 1.- L_C4G_C3G;
		LICMP_2[5][7] = 1.- L_C3G_C3G;
		LICMP_2[6][7] = 1.- L_C3G_C3G;
		LICMP_2[7][7] = 1.- L_C3G_C3G;

	}
#   endif
	
// 	write input setting to console
	DEBUG( DEBUG_INPUT_DATA,
		endl <<
		"INI -------------------------------------------------------------------------------------------------------------------------------" << endl <<
		"INI  INPUT PARAMETERS FOR SIMULATIONS" << endl <<
		"INI        longitude (x)                 " << longitude << endl <<
		"INI        latitude (y)                  " << latitude << endl <<
		"INI        years_to_run                  " << years_to_run << endl <<
		"INI        num_of_trees                  " << num_of_trees << endl <<
		"INI        with_fire                     " << with_fire << endl <<
		"INI        seed for random numbers       " << d_seed << endl <<
		"INI        scenario                      " << i_scen << endl <<
		"INI        output directory              " << output_directory << endl <<
		"INI -------------------------------------------------------------------------------------------------------------------------------")
#ifdef S_RAIN_FROM_FILE
// 	Get file name with rainfall data and open this file
	string FILE_file_name = IN_DATA_HOME + (string)"prec_sites/prec_" + (string)argv[3] + "__" + 
							(string)argv[4] + "__" + (string)argv[6] + ".dat";
	
	DEBUG( DEBUG_INPUT_DATA,
		"INI  USE RAINFALL FROM FILE              " << FILE_file_name << endl <<
		"INI -------------------------------------------------------------------------------------------------------------------------------")
	
	ifstream FILE_prec_file(FILE_file_name.c_str());
	if (!FILE_prec_file) {
		cout << endl << "INI ERROR: Cannot open " << FILE_file_name << ". Exit simulation." << endl << endl;
		return 1;
	}
#endif
#   ifdef GRAZING
    grazing_rate = (double) atof( argv[7] )/50.;
	DEBUG( DEBUG_INPUT_DATA,
		"INI  USE GRAZING RATE                    " << grazing_rate << endl <<
		"INI -------------------------------------------------------------------------------------------------------------------------------")
#   endif
#   ifdef S_CLIM_SENSITIVITY
	DEBUG( DEBUG_INPUT_DATA,
		"INI  USE MANIPULATED CLIMATE             " << endl <<
		"INI        CO2 concentration             " << d_params[0] << endl <<
		"INI        Temperature factor            " << d_params[1] << endl <<
		"INI        Precipitation, mean           " << d_params[2] << endl <<
		"INI        Precipitation, seasonality    " << d_params[3] << endl <<
// 		"INI        Fire occurance                " << d_params[4] << endl <<
		"INI -------------------------------------------------------------------------------------------------------------------------------")
#   endif
#ifdef S_FIRE_FROM_FILE
// 	Get file name with fire regime and open this file
	double FIRE_year;
	double FIRE_day;
	double FIRE_intensity;
	string FIRE_file_name = IN_DATA_HOME + (string)"fire_sites/fire_" + (string)argv[3] + "__" + 
							(string)argv[4] + "__" + (string)argv[6] + ".dat";
	
	DEBUG( DEBUG_INPUT_DATA,
		"INI  USE FIRE SCENARIO FROM FILE         " << FIRE_file_name << endl <<
		"INI -------------------------------------------------------------------------------------------------------------------------------")
	
	ifstream FIRE_file(FIRE_file_name.c_str());
	if (!FIRE_file) {
		cout << endl << "INI ERROR: Cannot open " << FIRE_file_name << ". Exit simulation." << endl << endl;
		return 1;
	}
	FIRE_file >> FIRE_year >> FIRE_day >> FIRE_intensity;   // get first fire
#endif

#ifdef S_ELEPHANTS
	DEBUG( DEBUG_INPUT_DATA,
				   "INI ELEPHANT MODEL                       " << endl <<
				   "INI        visit frequency               " << setw(7) << d_params[0]  << " times per year" << endl <<
				   "INI        utilization intensity         " << setw(7) << d_params[1]  << " times the biomass" << endl <<
				   "INI        bark stripping, mortality     " << setw(7) << d_params[2]  << " fraction of trees" << endl <<
				   "INI        maximum tree consumption      " << setw(7) << d_params[3]  << " fraction of tree biomass" << endl <<
				   "INI        proportion of males           " << setw(7) << d_params[4]  << " " << endl <<
				   "INI        wastage at utilization        " << setw(7) << d_params[5]  << " fraction" << endl <<
				   "INI        elephant number               " << setw(7) << d_params[6]  << " " << endl <<
				   "INI        diet partitioning             " << setw(7) << d_params[7]  << " " << endl <<
				   "INI        daily consumption per el.     " << setw(7) << d_params[8]  << " kg per day" << endl <<
				   "INI        grazing                       " << setw(7) << d_params[9]  << " kg per day" << endl <<
				   "INI        area of park                  " << setw(7) << d_params[10] << " km^2" << endl <<
				   "INI -------------------------------------------------------------------------------------------------------------------------------")
#endif
	
#ifdef OPTIM_GLOBALS
	DEBUG( DEBUG_INPUT_DATA,
				   "INI  OPTIMIZE GLOBAL PARAMETERS          " << endl <<
				   "INI        CLD_TREE                      "  << CLD_TREE << endl <<
				   "INI        K_CAN_EXT_TREE                " << K_CAN_EXT_TREE[0] << K_CAN_EXT_TREE[1] << endl <<
				   "INI        IGNITION_PROB                 " << IGNITION_PROB << endl <<
				   "INI        IGNITION_PAR_2                " << IGNITION_PAR_2 << endl <<
				   "INI        UPSILON_STEM_TREE             " << UPSILON_STEM_TREE << endl <<
				   "INI        UPSILON_ROOT_TREE             " << UPSILON_ROOT_TREE << endl <<
				   "INI        SIGMA_GROW_RESP_TREE          " << SIGMA_GROW_RESP_TREE[0] << SIGMA_GROW_RESP_TREE[1] << endl <<
				   "INI        SIGMA_GROW_RESP_GRASS         " << SIGMA_GROW_RESP_GRASS[0] << SIGMA_GROW_RESP_GRASS[1] << endl <<
				   "INI        PROB_ROOT_SUCKER              " << PROB_ROOT_SUCKER[0] << PROB_ROOT_SUCKER[1] << endl <<
				   "INI        BETA_STEM_TREE                " << BETA_STEM_TREE << endl <<
				   "INI        BETA_ROOT_TREE                " << BETA_ROOT_TREE << endl <<
				   "INI        SEED_GERM_PROB                " << SEED_GERM_PROB[0] << SEED_GERM_PROB[1] << endl <<
				   "INI        DEATH_PROB_FROST              " << DEATH_PROB_FROST[0] << DEATH_PROB_FROST[1] << endl <<
				   "INI        DEATH_PROB_CARBON             " << DEATH_PROB_CARBON[0] << DEATH_PROB_CARBON[1] << endl <<
				   "INI        DEATH_PROB_COMP               " << DEATH_PROB_COMP[0] << DEATH_PROB_COMP[1] << endl <<
				   "INI        LIGHT_COMP_TREE_TREE_1        " << LIGHT_COMP_TREE_TREE_1 << endl <<
				   "INI        LIGHT_COMP_GRASS_TREE_1       " << LIGHT_COMP_GRASS_TREE_1 << endl <<
				   "INI        LIGHT_COMP_GRASS_GRASS_1      " << LIGHT_COMP_GRASS_GRASS_1 << endl <<
				   "INI        LIGHT_COMP_TREE_GRASS_1       " << LIGHT_COMP_TREE_GRASS_1 << endl <<
				   "INI -------------------------------------------------------------------------------------------------------------------------------")
#endif
	
	double ca_par_preassure = CA_PAR_PREASSURE;
	
	// read in databases/shortlist and initialize data
	clInDataReader MyReader;
	clInData IData = MyReader.getInData( longitude, latitude );
	
//     for ( int kkk=0; kkk<12; kkk++ ) IData.reh_[kkk] *= 2.;
// 	cout << IData.soil_N_ << "  " << IData.soil_C_ << "  " << endl;
    
	if ( IData.soil_N_<0 || IData.reh_[0]<0 || IData.theta_wp_[0]<0 )
	{
		cout << "INI ERROR: Cannot simulate the ocean, abort simulation: " << longitude << " " << latitude << endl << endl;
		return 1;
	}
    
    double soil_carbon = IData.soil_C_[0];
//     clSoil mySoil( soil_carbon/1000. ); // initialize humus pool with soil C, from g/m^2 to kg/m^2
    clSoil mySoil( soil_carbon/5000. ); // initialize humus pool with soil C, from g/m^2 to kg/m^2

	
#   ifdef CLIM_SCEN
	DEBUG( DEBUG_INPUT_DATA,
				   "INI  RUN CLIMATE SCENARIO                " << i_scen)

	double cs_mean_prec;
	double cs_rbeta_orig[12];
	ifstream cs_co2_file;
	ifstream cs_tmp_file;
	ifstream cs_prec_file;
	string cs_co2_file_name;
	string cs_tmp_file_name;
	string cs_pre_file_name;
	
	for ( int cs_cnt=0; cs_cnt<12; cs_cnt++ ) cs_rbeta_orig[cs_cnt] = IData.rbeta_[cs_cnt];
	
	if      ( i_scen== 0 ) cs_co2_file_name = IN_DATA_HOME + "ClimateChange/co2_tod.txt";    // Constant co2, tmp and rain
	else if ( i_scen== 1 ) cs_co2_file_name = IN_DATA_HOME + "ClimateChange/co2_A1B.txt";    // SRES A1B with global mean CO2 and temp
	else if ( i_scen== 2 ) cs_co2_file_name = IN_DATA_HOME + "ClimateChange/co2_A2.txt";     // SRES A2 with CO2, temp and rain
	else if ( i_scen== 3 ) cs_co2_file_name = IN_DATA_HOME + "ClimateChange/co2_A1B.txt";    // SRES A1B with CO2, temp and rain
	else if ( i_scen== 4 ) cs_co2_file_name = IN_DATA_HOME + "ClimateChange/co2_B1.txt";     // SRES B1 with CO2, temp and rain
	else if ( i_scen== 5 ) cs_co2_file_name = IN_DATA_HOME + "ClimateChange/co2_A2.txt";     // SRES A2 with CO2, temp, without rain
	else if ( i_scen== 6 ) cs_co2_file_name = IN_DATA_HOME + "ClimateChange/co2_A1B.txt";    // SRES A1B with CO2, temp, without rain
	else if ( i_scen== 7 ) cs_co2_file_name = IN_DATA_HOME + "ClimateChange/co2_B1.txt";     // SRES B1 with CO2, temp, without rain
// 	else if ( i_scen==10 ) cs_co2_file_name = IN_DATA_HOME + "ClimateChange/vostok_co2_all.dat";    // SRES B1 with CO2, temp, without rain
	else if ( i_scen==10 ) cs_co2_file_name = IN_DATA_HOME + "ClimateChange/co2_yyy.txt";    // hysteresis loop for co2
	else if ( i_scen==11 ) cs_co2_file_name = IN_DATA_HOME + "ClimateChange/co2_yyy.txt";    // hysteresis loop for co2
	else if ( i_scen==20 ) cs_co2_file_name = IN_DATA_HOME + "ClimateChange/co2_xxx.txt";    // constant co2 at 387ppm
// 	else if ( i_scen==20 ) cs_co2_file_name = IN_DATA_HOME + "ClimateChange/co2_550.txt";
	else if ( i_scen==21 ) cs_co2_file_name = IN_DATA_HOME + "ClimateChange/co2_xxx.txt";    // constant co2 at 387ppm
//	else if ( i_scen==30 ) cs_co2_file_name = IN_DATA_HOME + "ClimateChange/co2_550.txt";    // constant co2 at 387ppm
 	else if ( i_scen==30 ) cs_co2_file_name = IN_DATA_HOME + "ClimateChange/co2_xxx.txt";    // constant co2 at 387ppm
			
	if ( i_scen==11 ) IData = ManipulateSeasonality( -5, IData.ralpha_, IData.rbeta_, IData.pwet_, IData.rdo_, IData );
	if ( i_scen==21 ) IData = ManipulateSeasonality( -5, IData.ralpha_, IData.rbeta_, IData.pwet_, IData.rdo_, IData );
	
	
	DEBUG( DEBUG_INPUT_DATA,
		"INI        CO2 file                      " << cs_co2_file_name )
	
	if      ( i_scen== 0 ) cs_tmp_file_name = IN_DATA_HOME + "ClimateChange/tmp_tod.txt";
	else if ( i_scen== 1 ) cs_tmp_file_name = IN_DATA_HOME + "ClimateChange/tmp_A1B.txt";
	else if ( i_scen== 7 ) cs_tmp_file_name = IN_DATA_HOME + "ClimateChange/tmp_B1.txt";
	else if ( i_scen==10 ) cs_tmp_file_name = IN_DATA_HOME + "ClimateChange/tmp_yyy.txt";   // constant temperature
	else if ( i_scen==11 ) cs_tmp_file_name = IN_DATA_HOME + "ClimateChange/tmp_yyy.txt";
	else if ( i_scen==20 ) cs_tmp_file_name = IN_DATA_HOME + "ClimateChange/tmp_xxx.txt";   // hysteresis loop for temperature
	else if ( i_scen==21 ) cs_tmp_file_name = IN_DATA_HOME + "ClimateChange/tmp_xxx.txt";   // hysteresis loop for temperature
	else if ( i_scen==30 ) cs_tmp_file_name = IN_DATA_HOME + "ClimateChange/tmp_yyy.txt";   // constant temperature
	else if ( i_scen== 2 || i_scen==3 || i_scen==4 || i_scen==5 || i_scen==6 )
	{
		int    cs_y_cnt = 0;
// 		double cs_x_coo = 52.5;
		double cs_x_coo = 168.75;
		double cs_y_coo = 38.238;
		
		while( fabs( cs_x_coo-longitude )>1.875/2. )   cs_x_coo -= 1.875;
		for( int cs_cnt=0; cs_cnt<42; cs_cnt++ )
		{
			if( fabs( cs_y_vec[cs_cnt]-latitude )<=0.935 ) cs_y_coo = cs_y_vec[cs_cnt];
		}
		
		char cs_tmp_x[10]; sprintf(cs_tmp_x, "%.3f", cs_x_coo);
		char cs_tmp_y[10]; sprintf(cs_tmp_y, "%.3f", cs_y_coo);
		
		if ( i_scen== 2 ) cs_tmp_file_name = IN_DATA_HOME + (string)"ClimateChange/temp_files/SRA2_tas_"
		                                   + (string)cs_tmp_x + "_" + (string)cs_tmp_y + ".dat";
		if ( i_scen== 3 ) cs_tmp_file_name = IN_DATA_HOME + (string)"ClimateChange/temp_files/SRA1B_tas_"
		                                   + (string)cs_tmp_x + "_" + (string)cs_tmp_y + ".dat";
		if ( i_scen== 4 ) cs_tmp_file_name = IN_DATA_HOME + (string)"ClimateChange/temp_files/SRB1_tas_"
		                                   + (string)cs_tmp_x + "_" + (string)cs_tmp_y + ".dat";
		if ( i_scen== 5 ) cs_tmp_file_name = IN_DATA_HOME + (string)"ClimateChange/temp_files/SRA2_tas_"
		                                   + (string)cs_tmp_x + "_" + (string)cs_tmp_y + ".dat";
		if ( i_scen== 6 ) cs_tmp_file_name = IN_DATA_HOME + (string)"ClimateChange/temp_files/SRA1B_tas_"
		                                   + (string)cs_tmp_x + "_" + (string)cs_tmp_y + ".dat";
	}
	
	DEBUG( DEBUG_INPUT_DATA,
		   "INI        Temperature file              " << cs_tmp_file_name )
	
	if      ( i_scen== 0 ) cs_pre_file_name = IN_DATA_HOME + "ClimateChange/pre_tod.txt";
	else if ( i_scen== 1 ) cs_pre_file_name = IN_DATA_HOME + "ClimateChange/pre_tod.txt";
	else if ( i_scen== 5 ) cs_pre_file_name = IN_DATA_HOME + "ClimateChange/pre_tod.txt";
	else if ( i_scen== 6 ) cs_pre_file_name = IN_DATA_HOME + "ClimateChange/pre_tod.txt";
	else if ( i_scen== 7 ) cs_pre_file_name = IN_DATA_HOME + "ClimateChange/pre_tod.txt";
	else if ( i_scen==10 ) cs_pre_file_name = IN_DATA_HOME + "ClimateChange/pre_yyy.txt";
	else if ( i_scen==11 ) cs_pre_file_name = IN_DATA_HOME + "ClimateChange/pre_yyy.txt";
	else if ( i_scen==20 ) cs_pre_file_name = IN_DATA_HOME + "ClimateChange/pre_xxx.txt";
	else if ( i_scen==21 ) cs_pre_file_name = IN_DATA_HOME + "ClimateChange/pre_xxx.txt";
	else if ( i_scen==30 ) cs_pre_file_name = IN_DATA_HOME + "ClimateChange/pre_zzz.txt";
	else if ( i_scen== 2 || i_scen==3 || i_scen==4 )
	{
		int    cs_y_cnt = 0;
// 		double cs_x_coo = 52.5;
		double cs_x_coo = 168.75;
		double cs_y_coo = 38.238;
		
		while( fabs( cs_x_coo-longitude )>1.875/2. )   cs_x_coo -= 1.875;
		for( int cs_cnt=0; cs_cnt<42; cs_cnt++ )
		{
			if( fabs( cs_y_vec[cs_cnt]-latitude )<=0.935 ) cs_y_coo = cs_y_vec[cs_cnt];
		}
		
		char cs_tmp_x[10]; sprintf(cs_tmp_x, "%.3f", cs_x_coo);
		char cs_tmp_y[10]; sprintf(cs_tmp_y, "%.3f", cs_y_coo);
		
		if ( i_scen==2 ) cs_pre_file_name = IN_DATA_HOME + (string)"ClimateChange/rain_files/SRA2_pr_"
											+ (string)cs_tmp_x + "_" + (string)cs_tmp_y + ".dat";
		if ( i_scen==3 ) cs_pre_file_name = IN_DATA_HOME + (string)"ClimateChange/rain_files/SRA1B_pr_"
											+ (string)cs_tmp_x + "_" + (string)cs_tmp_y + ".dat";
		if ( i_scen==4 ) cs_pre_file_name = IN_DATA_HOME + (string)"ClimateChange/rain_files/SRB1_pr_"
											+ (string)cs_tmp_x + "_" + (string)cs_tmp_y + ".dat";
	}
	
	DEBUG( DEBUG_INPUT_DATA,
		   "INI        Precipitation file            " << cs_pre_file_name )
	
	
	// open files
	
	cs_co2_file.open(cs_co2_file_name.c_str());
	
	if ( !cs_co2_file )
	{
		cout << "INI ERROR: Cannot open CO2 file " << cs_co2_file_name << ", abort simulation." << endl << endl;
		return 1;
	}
	
	
	cs_tmp_file.open(cs_tmp_file_name.c_str());
	if (!cs_tmp_file)
	{
		cout << "INI ERROR: Cannot open temperature file " << cs_tmp_file_name << ", abort simulation." << endl << endl;
		return 1;
	}
	
	
	cs_prec_file.open(cs_pre_file_name.c_str());
	if (!cs_prec_file)
	{
		cout << "INI ERROR: Cannot open precipitation file " << cs_pre_file_name << ", abort simulation." << endl << endl;
		return 1;
	}
	
	
	
	double cs_preassure;
	double cs_tmp_inc;
	
	cs_co2_file >> cs_preassure;
	ca_par_preassure = cs_preassure/10.;
	
	cs_tmp_file >> cs_tmp_inc;
	
	arry12 cs_tmp_orig;
	arry12 cs_tmp_day_orig;
	arry12 cs_tmp_min_orig;
	arry12 cs_tmp_max_orig;
	
	for ( int kk=0; kk<12; kk++ )
	{
		cs_tmp_orig[kk]     = IData.tmp_[kk];
		cs_tmp_day_orig[kk] = IData.tmp_day_[kk];
		cs_tmp_min_orig[kk] = IData.tmp_min_[kk];
		cs_tmp_max_orig[kk] = IData.tmp_max_[kk];
	}
	DEBUG( DEBUG_INPUT_DATA,
		"INI -------------------------------------------------------------------------------------------------------------------------------" << endl)

#   endif  //  CLIM_SCEN
	
	
	#ifdef S_CLIM_SENSITIVITY
	ca_par_preassure = d_params[0]/10.;
	
	for ( int cs_cnt=0; cs_cnt<12; cs_cnt++ ) {
		IData.tmp_[cs_cnt]     += d_params[1];
		IData.tmp_day_[cs_cnt] += d_params[1];
		IData.tmp_min_[cs_cnt] += d_params[1];
		IData.tmp_max_[cs_cnt] += d_params[1];
		IData.rbeta_[cs_cnt]   *= d_params[2];
	}
	
	for ( int i=0; i<12; i++ )  // original values
	{
		IData.ralpha_[i] += drand48()*1e-10;  // this ensures that all values are different,
		IData.rbeta_[i]  += drand48()*1e-10;  // this is required for the manipulation algorithm
		IData.pwet_[i]   += drand48()*1e-10;
	}
	
	IData = ManipulateSeasonality( d_params[3], IData.ralpha_, IData.rbeta_, IData.pwet_, IData.rdo_, IData );
	
	#endif
	
	#ifdef S_MANIP_SEASONALITY
	double MAN_SEAS_value = 0.;
	double MAN_SEAS_increment = 0.0007;
	double MAN_SEAS_ralpha[12];
	double MAN_SEAS_rbeta[12];
	double MAN_SEAS_pwet[12];
	double MAN_SEAS_rdo[12];
	
	for ( int i=0; i<12; i++ )  // original values
	{
		MAN_SEAS_ralpha[i]  = IData.ralpha_[i];
		MAN_SEAS_rbeta[i]   = IData.rbeta_[i];
		MAN_SEAS_pwet[i]    = IData.pwet_[i];
		MAN_SEAS_rdo[i]     = IData.rdo_[i];
		
		MAN_SEAS_ralpha[i] += drand48()*1e-10;  // this ensures that all values are different,
		MAN_SEAS_rbeta[i]  += drand48()*1e-10;  // this is required for the manipulation algorithm
		MAN_SEAS_pwet[i]   += drand48()*1e-10;
		
// 		cout << setw(14) << IData.ralpha_[i] << setw(14) << IData.rbeta_[i] << setw(14) << IData.pwet_[i] << setw(14) << IData.rdo_[i] << endl;
		
	}
	IData = ManipulateSeasonality( MAN_SEAS_value, MAN_SEAS_ralpha, MAN_SEAS_rbeta, MAN_SEAS_pwet, MAN_SEAS_rdo, IData );
	
	#endif   // S_MANIP_SEASONALITY
	
	#ifdef NO_SEASONALITY
	IData = ManipulateSeasonality( -5, IData.ralpha_, IData.rbeta_, IData.pwet_, IData.rdo_, IData );
	#endif
	
	IData.calcAtmospheric( latitude );
	IData.calcLeafPhoto( ca_par_preassure, soil_carbon );
	
	
	// calculate mean A0 and respiration, needed for output
	double A0C3_mean   = 0;
	double A0C4_mean   = 0;
	double Rl_mean     = 0;
	double tmp_mean    = 0;
	double sun_mean    = 0;
	double tmp_min     = 1000;
	double A0_max_C3   = 0;
	double hum_mean    = 0;
	
	for ( int ii=0; ii<12; ii++ )
	{
		if ( IData.A012C3_[ii] > A0_max_C3 )
		{
			A0_max_C3     = IData.A012C3_[ii];
			A0_max_C3_day = ii;
		}
		A0C3_mean += IData.A012C3_[ii];
		A0C4_mean += IData.A012C4_[ii];
		Rl_mean   += IData.RmL12C3_[ii];
		tmp_mean  += IData.tmp_day_[ii];
		tmp_min    = MyMin( IData.tmp_min_[ii], tmp_min );
		sun_mean  += IData.sun_[ii];
		hum_mean  += IData.reh_[ii];
//          cout << IData.A012C3_[ii] << "  " << IData.A012C4_[ii] << endl;
	}
	
	A0_max_C3_day *= 31;
	A0C3_mean     /= 12.;
	A0C4_mean     /= 12.;
	Rl_mean       /= 12.;
	tmp_mean      /= 12.;
	sun_mean      /= 12.;
	hum_mean      /= 12.;
	
	cout << A0C3_mean << " " << A0C4_mean << endl;
	
	GetNetRadiation( latitude, IData.sun_[month_in_year], IData.tmp_max_[month_in_year],
					 IData.tmp_min_[month_in_year], IData.eA12_[month_in_year],
					 day_in_year, &radnet_day, &sunhrs_day );
	
	mmsTOkgd = MMSTOKGD_HELPER*sunhrs_day;
	
	diff_fc_wp = new double[IData.soil_layers_];
	
	for ( int i=0; i<IData.soil_layers_; i++ )
	{
		if ( IData.theta_fc_[i] > 0 && IData.theta_wp_[i] > 0 )
			diff_fc_wp[i] = 1./(IData.theta_fc_[i]-IData.theta_wp_[i]);
		else
			diff_fc_wp[i] = 0;
	}
	
	srand48(1122554162+(int)d_seed);
	#ifdef S_FIXED_RAND_INIT
	srand48(1122554162+1);
	#endif
	
	clTreePop MyTreePop( IData.depth_[IData.soil_layers_-1] );
	clGrassPop MyGrassPop;
	
	int tree_type = TR_SAV;          // NEW CODE
	if ( drand48()<=PROB_FOREST_TREE ) tree_type = TR_FOR;
	MyTreePop.addFirstTree( 100., tree_type );
	
	if ( num_of_trees>1 )
	{
// 		MyTreePop.addTree( 100., 1 );
		for ( int count_trees=2; count_trees<=num_of_trees; count_trees++ )
		{
			tree_type = TR_SAV;
// 			#ifndef S_CLIM_NO_FOREST_TREES
			if ( drand48()<=PROB_FOREST_TREE ) tree_type = TR_FOR;
// 			#endif
			MyTreePop.addTree( 150.*drand48(), tree_type );
		}
	}
	
	MyGrassPop.addGrass( INIT_MASS_GRASS, GR_C4_SAV );  // C4 savanna tree sub-canopy
	MyGrassPop.addGrass( INIT_MASS_GRASS, GR_C4_OPN );  // C4 between canopy
	MyGrassPop.addGrass( INIT_MASS_GRASS, GR_C4_FOR );  // C4 forest tree sub-canopy
	MyGrassPop.addGrass( INIT_MASS_GRASS, GR_C3_SAV );  // C3 savanna tree sub-canopy
	MyGrassPop.addGrass( INIT_MASS_GRASS, GR_C3_OPN );  // C3 between canopy
	MyGrassPop.addGrass( INIT_MASS_GRASS, GR_C3_FOR );  // C3 forest tree sub-canopy
	
	
	
	
#ifdef W_SYSDATA
	string _SYSDATA_output_file = output_directory+"SysData_"+doubleToString(longitude)+"_"+doubleToString(latitude)+"_"+
							doubleToString(with_fire)+"_"+doubleToString(d_seed)+"_"+doubleToString(i_scen)+".dat";

	ofstream _SYSDATA_SysData(_SYSDATA_output_file.c_str());
	if(!_SYSDATA_SysData)
	{
		cerr << "INI ERROR: Can't open file: " << _SYSDATA_output_file << endl;
		exit(1);
	}
	
	
	double _SIZE_STRUC_size_struc[NUM_SIZE_CLASSES];				// vector to store the sizestructure
	for ( int i=0; i<NUM_SIZE_CLASSES; i++ ) _SIZE_STRUC_size_struc[i]=0;
	
	string _SIZE_STRUC_output_file = output_directory+"SizeData_"+doubleToString(longitude)+"_"+doubleToString(latitude)+"_"+
							doubleToString(with_fire)+"_"+doubleToString(d_seed)+"_"+doubleToString(i_scen)+".dat";
	
	ofstream _SIZE_STRUC_SizeStrucFile(_SIZE_STRUC_output_file.c_str());
	if(!_SIZE_STRUC_SizeStrucFile)
	{
		cerr << "INI ERROR: Can't open file: " << _SIZE_STRUC_SizeStrucFile << endl;
		exit(1);
	}
	
	string _FIREDATA_output_file = output_directory+"FireData_"+doubleToString(longitude)+"_"+doubleToString(latitude)+"_"+
							doubleToString(with_fire)+"_"+doubleToString(d_seed)+"_"+doubleToString(i_scen)+".dat";
	
	ofstream _FIREDATA_FireData(_FIREDATA_output_file.c_str());
	if(!_FIREDATA_FireData)
	{
		cerr << "INI ERROR: Can't open file: " << _FIREDATA_output_file << endl;
		exit(1);
	}
	
	string _SOILDATA_output_file = output_directory+"SoilData_"+doubleToString(longitude)+"_"+doubleToString(latitude)+"_"+
							doubleToString(with_fire)+"_"+doubleToString(d_seed)+"_"+doubleToString(i_scen)+".dat";
	
	ofstream _SOILDATA_SoilData(_SOILDATA_output_file.c_str());
	if(!_SOILDATA_SoilData)
	{
		cerr << "INI ERROR: Can't open file: " << _SOILDATA_output_file << endl;
		exit(1);
	}
#endif
	


#   ifdef S_ELEPHANTS
    string _ELFS_output_file;
    
    _ELFS_output_file = output_directory+"ElfData.dat";
    
    ofstream _ELFS_Data(_ELFS_output_file.c_str());
    if(!_ELFS_Data)
    {
        cerr << "INI ERROR: Can't open file: " << _ELFS_output_file << endl;
        exit(1);
    }
#   endif


// --------------------------------------------------------------------------------------------------------------
// --- START YEAR LOOP ------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------
	// loop over all years
	for ( count_years=1; count_years<=years_to_run; count_years++ )
	{
		if ( count_years%25==1 )
			DEBUG(  DEBUG_YEARLY_DATA, "YEAR_LGD " <<
					setw(5)  << "Year" <<
					setw(12) << "Gr_above" <<
					setw(12) << "Gr_below" <<
					setw(12) << "Tr_above" <<
					setw(12) << "Tr_below" <<
					setw( 6) << "Tr_num" <<
					setw( 6) << "Tr_new" <<
					setw( 6) << "Tr_dead" <<
					setw(12) << "Tr_cover" <<
					setw(12) << "Tr_max_hei" <<
					setw(12) << "Tr_mean_hei" <<
					setw(12) << "Tr_b_area" <<
					setw(12) << "MAP" <<
					setw( 6) << "Tr_sav" <<
					setw( 6) << "Tr_for" <<
					setw(12) << "Tr_sav_cov" <<
					setw(12) << "Tr_for_cov" <<
					setw(12) << "C3_C4_rat"
					)
		
		GLOB_YEAR = count_years;
		
// 		IData.calcLeafPhoto( ca_par_preassure, soil_carbon );
		
#		ifdef CLIM_SCEN
		cs_co2_file >> cs_preassure;
		ca_par_preassure = cs_preassure/10.;
		
		IData.calcLeafPhoto( ca_par_preassure, soil_carbon );
		
		A0_max_C3     = IData.A012C3_[0];
		A0_max_C3_day = 0;
		A0C3_mean     = 0;
		A0C4_mean     = 0;
		Rl_mean       = 0;
		tmp_mean      = 0;
		
		cs_tmp_file >> cs_tmp_inc;
		for ( int kk=0; kk<12; kk++ )
		{
			if ( IData.A012C3_[kk] > A0_max_C3 )
			{
				A0_max_C3     = IData.A012C3_[kk];
				A0_max_C3_day = kk;
			}
			
			IData.tmp_[kk]     = cs_tmp_inc + cs_tmp_orig[kk];
			IData.tmp_day_[kk] = cs_tmp_inc + cs_tmp_day_orig[kk];
			IData.tmp_min_[kk] = cs_tmp_inc + cs_tmp_min_orig[kk];
			IData.tmp_max_[kk] = cs_tmp_inc + cs_tmp_max_orig[kk];

			A0C3_mean += IData.A012C3_[kk];
			A0C4_mean += IData.A012C4_[kk];
			Rl_mean   += IData.RmL12C3_[kk];
			tmp_mean  += IData.tmp_day_[kk];
		}
		A0_max_C3_day *= 31;
		A0C3_mean     /= 12.;
		A0C4_mean     /= 12.;
		Rl_mean       /= 12.;
		tmp_mean      /= 12.;
		
// 		if ( ( i_scen==2 || i_scen==3 || i_scen==4 || i_scen>9 ) && count_years>=239 ) {
 		if ( count_years>=239 )
		{
			for ( int cs_cnt=0; cs_cnt<12; cs_cnt++ )
			{
				cs_prec_file >> cs_mean_prec;
// 				cout << cs_mean_prec << endl;
				IData.rbeta_[cs_cnt] = cs_rbeta_orig[cs_cnt]*cs_mean_prec;
			}
 		}
		
		#endif  // CLIM_SCEN
		
		#ifdef S_MANIP_FIRE
		if ( count_years>=  201 && count_years<  550 ) glob_fire_param += 0.01074499;
		if ( count_years>=  551 && count_years<  900 ) glob_fire_param -= 0.01074499;
		if ( count_years>=  901 && count_years< 1250 ) glob_fire_param += 0.01074499;
		if ( count_years>= 1251 && count_years< 1600 ) glob_fire_param -= 0.01074499;
		// 		cout << glob_fire_param << endl;
		#endif
		
		
		#ifdef S_MANIP_SEASONALITY
		if ( count_years>= 201 && count_years< 550 ) MAN_SEAS_value += MAN_SEAS_increment;
		if ( count_years>= 551 && count_years< 900 ) MAN_SEAS_value -= MAN_SEAS_increment;
		if ( count_years>= 901 && count_years<1250 ) MAN_SEAS_value += MAN_SEAS_increment;
		if ( count_years>=1251 && count_years<1600 ) MAN_SEAS_value -= MAN_SEAS_increment;
		IData = ManipulateSeasonality( MAN_SEAS_value, MAN_SEAS_ralpha, MAN_SEAS_rbeta, MAN_SEAS_pwet, MAN_SEAS_rdo, IData );
		#endif
		
		
		// generate daily rain sequence and store it in Rain
		RainFallYear( DAYS_IN_YEAR, IData.rdo_, IData.pwet_, IData.ralpha_, IData.rbeta_, Rain );
		
		
		
#		ifdef S_RAIN_FROM_FILE
		for ( int nnn=0; nnn<365; nnn++ ) {
			FILE_prec_file >> Rain[nnn];
// 			Rain[nnn] *= ZZZ/10000.;
		}
#		endif
		
#      ifdef ELEPHANTS_FUT
        for ( int nnn=0; nnn<365; nnn++ )
			Rain[nnn] *= 0.939627;  // tsavo
// 			Rain[nnn] *= 0.707084;  // chobe
// 			Rain[nnn] *= 0.8529412;  // kruger
		#     endif
		
#		ifdef S_DEF_RAIN
		for ( int a_rain_counter=0; a_rain_counter<365; a_rain_counter++ )
			Rain[a_rain_counter] = 5.*(sin(2.*3.1415/365.*(double)a_rain_counter)+1.);
#		endif
		
		
#       ifdef S_RAIN_MOD_FEEDBACK
        for ( int a_rain_counter=0; a_rain_counter<365; a_rain_counter++ )
        {
            if      (count_years>=1000 & count_years<4000)
                Rain[a_rain_counter] += (0.0002739726*(count_years-1000));
            else if (count_years>=4000 & count_years<5000)
                Rain[a_rain_counter] += (0.8219178);
            else if (count_years>=5000 & count_years<8000)
                Rain[a_rain_counter] += (0.0002739726*(8000-count_years));
            else ;
        }
#       endif
#       ifdef S_RAIN_MOD_VEG_FEEDBACK
        for ( int a_rain_counter=0; a_rain_counter<365; a_rain_counter++ )
        {
            double tmp_tre_c  = MyTreePop.getpCanopy();
            double tmp_gr_bio = tmp_tre_c*MyGrassPop.getLeafBiomassCan()+(1.-tmp_tre_c)*MyGrassPop.getLeafBiomassOpn();

            if      (count_years>=1000 & count_years<4000)   // 0.2739726 assumes that trees and grasses can increase prec. by 400mm/year
                Rain[a_rain_counter] += (0.0002739726*(count_years-1000) + 1.095890*tmp_tre_c + 1.095890*tmp_gr_bio);
            else if (count_years>=4000 & count_years<5000)
                Rain[a_rain_counter] += (0.8219178 + 1.095890*tmp_tre_c + 1.095890*tmp_gr_bio);
            else if (count_years>=5000 & count_years<8000)
                Rain[a_rain_counter] += (0.0002739726*(8000-count_years) + 1.095890*tmp_tre_c + 1.095890*tmp_gr_bio);
            else ;
        }
#       endif

		#ifdef S_ELEPHANTS
		if ( with_fire>=1  && count_years>150 ) GetIgnitions( Ignitions, MyTreePop.getpCanopy() );
		#elif defined S_CLIM_FIRE_FROM_150
		if ( with_fire>=1  && count_years>150 ) GetIgnitions( Ignitions, MyTreePop.getpCanopy() );
		#elif defined S_CLIM_FIRE_FROM_600
		if ( with_fire>=1  && count_years>600 ) GetIgnitions( Ignitions, MyTreePop.getpCanopy() );
		#elif defined S_CLIM_FIRE_FROM_50
		if ( with_fire>=1  && count_years>50 ) GetIgnitions( Ignitions, MyTreePop.getpCanopy() );
		#else
		if ( with_fire==1  && count_years> 30 ) GetIgnitions( Ignitions, MyTreePop.getpCanopy() );
		if ( with_fire==2  && count_years>193 ) GetIgnitions( Ignitions, 5. );   // fire suppression after 1954 EBP experiment
#       endif
		
		MyTreePop.setBornToZero();
		MyTreePop.setDeadToZero();
		
		evapo_sum_year		= 0;
		rain_sum_year		= 0;
		fire_flag           = 0;
		
		// --------------------------------------------------------------------------------------------------------------
		// --- START DAY LOOP -------------------------------------------------------------------------------------------
		// --------------------------------------------------------------------------------------------------------------
		// loop over all days in a year
		for ( day_in_year=0; day_in_year<DAYS_IN_YEAR; day_in_year++ )
		{
			
			month_in_year = (int) floor( ((double) (day_in_year+1))/30.42 );
// 			cout << "PREC " << setw(14) << month_in_year << setw(14) << Rain[day_in_year] << endl;
			
			T_fac = MaintRespTmpFac(IData.tmp_min_[month_in_year]);   // tmp function needed for respiration
			
			tmp_min_month = IData.tmp_min_[month_in_year];
			
			GetNetRadiation( latitude, IData.sun_[month_in_year], IData.tmp_max_[month_in_year],
							 IData.tmp_min_[month_in_year], IData.eA12_[month_in_year],
							 day_in_year, &radnet_day, &sunhrs_day );
 			
// 			radnetsum += radnet_day;
			
			mmsTOkgd = MMSTOKGD_HELPER*sunhrs_day;
			
			if ( drand48()<IData.frost_[month_in_year] ) frost = 1;
			else frost = 0;
			
			MyTreePop.RunDeathProcess(frost);
			
			// compute reference evapotranspiration for soil decomposition (Yasso) ----------------------------------------
			EtSiteRef = FOAPenManRef( IData.tmp_[month_in_year], IData.tmp_[(month_in_year-1)%12],
									  IData.s12_[month_in_year], radnet_day, IData.gama12_[month_in_year],
// 									  WindAtHeight(2.,IData.wnd_[month_in_year]),
									  WindAtHeight(MyTreePop.getMeanHeight(),IData.wnd_[month_in_year]),
									  IData.eS12_[month_in_year], IData.eA12_[month_in_year] );
			
			
			// Evapotranspiration of soil ---------------------------------------------------------------------------------
			if ( MyGrassPop.getAllDormant()<=3  )   // no respiration if 1/2 of grasses are active
				EtSiteGround = 0;
			else {
				double gb_soil = GetgbCANOPY( IData.wnd_[month_in_year], 2.);
				EtSiteGround   = 0.0864/2.45*exp(-4.28+11.97*MyMin(0.35,IData.theta_[0]))*IData.VPD12_[month_in_year]/gb_soil;
			}
			
			
			
			// Evapotranspiration of grasses --------------------------------------------------------------------------------
			EtSiteGrass = MyGrassPop.getEt(  IData.tmp_day_[month_in_year],
						  IData.atm_press_, radnet_day, IData.s12_[month_in_year],
						  SP_HEAT, IData.gama12_[month_in_year],
						  IData.rho12_[month_in_year], IData.VPD12_[month_in_year],
						  MyTreePop.getpCanopy() );  //mm/day
			
			// Evapotranspiration of trees --------------------------------------------------------------------------------
			EtSiteTrees = MyTreePop.getEt( IData.tmp_day_[month_in_year],
						  IData.atm_press_, radnet_day, IData.s12_[month_in_year],
						  SP_HEAT, IData.gama12_[month_in_year],
						  IData.rho12_[month_in_year], IData.VPD12_[month_in_year] );  //mm/day
			
			//-----------------------------------------------------------------------------------------
// 			cout << "ETR " << EtSiteGrass << " " << EtSiteTrees << endl;
			// Evapotranspiration total
			EtSite = EtSiteGround + EtSiteTrees + EtSiteGrass;
			
			
			rain_sum_year        += Rain[day_in_year];
			evapo_sum_year       += EtSite;
			rain_sum             += Rain[day_in_year];
			evapo_sum            += EtSite;
			evapo_grass_sum      += EtSiteGrass;
			evapo_soil_sum       += EtSiteGround;
			evapo_ref_sum        += EtSiteRef;
			
			if ( Rain[day_in_year]>EtSite )
			{
				BucketIn(  Rain[day_in_year]-EtSite, IData.theta_, IData.theta_fc_, IData.thickness_, IData.soil_layers_ );
			}
			if ( Rain[day_in_year]<=EtSite )
			{
				BucketOut( EtSite-Rain[day_in_year], IData.theta_, IData.theta_wp_, IData.thickness_, IData.soil_layers_ );
			}
			
			
			GetSoilTheta( diff_fc_wp, IData.theta_, IData.theta_wp_, IData.g_theta_, IData.soil_layers_ );
			
			// Run tree and grass physiology.
			MyTreePop.RunPhysiology(IData.A012C3_[month_in_year], IData.RmL12C3_[month_in_year], mmsTOkgd,
									IData.tmp_day_[month_in_year], IData.wnd_[month_in_year],
									IData.g_theta_, ca_par_preassure, IData.atm_press_,
// 									IData.reh_[month_in_year], MyGrassPop.getHeightOpn(),       MyGrassPop.getHeightOpn(),
									IData.reh_[month_in_year], MyGrassPop.getHeight(GR_C4_OPN), MyGrassPop.getHeight(GR_C3_OPN),
									day_in_year, month_in_year, IData.theta_[0], IData.theta_fc_[1], IData.theta_wp_[1],
									T_fac, frost, MyGrassPop.getC34Ratio(), IData.thickness_, IData.soil_layers_  );
			
			
#ifndef S_NOGRASS
			MyGrassPop.RunPhysiology( IData.A012C4_[month_in_year],  IData.A012C3_[month_in_year],
									  IData.RmL12C4_[month_in_year], IData.RmL12C3_[month_in_year], mmsTOkgd,
									  IData.tmp_day_[month_in_year], IData.wnd_[month_in_year],
									  IData.g_theta_, ca_par_preassure, IData.atm_press_,
									  IData.reh_[month_in_year], MyTreePop.getMeanSavHeight(),
									  MyTreePop.getMeanForHeight(), T_fac, frost, MyTreePop.getpCanopySav(),
									  MyTreePop.getpCanopyFor(), IData.thickness_, day_in_year );
									  
									  
		 void RunPhysiology( double p_can_sav, double p_can_for ); //NOTE what is that???
														  
			#endif
			
			
			
			
// 			cout << setw(14) << MyGrassPop.getFuelMoisture()*IData.reh_[month_in_year]/100. << setw(14) << MyGrassPop.getFuelMoisture()*(IData.g_theta_[0]+IData.g_theta_[1]+IData.g_theta_[2]+IData.g_theta_[3])/4. << setw(14) << MyGrassPop.getFuelMoisture()*(IData.pwet_[month_in_year]+IData.reh_[month_in_year]/100.) << setw(14) << MyGrassPop.getFuelMoisture() << endl;
			
#			ifdef S_FIRE_FROM_FILE
			if ( count_years==FIRE_year && day_in_year==FIRE_day )
#           else
			if ( Ignitions[day_in_year]>0 )
#           endif
			{
				double fire_intensity;
				
				double dead_fuel = MyGrassPop.getDryBiomassForFire();   // in kg/m^2
								 // + MyTreePop.getDryBiomassForFire();   // in kg/m^2
				
				double live_fuel  = MyGrassPop.getWetBiomassForFire();
				
				double live_fuel_moisture = IData.reh_[month_in_year]/100.+IData.pwet_[month_in_year];
				
				double dead_fuel_moisture = MyGrassPop.getFuelMoisture()*live_fuel_moisture;
				
				
#				ifdef S_FIRE_FROM_FILE
				fire_intensity = FIRE_intensity;
#				else
				fire_intensity = LightFire( dead_fuel, live_fuel, dead_fuel_moisture, live_fuel_moisture,
											IData.wnd_[month_in_year], day_in_year );
#				endif
				
				if ( fire_intensity>0 )
				{
					fire_flag = 1;
					double patchiness    = Patchiness( fire_intensity );
					double scorch        = ScorchHeight( fire_intensity );
					double cc_fine       = CombComplFine( scorch );
					double cc_coarse     = CombComplCoarse( scorch );
					double cc_heavy      = CombComplHeavy( scorch );
					double cc_tk_helper  = CombComplTopkillHelper( scorch );
					
					MyTreePop.setStateAfterFire ( fire_intensity, patchiness, cc_fine, cc_coarse, cc_heavy, cc_tk_helper );
					MyGrassPop.setStateAfterFire( patchiness, cc_fine );
					
					fire_num++;
					total_fire_intensity += fire_intensity;
					
					#ifdef W_SYSDATA
					#include "fire_output_variables.h"
					#endif
					
					DEBUG( DEBUG_FIRE_DATA, "FIRE_LGD "
							<< setw(4)  << "fn"
							<< setw(4)  << "yr"
							<< setw(4)  << "day"
							<< setw(13) << "dead_f"
							<< setw(13) << "live_f"
							<< setw(13) << "dead_mois"
							<< setw(13) << "live_mois"
							<< setw(13) << "fuel"
							<< setw(13) << "moisture"
							<< setw(13) << "intensity"
							<< setw(13) << "patchiness"
							<< setw(13) << "scorch"
							<< setw(13) << "cc_fine"
							<< setw(13) << "cc_coarse"
							<< setw(13) << "cc_heavy"
							<< setw(13) << "cc_tk_helper" )

					DEBUG( DEBUG_FIRE_DATA, "FIRE_DBG "
							<< setw(4)  << fire_num
							<< setw(4)  << count_years
							<< setw(4)  << day_in_year
							<< setw(13) << dead_fuel
							<< setw(13) << live_fuel
							<< setw(13) << dead_fuel_moisture
							<< setw(13) << live_fuel_moisture
							<< setw(13) << dead_fuel+live_fuel
							<< setw(13) << (live_fuel*live_fuel_moisture + dead_fuel*dead_fuel_moisture)/(live_fuel+dead_fuel)
							<< setw(13) << fire_intensity
							<< setw(13) << patchiness
							<< setw(13) << scorch
							<< setw(13) << cc_fine
							<< setw(13) << cc_coarse
							<< setw(13) << cc_heavy
							<< setw(13) << cc_tk_helper
						 )
					
#					ifdef S_FIRE_FROM_FILE
					FIRE_file >> FIRE_year >> FIRE_day >> FIRE_intensity;
#					endif
					
					total_nee -= ( MyGrassPop.getLeafLiveCombustion()
								 + MyGrassPop.getLeafDeadStCombustion()
								 + MyGrassPop.getLeafDeadLyCombustion()
								 + MyTreePop.getStemLiveCombustion()
								 + MyTreePop.getLeafLiveCombustion()
								 + MyTreePop.getStemDeadStCombustion()
								 + MyTreePop.getStemDeadLyCombustion()
								 + MyTreePop.getLeafDeadStCombustion()
								 + MyTreePop.getLeafDeadLyCombustion() )*0.44;   // in kg C/m^2
				
				}  // if ( fire_intensity>0 )
				
			}  // if ( Ignitions[day_in_year]>0 )
			
			
			
#           ifdef S_ELEPHANTS

            if ( count_years>150 )
//             if ( count_years>100 )
            {

                myElephants.GenerateVisitSequence( day_in_year, count_years );

                if ( myElephants.getVisitSequence(day_in_year) > 0 )   // elephant impact, not daily
                {
                    myElephants.getGrassFraction( IData.g_theta_[0] );

                    double biomass_left = MyTreePop.setStateAfterElephants(
                        myElephants.getTreeConsumption( day_in_year ),
                        myElephants.getMaleProb(),
                        myElephants.getTreeMortality(),
                        myElephants.getMaxTreeCons() );

                    biomass_left = MyGrassPop.setStateAfterElephants(
                        myElephants.getGrassConsumption( day_in_year )+biomass_left,
                        MyTreePop.getpCanopy() );
                }
            }
            // daily grass biomass removal by grazers
            MyGrassPop.setStateAfterElephants( myElephants.getGrazingBm(), MyTreePop.getpCanopy() );

            if ( day_in_year == 0 )
            {
                int ELFS_small_trees = 0;
                int ELFS_large_trees = 0;
            
                for ( int i=0; i<MyTreePop.getPopSize(); i++ )
                {
                    if( MyTreePop.getTreeHeight(i) < 1. ) ELFS_small_trees++;
                    if( MyTreePop.getTreeHeight(i) > 5. ) ELFS_large_trees++;
                }
            
                _ELFS_Data  << setw(14) << d_params[0]   // 1
                            << setw(14) << d_params[1]   // 2
                            << setw(14) << d_params[2]   // 3
                            << setw(14) << d_params[3]   // 4
                            << setw(14) << d_params[4]   // 5
                            << setw(14) << d_params[5]   // 6
                            << setw(14) << d_params[6]   // 7
                            << setw(14) << d_params[7]   // 8
                            << setw(14) << d_params[8]   // 9
                            << setw(14) << d_params[9]   // 10
                            << setw(14) << MyTreePop.getBlYearMax()*10.             // 11  in t/ha
                            << setw(14) << MyTreePop.getBsYearMax()*10.             // 12  in t/ha
                            << setw(14) << MyTreePop.getBrYearMax()*10.             // 13  in t/ha
                            << setw(14) << MyTreePop.getMaxBasalAreaYearMean()      // 14
                            << setw(14) << MyTreePop.getPopSize()                   // 15
                            << setw(14) << ELFS_small_trees                         // 16
                            << setw(14) << ELFS_large_trees                         // 17
//                             << setw(14) << mean_gw/MyTreePop.getPopSize()           // 18
//                             << setw(14) << mean_ci/MyTreePop.getPopSize()           // 19
//                             << setw(14) << mean_qi/MyTreePop.getPopSize()           // 20
//                             << setw(14) << mean_qs/MyTreePop.getPopSize()           // 21
//                             << setw(14) << mean_acs/MyTreePop.getPopSize()          // 22
                            << endl;
            }
#           endif   // S_ELEPHANTS


			// Update soil pools, biomasses given in kg/m^2
			mySoil.UpdateCarbonPools(	MyTreePop.getSoilNWL()+		// non woody litter    - tree fine roots and dead leaf
										MyGrassPop.getSoilNWL(),	// non woody litter    - grass roots and dead leaf
										MyTreePop.getSoilFWL(),		// fine woody litter   - tree coarse roots and fine stem
										MyTreePop.getSoilCWL(),		// coarse woody litter - tree stem coarse and heavy
										tmp_mean,
										(1.-MyTreePop.getDormant(0))*(Rain[day_in_year]-EtSiteRef )  );
			
// 			soil_carbon = mySoil.GetCarbonStored()*1000.;  // trasformed to g/m^2, needed for photosynthesis


#           ifdef W_SYSDATA
            if ( day_in_year%TIMESTEP_DAILY_DATA==0 )
                _SOILDATA_SoilData
                    << setw(14) << mySoil.GetXfwl()*10.     // in t/ha
                    << setw(14) << mySoil.GetXcwl()*10.     // in t/ha
                    << setw(14) << mySoil.GetXext()*10.     // in t/ha
                    << setw(14) << mySoil.GetXcel()*10.     // in t/ha
                    << setw(14) << mySoil.GetXlig()*10.     // in t/ha
                    << setw(14) << mySoil.GetXhu1()*10.     // in t/ha
                    << setw(14) << mySoil.GetXhu2()*10.     // in t/ha
                    << setw(14) << mySoil.GetUnwl()*10.     // in t/ha
                    << setw(14) << mySoil.GetUfwl()*10.     // in t/ha
                    << setw(14) << mySoil.GetUcwl()*10.     // in t/ha
                    << setw(14) << mySoil.GetRext()*10.     // in t/ha
                    << setw(14) << mySoil.GetRcel()*10.     // in t/ha
                    << setw(14) << mySoil.GetRlig()*10.     // in t/ha
                    << setw(14) << mySoil.GetRhu1()*10.     // in t/ha
                    << setw(14) << mySoil.GetRhu2()*10.     // in t/ha
                    << endl;
#           endif
			
			total_npp +=  (( MyGrassPop.getGPP() - MyGrassPop.getRma() - MyGrassPop.getRgr()
						    + MyTreePop.getGPP()  - MyTreePop.getRma()  - MyTreePop.getRgr() )*0.44 ); // in kg C/m^2
			total_nee += ((( MyGrassPop.getGPP() - MyGrassPop.getRma() - MyGrassPop.getRgr()
						    + MyTreePop.getGPP()  - MyTreePop.getRma()  - MyTreePop.getRgr() )*0.44 )
						    - mySoil.GetCarbonRelease() ); // in kg C/m^2;
							
							
// ==========================================================================================================
//
//      IT FOLLOWS ONLY DATA OUTPUT AND THE ENDS OF THE DAY LOOP AND THE YEAR LOOP
//
// ==========================================================================================================
			
			
			
#       ifdef W_SYSDATA
#       include "daily_output_variables.h"
#       endif

// ---------------------------------------------------------------------------------
// 			write input data for aDGVM.2
			#ifdef GET_INPUT_DATA
			cout << setw(14) << IData.Q012_[(int)floor((double)day_in_year/30.4)]/10.
				 << setw(14) << radnet_day
				 << setw(14) << IData.tmp_day_[(int)floor((double)day_in_year/30.4)]
				 << setw(14) << IData.reh_[(int)floor((double)day_in_year/30.4)]
				 << setw(14) << IData.wnd_[(int)floor((double)day_in_year/30.4)]
				 << setw(14) << IData.atm_press_
				 << setw(14) << ca_par_preassure*10.
				 << setw(14) << Rain[day_in_year]
				 << setw(14) << IData.theta_wp_[0]
				 << setw(14) << IData.theta_fc_[0]
				 << setw(14) << IData.sun_[(int)floor((double)day_in_year/30.4)]
				 << setw(14) << IData.soil_C_
				 << setw(14) << IData.soil_N_
				 << endl;
			#endif
// ---------------------------------------------------------------------------------
		
		} // End of day loop
		
		if ( count_years>100 )
			total_basal_area += MyTreePop.getMaxBasalAreaYearMean();
		
// 		cout << "GSL" << setw(14) << MyTreePop.getActiveDays() << setw(14) << MyGrassPop.getActiveDays() << setw(14) << rain_sum_year << endl;
		
		int ELFS_small_trees = 0;
		int ELFS_large_trees = 0;

		#ifdef S_ELEPHANTS
		for ( int i=0; i<MyTreePop.getPopSize(); i++ )
		{
			if( MyTreePop.getTreeHeight(i) < 1. ) ELFS_small_trees++;
			if( MyTreePop.getTreeHeight(i) > 5. ) ELFS_large_trees++;
		}
		#endif
		
		DEBUG(  DEBUG_YEARLY_DATA, "YEAR_DBG " <<
						setw( 5) << count_years <<
						setw(12) << MyGrassPop.getLeafBmLive()*10. <<  // t/ha
						setw(12) << MyGrassPop.getRootBmLive()*10. <<  // t/ha
						setw(12) << MyTreePop.getStemBmLive()*10.+MyTreePop.getLeafBmLive()*10.  <<   // t/ha
						setw(12) << MyTreePop.getRootBmLive()*10. <<  // t/ha
						setw( 6) << MyTreePop.getPopSize() <<
						setw( 6) << MyTreePop.getNewBorn() <<
						setw( 6) << MyTreePop.getDeadTrees() <<
						setw(12) << MyTreePop.getpCanopy() <<
						setw(12) << MyTreePop.getMaxHeight() <<
						setw(12) << MyTreePop.getMeanHeight() <<
						setw(12) << MyTreePop.getBasalArea() <<
						setw(12) << rain_sum_year <<
						setw( 6) << MyTreePop.getSavTreeNum() <<
						setw( 6) << MyTreePop.getForTreeNum() <<
						setw(12) << MyTreePop.getpCanopySav() <<
						setw(12) << MyTreePop.getpCanopyFor() <<
						setw(12) << MyGrassPop.getC34Ratio() 
		)
		
		
		
#       ifdef W_EBP      // data output in 1954 and 1997 (EBP experiments) 
        if ( count_years==193 || count_years==235 )
        {
#	       include "yearly_output_variables.h"    // file defines variables for yearly output
            myYearlyData.setYearlyData( yearly_output_data );
            myYearlyData.printYearlyDataToFile( output_directory );
        }
#       endif


		#ifdef S_CLIM_SENSITIVITY  // allows to set CO2, temperature and rainfall
		if ( count_years%50==0 )
		{
			#include "yearly_output_variables.h"    // file defines variables for yearly output
			myYearlyData.setYearlyData( yearly_output_data );
			myYearlyData.printYearlyDataToFile( output_directory );
		}
		#endif



		#ifdef S_CLIM_NO_C4
		if ( count_years==150 || count_years==200 || count_years==250  || count_years==300 || count_years==450 || count_years==600 )
		{
			#include "yearly_output_variables.h"    // file defines variables for yearly output
			myYearlyData.setYearlyData( yearly_output_data );
			myYearlyData.printYearlyDataToFile( output_directory );
		}
		#endif
		
		
		#ifdef S_CLIM_C4_FROM_150
		if ( count_years==150 || count_years==200 || count_years==250  || count_years==300 || count_years==450 || count_years==600 )
		{
			#include "yearly_output_variables.h"    // file defines variables for yearly output
			myYearlyData.setYearlyData( yearly_output_data );
			myYearlyData.printYearlyDataToFile( output_directory );
		}
		#endif
		
		#ifdef S_CLIM_C4_FROM_600
		if ( count_years==150 || count_years==200 || count_years==250  || count_years==300 || count_years==450 || count_years==600 )
		{
			#include "yearly_output_variables.h"    // file defines variables for yearly output
			myYearlyData.setYearlyData( yearly_output_data );
			myYearlyData.printYearlyDataToFile( output_directory );
		}
		#endif
		
//		#ifdef CLIM_SCEN
// 		output after 100 years (spin-up) and in 2011
//		if ( count_years==years_to_run )
//        {
//			#include "yearly_output_variables.h"    // file defines variables for yearly output
//            myYearlyData.setYearlyData( yearly_output_data );
//            myYearlyData.printYearlyDataToFile( output_directory );
//        }
//		#endif

// 		output at end of simulation
		if ( count_years==years_to_run )
		{
#   		include "yearly_output_variables.h"    // file defines variables for yearly output
			myYearlyData.setYearlyData( yearly_output_data );
			myYearlyData.printYearlyDataToFile( output_directory );
		}
		
		
		#ifdef S_FIXED_RAND_INIT
		if ( count_years==100) srand48(1122554162+(int)d_seed);
		#endif
		
		

    } // end of year loop




#   ifdef W_SYSDATA
	_SYSDATA_SysData.close();
	_SIZE_STRUC_SizeStrucFile.close();
	_FIREDATA_FireData.close();
	_SOILDATA_SoilData.close();
#   endif

#   ifdef S_ELEPHANTS
    _ELFS_Data.close();
#   endif

#	ifdef S_RAIN_FROM_FILE
	FILE_prec_file.close();
#	endif
	
#   ifdef S_FIRE_FROM_FILE
    FIRE_file.close();
#   endif

#   ifdef CLIM_SCEN
	cs_co2_file.close();
	cs_tmp_file.close();
#	endif
	
	delete[] diff_fc_wp;
	IData.deleteMemory();
	
	
// 	cout << longitude << "  " << latitude << " " << radnetsum/365. << endl;
	
	
	return 0;
}














