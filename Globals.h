#ifndef Globals_h___
#define Globals_h___

int const NUM_TR_TYPES  = 3;
int const NUM_GR_TYPES  = 8;

int const TR_SAV        = 0;
int const TR_FOR        = 1;
int const TR_BBS        = 2;

int const GR_C4OPN      = 0;
int const GR_C3OPN      = 1;

int const GR_C4SAV      = 2;
int const GR_C3SAV      = 3;

int const GR_C4FOR      = 4;
int const GR_C3FOR      = 5;

int const GR_C4BBS      = 6;
int const GR_C3BBS      = 7;


#ifndef OPTIM_GLOBALS
const double DEATH_PROB_FROST[NUM_TR_TYPES]  = { 0.9999, 0.2848697, 1e-04 };
const double DEATH_PROB_CARBON[NUM_TR_TYPES] = { 0.3644492, 0.2074071, 1e-04 };
const double DEATH_PROB_COMP[NUM_TR_TYPES]   = { 1e-04, 0.7863725, 0.9999 };
const double IGNITION_PROB  = 1e-04;
const double IGNITION_PAR_2 = 0.06113167;
const double L_SAV_SAV = 0.835574;
const double L_SAV_FOR = 1e-04;
const double L_SAV_BBS = 0.603885;
const double L_FOR_SAV = 0.6586143;
const double L_FOR_FOR = 1e-04;
const double L_FOR_BBS = 0.7283344;
const double L_BBS_SAV = 0.3798453;
const double L_BBS_FOR = 0.9999;
const double L_BBS_BBS = 0.02925457;
const double L_C4G_SAV = 0.9999;
const double L_C4G_FOR = 0.5847419;
const double L_C4G_BBS = 0.9999;
const double L_C3G_SAV = 1e-04;
const double L_C3G_FOR = 0.1586156;
const double L_C3G_BBS = 1e-04;
const double L_C4G_C4G = 0.2287690;
const double L_C3G_C3G = 0.4015056;
const double L_SAV_C4G = 0.1337562;
const double L_SAV_C3G = 0.4083678;
const double L_FOR_C4G = 0.1232729;
const double L_FOR_C3G = 0.4748842;
const double L_BBS_C4G = 0.694666;
const double L_BBS_C3G = 1e-04;

// NOTE these are the "default" values
// const double DEATH_PROB_FROST[NUM_TR_TYPES]		= { 0.001,  0.001, 0.001 };
// const double DEATH_PROB_CARBON[NUM_TR_TYPES]	= { 0.001, 0.001, 0.0001 };
// const double DEATH_PROB_COMP[NUM_TR_TYPES]		= { 0.001, 0.001, 0.0001 };
// const double DEATH_PROB_DROUGHT[NUM_TR_TYPES]	= { 0.0000,  0.0000, 0.00 };
// 
// const double DEATH_PROB_FROST[NUM_TR_TYPES]		= { 0.001,  0.001, 0.00001 };
// // const double DEATH_PROB_CARBON[NUM_TR_TYPES]	= { 0.0005, 0.001, 0.25 };
// const double DEATH_PROB_CARBON[NUM_TR_TYPES]	= { 0.001, 0.001, 0.001 };
// const double DEATH_PROB_COMP[NUM_TR_TYPES]		= { 0.001, 0.001, 0.0001 };
// const double DEATH_PROB_DROUGHT[NUM_TR_TYPES]	= { 0.00001,  0.00001, 0.001 };
// const double IGNITION_PROB					= 0.01;				// probability of fire ignition
// const double IGNITION_PAR_2					= 0.1;				// parameter describing fire ignition sequence, treecover
// 
// // trees shade trees
// const double L_SAV_SAV = 0.1;
// const double L_SAV_FOR = 0.1;
// const double L_SAV_BBS = 0.1;
// 
// const double L_FOR_SAV = 0.3;
// const double L_FOR_FOR = 0.3;
// const double L_FOR_BBS = 0.3;
// 
// const double L_BBS_SAV = 0.99;
// const double L_BBS_FOR = 0.99;
// const double L_BBS_BBS = 0.99;
// 
// // grasses shade trees
// const double L_C4G_SAV = 0.1;
// const double L_C4G_FOR = 0.1;
// const double L_C4G_BBS = 0.1;
// 
// const double L_C3G_SAV = 0.1;
// const double L_C3G_FOR = 0.1;
// const double L_C3G_BBS = 0.1;
// 
// // grasses shade grasses
// const double L_C4G_C4G = 0.2;
// const double L_C3G_C3G = 0.2;
// 
// // trees shade grasses
// const double L_SAV_C4G = 0.2;
// const double L_SAV_C3G = 0.2;
// 
// const double L_FOR_C4G = 0.5;
// const double L_FOR_C3G = 0.5;
// 
// const double L_BBS_C4G = 0.999;
// const double L_BBS_C3G = 0.999;



const double LC_TR_TR_1[NUM_TR_TYPES][NUM_TR_TYPES] = {
	L_SAV_SAV,    L_SAV_FOR,    L_SAV_BBS,
	L_FOR_SAV,    L_FOR_FOR,    L_FOR_BBS,
	L_BBS_SAV,    L_BBS_FOR,    L_BBS_BBS };
	
const double LC_TR_TR_2[NUM_TR_TYPES][NUM_TR_TYPES] = { 
	1.-L_SAV_SAV, 1.-L_SAV_FOR, 1.-L_SAV_BBS,
	1.-L_FOR_SAV, 1.-L_FOR_FOR, 1.-L_FOR_BBS,
	1.-L_BBS_SAV, 1.-L_BBS_FOR, 1.-L_BBS_BBS };

const double LC_C3_TR_1[NUM_TR_TYPES] = {    L_C3G_SAV,    L_C3G_FOR,    L_C3G_BBS };
const double LC_C4_TR_1[NUM_TR_TYPES] = {    L_C4G_SAV,    L_C4G_FOR,    L_C4G_BBS };

const double LC_C3_TR_2[NUM_TR_TYPES] = { 1.-L_C3G_SAV, 1.-L_C3G_FOR, 1.-L_C3G_BBS };
const double LC_C4_TR_2[NUM_TR_TYPES] = { 1.-L_C4G_SAV, 1.-L_C4G_FOR, 1.-L_C4G_BBS };
	
const double LC_GR_GR_1[NUM_GR_TYPES] = {    L_C4G_C4G,    L_C3G_C3G,
	L_C4G_C4G,    L_C3G_C3G,
	L_C4G_C4G,    L_C3G_C3G,
	L_C4G_C4G,    L_C3G_C3G };
const double LC_GR_GR_2[NUM_GR_TYPES] = { 1.-L_C4G_C4G, 1.-L_C3G_C3G,
	1.-L_C4G_C4G, 1.-L_C3G_C3G,
	1.-L_C4G_C4G, 1.-L_C3G_C3G,
	1.-L_C4G_C4G, 1.-L_C3G_C3G };
	
const double LC_TR_GR_1[NUM_TR_TYPES][NUM_GR_TYPES] = {
	L_SAV_C4G,    L_SAV_C3G,    L_SAV_C4G,    L_SAV_C3G,    L_SAV_C4G,    L_SAV_C3G,    L_SAV_C4G,    L_SAV_C3G,
	L_FOR_C4G,    L_FOR_C3G,    L_FOR_C4G,    L_FOR_C3G,    L_FOR_C4G,    L_FOR_C3G,    L_FOR_C4G,    L_FOR_C3G,
	L_BBS_C4G,    L_BBS_C3G,    L_BBS_C4G,    L_BBS_C3G,    L_BBS_C4G,    L_BBS_C3G,    L_BBS_C4G,    L_BBS_C3G };
	
const double LC_TR_GR_2[NUM_TR_TYPES][NUM_GR_TYPES] = {
	1.-L_SAV_C4G, 1.-L_SAV_C3G, 1.-L_SAV_C4G, 1.-L_SAV_C3G, 1.-L_SAV_C4G, 1.-L_SAV_C3G, 1.-L_SAV_C4G, 1.-L_SAV_C3G,
	1.-L_FOR_C4G, 1.-L_FOR_C3G, 1.-L_FOR_C4G, 1.-L_FOR_C3G, 1.-L_FOR_C4G, 1.-L_FOR_C3G, 1.-L_FOR_C4G, 1.-L_FOR_C3G,
	1.-L_BBS_C4G, 1.-L_BBS_C3G, 1.-L_BBS_C4G, 1.-L_BBS_C3G, 1.-L_BBS_C4G, 1.-L_BBS_C3G, 1.-L_BBS_C4G, 1.-L_BBS_C3G };

#endif






// Model directories

// Home directory of model, contains the source code and the executable
const string MODEL_HOME					= "";

// Home directory of input data, contains
//    filenames.txt        Contain names of global databases,
//                         these databases must be in the same directory
//    shortlists.txt       Contains the names of shortlists that have already been read
//                         for study sites
//    shortlists           Directory that contains shortlists that have already been read
//    ClimateChange        Directory that contains files with climate change scenarios
//                         (CO2, precipitation, temperature)
//    fire_sites           Directory with fire scenarios for specific study sites, these
//                         scenarios are used when compiler flag S_FIRE_FROM_FILE is set
//    prec_sites           Directory with precipitation data for specific study sites, these
//                         data are used when compiler flag S_RAIN_FROM_FILE is set
const string IN_DATA_HOME				= "/Users/glennmoncrieff/adgvm/data/InputData/";

// Home directory of output data, this directory contains sub-directories (experiment_xyz).
// Output data are written in OUT_DATA_HOME/(experiment_xyz),
// the directory (experiment_xyz) is specified in the model command line.
const string OUT_DATA_HOME				= "/Users/glennmoncrieff/adgvm/data/OutputData/aDGVM121";




// Model constants

// Start constants for the leaf photosynthesis model

const double ABS_PHOTONS_C3				= 0.86;			// absorbtance to incident flux of photons, Collatz C3
const double ALPHA_C3					= 0.08;			// intrinsic quantum efficiency for C02 uptake, Collatz C3
const double ALPHAR_F_C4				= 0.067;		// product of alpha_r (intrinsic quantum yield of C3
const double ABS_PHOTONS_C4				= 0.80;			// absorbtance to incident flux of photons, Collatz C4
														// photosynthesis) and f (fraction of absorbed photons used
														// by C3 reactions), FJiC4, Collatz C4 units mol/mol
const double KAPPA_C4					= 0.7*1e6;		// initial slope of photosynthetic CO2 resp for C4, 0.7 mumol/m2/s
														// FC4Jc, Collatz C4



#ifdef ELEPHANTS_FUT
const double CA_PAR_PREASSURE           = 70.09;        // Atmospheric partial pressure C02
#elif defined PLUS_CO2
const double CA_PAR_PREASSURE           = 70.09;        // Atmospheric partial pressure C02

#elif defined S_GCM_180
const double CA_PAR_PREASSURE           = 18.0;        // Atmospheric partial pressure C02
#elif defined S_GCM_280
const double CA_PAR_PREASSURE           = 28.0;        // Atmospheric partial pressure C02
#elif defined S_GCM_400
const double CA_PAR_PREASSURE           = 40.0;        // Atmospheric partial pressure C02

#elif defined S_CO2_170PPM
const double CA_PAR_PREASSURE           = 17.0;        // Atmospheric partial pressure C02
#elif defined S_CO2_200PPM
const double CA_PAR_PREASSURE           = 20.0;        // Atmospheric partial pressure C02
#elif defined S_CO2_400PPM
const double CA_PAR_PREASSURE           = 40.0;        // Atmospheric partial pressure C02
#elif defined S_CO2_700PPM
const double CA_PAR_PREASSURE           = 70.0;        // Atmospheric partial pressure C02

#else
const double CA_PAR_PREASSURE			= 38.70;		// Atmospheric partial pressure C02
#endif















const double OI_PAR_PREASSURE			= 21.0;			// O2 Partial pressure KPa, Collatz

// const double MTR = 9.;

const double B_TREE[NUM_TR_TYPES]       = { 0.01, 0.01, 0.01 };   // ball berry constant conductance units
const double M_TREE[NUM_TR_TYPES]       = { 9.,   9.,   9.   };   // ball berry constant dimensionless

// const double B_GRASS[NUM_GR_TYPES]      = { 0.04, 0.01, 0.04, 0.01, 0.04, 0.01, 0.04, 0.01 }; // ball berry constant conductance units
// const double M_GRASS[NUM_GR_TYPES]      = { 4.,   9.,   4.,   9.,   4.,   9.,   4.,   9.   }; // ball berry constant dimensionless
const double B_GRASS[NUM_GR_TYPES]      = { 0.04, 0.01, 0.04, 0.01, 0.04, 0.01, 0.04, 0.01 }; // ball berry constant conductance units
const double M_GRASS[NUM_GR_TYPES]      = { 5.,   9.,   5.,   9.,   5.,   9.,   5.,   9.   }; // ball berry constant dimensionless

const double KC_K25_C3					= 30.0; 	   	// c3 value of Kc at 25degC, Pa, Collatz C3
const double KC_Q10_C3					= 2.1;   		// c3 rate of change of Kc with temp, Collatz C3
const double KO_K25_C3					= 30.0;    		// c3 value of Ko at 25degC, KPa, Collatz C3
const double KO_Q10_C3					= 1.2;   		// c3 rate of change of Ko with temp, Collatz C3
const double TAU_K25_C3					= 2600.0;  		// c3 value of tau at 25degC,  Collatz C3
const double TAU_Q10_C3					= 0.57;  		// c3 rate of change of tau with temp, Collatz C3

const double KC_K25_C4					= 140.0;	   	// c4 value of Kc at 25degC, Pa, Collatz C4
const double KC_Q10_C4					= 2.1;   		// c4 rate of change of Kc with temp, Chen
const double KO_K25_C4					= 34.0; 	   	// c4 value of Ko at 25degC, KPa, Collatz C4
const double KO_Q10_C4					= 1.2;   		// c4 rate of change of Ko with temp, Chen
const double TAU_K25_C4					= 2600.0;  		// c4 value of tau at 25degC,  Collatz C4
const double TAU_Q10_C4					= 0.67;  		// c4 rate of change of tau with temp, Chen

// end constants for the leaf photosynthesis model


// start constants needed for leaf boundary layer conductance

const double CLD_TREE					= 0.02;			// charac leaf dimension tree (assumes 2 cm leaf width)
const double CLD_GRASS					= 0.005;  		// charac leaf dimesnion grass (assumes 5 mm leaf width)

// end constants needed for leaf boundary layer conductance


// start constants needed for (canopy) boundary layer conductance
const double Z_WIND						= 10.;			// height at which reference windspeed is measured, New et al. 2000
const double ZD_CONST					= 0.86;			// constant for estimating displacement height, Jones 1992
const double Z0_CONST					= 0.06;			// constant for estimating roughness length, Jones 1992
const double VEGETATION_HEIGHT			= 1.5;			// H_bar average vegetation height in m (needed to calculate wind)
const double ROUGHNESS_LENGTH			= Z0_CONST*VEGETATION_HEIGHT;
const double DISPLACEMENT_HEIGHT		= ZD_CONST*VEGETATION_HEIGHT;
const double REF_HEIGHT_Z				= 10.;			// Height, at which wind is meassured
const double WIND_AT_HEIGHT_HELPER      = 1./log( (REF_HEIGHT_Z-DISPLACEMENT_HEIGHT)/ROUGHNESS_LENGTH );

const double KARMAN_CONST				= 0.41;			// von Karman constant, for canopy boundary layer conductance
const double KARMAN_CONST_QUAD			= pow(KARMAN_CONST,2.);	// quadrat of von Karman constant

// Start CANOPY constants
const double K_CAN_EXT_TREE[NUM_TR_TYPES]			= { 0.5, 0.5, 0.4 };
const double K_CAN_EXT_TREE_INVERSE[NUM_TR_TYPES]	= { 1./K_CAN_EXT_TREE[0], 1./K_CAN_EXT_TREE[1], 1./K_CAN_EXT_TREE[2] };

const double K_CAN_EXT_GRASS[NUM_GR_TYPES]         = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 }; // canopy extinction coefficient
const double K_CAN_EXT_GRASS_INVERSE[NUM_GR_TYPES] = { 1./K_CAN_EXT_GRASS[0], 1./K_CAN_EXT_GRASS[1], 1./K_CAN_EXT_GRASS[2], 
                                                       1./K_CAN_EXT_GRASS[3], 1./K_CAN_EXT_GRASS[4], 1./K_CAN_EXT_GRASS[5], 
                                                       1./K_CAN_EXT_GRASS[6], 1./K_CAN_EXT_GRASS[7] };



// Start respiration constants
const double R_MAINT_RESP_TR[NUM_TR_TYPES]  = { 0.015, 0.015, 0.015 }; // tree RmL leaf  maint respiration as proportion of Vm, value from Collatz  original
const double R_MAINT_RESP_GR[NUM_GR_TYPES]  = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025 }; // grass RmL leaf  maint respiration as proportion of Vm

// Start allometric constants, Table 5
// const double HEIGHT_C1_TREE[NUM_TR_TYPES]   = { 1.295720, 1.295720, 3. }; //1.6 };    // parameters for height(StemBiomass), Higgins2007
// const double HEIGHT_C2_TREE[NUM_TR_TYPES]   = { 0.392157, 0.392157, 0.1 }; //0.25 };    // parameters for height(StemBiomass), Higgins2007

// const double HEIGHT_C1_GRASS[NUM_GR_TYPES]  = { 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5 };  // parameters for height(StemBiomass), Arora2005
// const double HEIGHT_C2_GRASS[NUM_GR_TYPES]  = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };  // parameters for height(StemBiomass), Arora2005

// parameters for alternative allometric equation
const double HEIGHT_C1_TREE[NUM_TR_TYPES]   = { 0.6604492, 0.1604492, 0.9604492 };   // parameters for height(StemBiomass), Higgins2007
const double HEIGHT_C2_TREE[NUM_TR_TYPES]   = { 2.5499324, 2.3499324, 4.9499324 };   // parameters for height(StemBiomass), Higgins2007

const double HEIGHT_C1_GRASS[NUM_GR_TYPES]  = { 2.505441, 2.505441, 2.505441, 2.505441, 2.505441, 2.505441, 2.505441, 2.505441 };  // parameters for height
const double HEIGHT_C2_GRASS[NUM_GR_TYPES]  = { 1.999963, 1.999963, 1.999963, 1.999963, 1.999963, 1.999963, 1.999963, 1.999963 };  // parameters for height



const double STEM_AREA_C1[NUM_TR_TYPES]     = { 2.797, 2.797, 2.797 };        // parameters for basal area, Higgins (2007) Ecology
const double STEM_AREA_C2[NUM_TR_TYPES]     = { 200.,  200.,  200.  };        // parameters for basal area, Higgins (2007) Ecology
const double STEM_AREA_HELPER[NUM_TR_TYPES] = { STEM_AREA_C1[0]/STEM_AREA_C2[0], STEM_AREA_C1[1]/STEM_AREA_C2[1], STEM_AREA_C1[2]/STEM_AREA_C2[2] };

const double GAMMA_CANOPY_TREE[NUM_TR_TYPES]  = { 0.37, 0.42, 0.45 };        // ratio of canopy radius to height for trees
const double GAMMA_CANOPY_GRASS[NUM_GR_TYPES] = { 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4 };  // ratio of canopy radius to height for grasses (estimated)

const double SLA_TREE[NUM_TR_TYPES]           = { 10., 12., 5. };         // specific leaf area tree, m2/kg from ???
const double SLA_GRASS[NUM_GR_TYPES]          = { 10.9, 10.9, 10.9, 10.9, 10.9, 10.9, 10.9, 10.9 }; // specific leaf area grass, m2/kg from Scholes&Walker



const double MAX_ROOT_DEP_TREE			= 2.;			// maximum rooting depth (m) trees -see Schenk and Jackson papers for estimates
const double ROOT_DENSITY_TREE			= 0.1e03;		// density of root biomass (kg/m3) assumes root biomass is a light wood
const double MIN_ROOT_RAD_TREE			= 0.015;		// minimum root radius, m

const double MAX_ROOT_DEP_GRASS			= 0.3;			// maximum rooting depth (m) grasses
const double ROOT_DENSITY_GRASS			= 0.1e03;		// density of root biomass (kg/m3) assumes root biomass is a light wood
const double MIN_ROOT_RAD_GRASS			= 0.005;		// minimum root radius
// end allometric constants



// Start constants describing respiration, Table 7
const double BETA_N						= 0.218;  		// respiration rate for roots and stems kg C per Kg N, Arora et al. 2003 
// (after Keyser et al. 2000)

const double BETA_STEM_TREE				= 0.025;		// fraction of respirating tissue
const double BETA_ROOT_TREE				= 0.025;		// fraction of respirating tissue
const double BETA_ROOT_GRASS			= 0.15;			// fraction of respirating tissue


const double UPSILON_STEM_TREE[NUM_TR_TYPES]       = { 150., 150., 150. };			// C:N for tree stems note this needs to be higher than root value
const double UPSILON_ROOT_TREE[NUM_TR_TYPES]       = {  60.,  60.,  60. };			// C:N for tree roots
const double UPSILON_STEM_GRASS[NUM_GR_TYPES]      = { 120., 120., 120., 120., 120., 120., 120., 120. }; // C:N ratio, C4 have higher C:N than C3
const double UPSILON_ROOT_GRASS[NUM_GR_TYPES]      = { 120., 120., 120., 120., 120., 120., 120., 120. }; // C:N ratio

const double SIGMA_GROW_RESP_TREE[NUM_TR_TYPES]    = { 0.35, 0.35, 0.35 };		// growth respiration constant (Arora =0.35)
const double SIGMA_GROW_RESP_GRASS[NUM_GR_TYPES]   = { 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35 }; // growth respiration constant (Arora =0.35)
// End constants describing respiration




// Start allocation constants, Table 8
const double A0_ROOT_GRASS[NUM_GR_TYPES]        = { 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4 };   //optimal resource allocation proportion to root
const double A0_STEM_GRASS[NUM_GR_TYPES]        = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };   //optimal resource allocation proportion to stem
const double A0_LEAF_GRASS[NUM_GR_TYPES]        = { 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6 };   //optimal resource allocation proportion to leaf

const double A0_ROOT_TREE[NUM_TR_TYPES]         = { 0.35, 0.15, 0.3  };        // optimal resource allocation proportion to root
const double A0_STEM_TREE[NUM_TR_TYPES]         = { 0.35, 0.50, 0.45 };        // optimal resource allocation proportion to stem
const double A0_LEAF_TREE[NUM_TR_TYPES]         = { 0.30, 0.35, 0.25 };        // optimal resource allocation proportion to leaf

const double A0_LEAF_TREE_INV[NUM_TR_TYPES]		= { 1./A0_LEAF_TREE[0], 1./A0_LEAF_TREE[1], 1./A0_LEAF_TREE[2] };
const double A0_ROOT_TREE_HELPER[NUM_TR_TYPES]	= { 1.+A0_ROOT_TREE[0], 1.+A0_ROOT_TREE[1], 1.+A0_ROOT_TREE[2] };
const double A0_STEM_TREE_HELPER[NUM_TR_TYPES]	= { 1.+A0_STEM_TREE[0], 1.+A0_STEM_TREE[1], 1.+A0_STEM_TREE[2] };

const double ALLOC_DENOM_HELPER[NUM_TR_TYPES]	= { 3.+A0_ROOT_TREE[0]+A0_STEM_TREE[0],
                                                    3.+A0_ROOT_TREE[1]+A0_STEM_TREE[1],
                                                    3.+A0_ROOT_TREE[2]+A0_STEM_TREE[2] };
// End allocation constants



// Start phenology constants, Table 9
// const double SLA_TREE[NUM_TR_TYPES]           = { 10., 12., 5. };         // specific leaf area tree, m2/kg from ???
const double LEAF_RETURN[NUM_TR_TYPES]        = { -3., -4., -10. };
// const double LEAF_RETURN[NUM_TR_TYPES]        = { 0., 0., 0. };

const double STRESS_INDEX_TREE[NUM_TR_TYPES]    = { 0., 0., -10. };    // 3. stress index for trees
const double STRESS_INDEX_GRASS[NUM_GR_TYPES]   = { 0., 0., 0., 0., 0., 0., 0., 0. };        // 5. stress index for grass

const double REM_BM_GRASS						= 0.1;			// proportion of grass leaf biomass kept after leaf abscission
const double REM_BM_TREE[NUM_TR_TYPES]			= { 0.01, 0.01, 0.01 };	// proportion of tree leaf biomass kept after leaf abscission

const int    D_POS_TREE[NUM_TR_TYPES]           = { 10, 10, 10 };    // 5 days of positive carbon gain while dormant before growth restarts (days)
const int    D_NEG_TREE[NUM_TR_TYPES]           = {  7,  7,  7 };    // 7 days of negative carbon gain before leaf fall (days)

const int    D_POS_GRASS[NUM_GR_TYPES]          = { 7, 7, 7, 7, 7, 7, 7, 7 };  // days of positive carbon gain while dormant before growth restarts (days)
const int    D_NEG_GRASS[NUM_GR_TYPES]          = { 5, 5, 5, 5, 5, 5, 5, 5 };  // days of negative carbon gain before leaf fall (days)

const double TI_CONST[NUM_TR_TYPES]             = { 15., 15., 15. };
const double TI_CONST_GR[NUM_GR_TYPES]          = { 15., 15., 15., 15., 15., 15., 15., 15. };
// End phenology constants


// Start turnover parameters, Table 10
// trees
const double EXTINC_EXP_TREE				= -0.909;    // exponent for leaf longevity
const double EXTINC_FAC_TREE				= 29664.;    // factor for leaf longevity
const double OMEGA_LEAF_TREE[NUM_TR_TYPES]	= { EXTINC_FAC_TREE*pow(10.*SLA_TREE[0],EXTINC_EXP_TREE),
                                                EXTINC_FAC_TREE*pow(10.*SLA_TREE[1],EXTINC_EXP_TREE),
                                                EXTINC_FAC_TREE*pow(10.*SLA_TREE[2],EXTINC_EXP_TREE)};  // factor 10 to transform units
const double OMEGA_ROOT_TREE				= 1e6;       // root longevity in days 
const double OMEGA_STEM_TREE				= 1e6;       // stem longevity in days 
// const double OMEGA_ROOT_TREE				= 9125.;     // (root longevity in days (Gill and Jackson 2000)
// const double OMEGA_STEM_TREE				= 9125.;     // stem longevity in days (no source)
const double OMEGA_LEAF_TR_TREE				= 30.;       // transition from hanging to lying dead leaf
const double OMEGA_STEM_TR_TREE				= 300.;      // transition from standing to lying dead stem biomass

// grass
const double EXTINC_EXP_GRASS				= -0.909;    // exponent for leaf longevity
const double EXTINC_FAC_GRASS				= 29664.;    // factor for leaf longevity
const double OMEGA_LEAF_GRASS[NUM_GR_TYPES]	= { EXTINC_FAC_GRASS*pow(10.*SLA_GRASS[0],EXTINC_EXP_GRASS),
                                                EXTINC_FAC_GRASS*pow(10.*SLA_GRASS[1],EXTINC_EXP_GRASS),
                                                EXTINC_FAC_GRASS*pow(10.*SLA_GRASS[2],EXTINC_EXP_GRASS),
                                                EXTINC_FAC_GRASS*pow(10.*SLA_GRASS[3],EXTINC_EXP_GRASS),
                                                EXTINC_FAC_GRASS*pow(10.*SLA_GRASS[4],EXTINC_EXP_GRASS),
                                                EXTINC_FAC_GRASS*pow(10.*SLA_GRASS[5],EXTINC_EXP_GRASS),
                                                EXTINC_FAC_GRASS*pow(10.*SLA_GRASS[6],EXTINC_EXP_GRASS),
                                                EXTINC_FAC_GRASS*pow(10.*SLA_GRASS[7],EXTINC_EXP_GRASS) };


const double OMEGA_ROOT_GRASS				= 384.;      // (root longevity in days (Gill and Jackson 2000)
const double OMEGA_STEM_GRASS				= 384.;      // stem longevity in days (no source)
const double OMEGA_LEAF_TR_GRASS			= 150.;     // transition rate from standing to lying dead leaf biomass


// helper used for calculations


const double HELPER_BLD_TREE[NUM_TR_TYPES]	= {    1./OMEGA_LEAF_TREE[0],    1./OMEGA_LEAF_TREE[1],    1./OMEGA_LEAF_TREE[2] };		// for turnover
const double HELPER_BL__TREE[NUM_TR_TYPES]	= { 1.-1./OMEGA_LEAF_TREE[0], 1.-1./OMEGA_LEAF_TREE[1], 1.-1./OMEGA_LEAF_TREE[2] };		// for turnover

	
const double HELPER_BSD_TREE				=    1./OMEGA_STEM_TREE;		// for turnover
const double HELPER_BS__TREE				= 1.-1./OMEGA_STEM_TREE;		// for turnover

const double HELPER_BRD_TREE				=    1./OMEGA_ROOT_TREE;		// for turnover
const double HELPER_BR__TREE				= 1.-1./OMEGA_ROOT_TREE;		// for turnover

const double HELPER_BLTL_TREE				=    1./OMEGA_LEAF_TR_TREE;		// for turnover
const double HELPER_BLT__TREE				= 1.-1./OMEGA_LEAF_TR_TREE;		// for turnover

const double HELPER_BSTL_TREE				=    1./OMEGA_STEM_TR_TREE;		// for turnover
const double HELPER_BST__TREE				= 1.-1./OMEGA_STEM_TR_TREE;		// for turnover

const double HELPER_BLS_GRASS[NUM_GR_TYPES]	= {    1./OMEGA_LEAF_GRASS[0],    1./OMEGA_LEAF_GRASS[1],    1./OMEGA_LEAF_GRASS[2],
                                                   1./OMEGA_LEAF_GRASS[3],    1./OMEGA_LEAF_GRASS[4],    1./OMEGA_LEAF_GRASS[5], 
                                                   1./OMEGA_LEAF_GRASS[6],    1./OMEGA_LEAF_GRASS[7] };    // for turnover
const double HELPER_BL__GRASS[NUM_GR_TYPES]	= { 1.-1./OMEGA_LEAF_GRASS[0], 1.-1./OMEGA_LEAF_GRASS[1], 1.-1./OMEGA_LEAF_GRASS[2],
                                                1.-1./OMEGA_LEAF_GRASS[3], 1.-1./OMEGA_LEAF_GRASS[4], 1.-1./OMEGA_LEAF_GRASS[5],
                                                1.-1./OMEGA_LEAF_GRASS[6], 1.-1./OMEGA_LEAF_GRASS[7] };    // for turnover

const double HELPER_BSS_GRASS				=    1./OMEGA_STEM_GRASS;		// for turnover
const double HELPER_BS__GRASS				= 1.-1./OMEGA_STEM_GRASS;		// for turnover

const double HELPER_BRD_GRASS				=    1./OMEGA_ROOT_GRASS;		// for turnover
const double HELPER_BR__GRASS				= 1.-1./OMEGA_ROOT_GRASS;		// for turnover

const double HELPER_BLTL_GRASS				=    1./OMEGA_LEAF_TR_GRASS;	// for turnover
const double HELPER_BLT__GRASS				= 1.-1./OMEGA_LEAF_TR_GRASS;	// for turnover
// End turnover parameters



// Start constant for transitions from dead (lying) biomass pools into Yasso soil model pools
const double OMEGA_LD_NWL_TREE				= 200.;		// dead leaves to non woody litter
const double OMEGA_RF_NWL_TREE				= 5.;		// fine roots to non woody litter
const double OMEGA_RC_FWL_TREE				= 5.;		// coarse roots to fine woody litter
const double OMEGA_SF_FWL_TREE				= 200.;		// fine stem biomass to fine woody litter
const double OMEGA_SC_CWL_TREE				= 500.;		// coarse/heavy stem biomas to coarse woody litter

const double OMEGA_LD_NWL_GRASS				= 200.;		// dead leaves to non woody litter
const double OMEGA_RF_NWL_GRASS				= 5.;		// fine roots to non woody litter

const double HELPER_LD_NWL_TREE				= 1./OMEGA_LD_NWL_TREE;
const double HELPER_RF_NWL_TREE				= 1./OMEGA_RF_NWL_TREE;
const double HELPER_RC_FWL_TREE				= 1./OMEGA_RC_FWL_TREE;
const double HELPER_SF_FWL_TREE				= 1./OMEGA_SF_FWL_TREE;
const double HELPER_SC_CWL_TREE				= 1./OMEGA_SC_CWL_TREE;

const double HELPER_LD_NWL_GRASS			= 1./OMEGA_LD_NWL_GRASS;
const double HELPER_RF_NWL_GRASS			= 1./OMEGA_RF_NWL_GRASS;
// End constant for transitions from dead biomass pools to Yasso soil model pools



// Start parameters for water index in Aindex
const double WATER_IND_TREE[NUM_TR_TYPES]		= { 0.66, 0.66, 0.66 };
const double WATER_IND_TREE_1[NUM_TR_TYPES]		= { 1.-WATER_IND_TREE[0], 1.-WATER_IND_TREE[1], 1.-WATER_IND_TREE[2] };
const double WATER_IND_GRASS[NUM_GR_TYPES]		= { 0.66, 0.66, 0.66, 0.66, 0.66, 0.66, 0.66, 0.66 };
const double WATER_IND_GRASS_1[NUM_GR_TYPES]	= { 1.-WATER_IND_GRASS[0], 1.-WATER_IND_GRASS[1], 1.-WATER_IND_GRASS[2], 
                                                    1.-WATER_IND_GRASS[3], 1.-WATER_IND_GRASS[4], 1.-WATER_IND_GRASS[5],
                                                    1.-WATER_IND_GRASS[6], 1.-WATER_IND_GRASS[7] };
// End parameters for water index in Aindex



// Start constants for fire model, Table 14
// constants for fire intensity (Higgins 2008)
const double FIRE_H							= 16890.;
const double FIRE_C							= 301.;
const double FIRE_AW						= 119.7;
const double FIRE_QM						= 2.6e6;
const double FIRE_QV						= 160749.;

// constants for topkill probability (Higgins et al. 2000)
// +- arbitrary values for forest trees that ensure that P(topkill)=1
const double TOP_KILL_CONST[NUM_TR_TYPES]	= { 4.3,      6.3,      6.3      };
const double TOP_KILL_H[NUM_TR_TYPES]		= { 5.003,    3.003,    3.003    };
const double TOP_KILL_I[NUM_TR_TYPES]		= { 0.004408, 0.006408, 0.006408 };


// constants for fire ignition
const double IGNITION_MIN_INT				= 300.;				// minimum fire intensity for ignition
// const double IGNITION_MIN_INT				= 150.;				// minimum fire intensity for ignition
// #ifndef OPTIM_GLOBALS
// const double IGNITION_PROB					= 0.01;				// probability of fire ignition
// const double IGNITION_PAR_2					= 0.1;				// parameter describing fire ignition sequence, treecover
// #endif

const double IGNITION_PAR_1					= -0.083333;		// parameter describing fire ignition sequence, treecover

const double SEED_GERM_PROB[NUM_TR_TYPES]	= { 0.25, 0.25, 0.5   };
// const double SEED_GERM_PROB[NUM_TR_TYPES]	= { 0.25, 0.25, 0.5   };
const double PROB_ROOT_SUCKER[NUM_TR_TYPES]	= { 0.,   0.,   0.    };

const double RESPROUTING_PROB[NUM_TR_TYPES]	= { 0.99, 0.5, 0. };

const double DESIC_COEFF[NUM_GR_TYPES]		= { 0.9, 0.999, 0.9, 0.999, 0.9, 0.999, 0.9, 0.999 };  // C4 dries out faster (Bond 2008)

const double BBS_FIRE_COVER					= 0.6;
// End constants for fire model

// Start constants for reproduction and mortality, Table 13
const double SEED_WEIGHT[NUM_TR_TYPES]		= { 0.001, 0.001, 0.001 };    // weight of a single seed in kg (Hovstadt 1999)
// const double SEED_FRAC_DAY[NUM_TR_TYPES]	= { 0.1, 0.1, 0.01 };
const double SEED_FRAC_DAY[NUM_TR_TYPES]	= { 0.1, 0.1, 0. };
// const double SEED_DECAY_RATE[NUM_TR_TYPES]	= { 0.9967, 0.9967, 0.9967 };       // Higgins 2000
const double SEED_DECAY_RATE[NUM_TR_TYPES]	= { 0.995, 0.995, 0.999 };       // Higgins 2000
const int    REPROD_AGE[NUM_TR_TYPES]       = { 3650, 3650, 365*3 };  // 10 years and 5 years

// End constants for reproduction and mortality, Table 13


// Start constants for soil, Table 15
// const int    SOIL_LAYERS        = 12;    		// Number of soil layers
// const double DEPTH[SOIL_LAYERS] = {      0.05,     0.1,    0.2,    0.3,    0.4,    0.6,    0.8,    1.0,     1.25,     1.5,     1.75,     2.0 };
// const double THICK[SOIL_LAYERS] = { 0.05,     0.05,    0.1,    0.1,    0.1,    0.2,    0.2,    0.2,    0.25,     0.25,    0.25,     0.25,    };
// end constants for soil



// parameters for competitor probability
const double COMP_PAR_1		= 1.1;
const double COMP_PAR_2		= 1.5;
const double COMP_PAR_3		= 1e-4;
// end parameters for competitor probability



// constants to avoid zero in denominator
const double MIN_LAI					= 0.000001;	// if the LAI of a plant is 0, we have a problem with our formulas
const double MIN_HEIGHT					= 0.000001;


// other constants
const int		DAYS_IN_YEAR			= 365;			// Days in a year
const int		MONTH_IN_YEAR			= 12;
const double	GRID_SIZE				= 10000.;		// 1ha


//Penman Monteith constants, source: mostly FAO
const double	SGC						= 0.287;      	// specific gas constant, KJ/kg/K, FOArho}
const double	SP_HEAT					= 1.013e-3;   	// MJ/kg/degC = J/g/degC from FOA specific heat of moist air}
const double	LAMBDA					= 2.45;       	// MJ/kg from FOA latent heat of air at 20degC}
const double	GSC						= 0.082;       	// solar constant - extrateresstrial solar radiation MJ/m2/minute
// const double	GSC						= 0.082*24.;       	// solar constant - extrateresstrial solar radiation MJ/m2/day
const double	ANGSTRONG_A				= 0.25;        	// constant defining prop of extrat. radiation reaching 
// earth on overcast days
const double	ANGSTRONG_B				= 0.50;        	// constant defining the additional prop of extrat. radiation 
// reaching earth on clear days
const double	ALBEDO					= 0.23;        	// canopy reflection coefficient estimate for a grass 
// reference crop (dimensionless)
const double	SBC						= 4.903e-09;   	// stephan bolzman constant MJ.k^-4.day^-1  



// globals from fxn.h
// const double	E25H_CONST_1			= 0.807;
// const double	E25H_CONST_2			= 2.137;
const double	CS_FACTOR				= 1.4;
const double	MAIN_RESP_TRANS			= 3.22;
const double	MAIN_RESP_FAC			= 0.046;
const double	MAIN_RESP_TEMP_TRANS	= 20.;
const double	MAIN_RESP_TEMP_FAC		= 10.;
const double	MAIN_RESP_CARB_PROPORT	= 1.;


// initial values for trees and grass
const double	INIT_MASS_GRASS				= 0.001;				// 10 g per m^2
const double	INIT_MASS_TREE				= 0.001;				// 10 g per plant



// constants to compute moisture of wet biomass for fire, Dennison et al.
const double	MOIST_CONST_A				= 0.24;
const double	MOIST_CONST_B				= 1.09;
const double	MOIST_CONST_C				= 14.2;
const double	MOIST_CONST_D				= 0.78;	
// 
// converts micromol/m^2/s into kg/day/m^2
const double    MMSTOKGD_HELPER             = 1e-6*1e-3*44.*(12./44.)*(1.544)*3600.;


// constants needed for data input
const double	INPUT_GRID_SIZE				= 0.1666667/2.;		// half of grid size of input data ( 10min )
const int		DATA_NUM					= 104;


const int NUM_SIZE_CLASSES = 70;

const int YEARLY_DATA_LENGTH = 56;

const int INPUT_PARAMS = 34;

double const cs_y_vec[42] = { -38.238, -36.372, -34.507, -32.642, -30.777, -28.911, -27.046,
                              -25.181, -23.316, -21.450, -19.585, -17.720, -15.855, -13.989, -12.124, -10.259,
                               -8.394,  -6.528,  -4.663,  -2.798,  -0.933,   0.933,   2.798,   4.663,   6.528,
                                8.394,  10.259,  12.124,  13.989,  15.855,  17.720,  19.585,  21.450,  23.316,
                               25.181,  27.046,  28.911,  30.777,  32.642,  34.507,  36.372,  38.238 };


double const	SFRAC_FINE					= 0.34;		// fraction of stem biomass that is fine wood
double const	SFRAC_COARSE				= 0.25;		// fraction of stem biomass that is coarse wood
double const	SFRAC_HEAVY					= 0.41;		// fraction of stem biomass that is heavy wood

// const int		TIMESTEP_DAILY_DATA			= 7;		// write daily data every NN days
const int		TIMESTEP_DAILY_DATA			= 92;		// write daily data every NN days
// const int		TIMESTEP_DAILY_DATA			= 1;		// write daily data every NN days
// const int		TIMESTEP_DAILY_DATA			= 385;		// write daily data every NN days


#ifdef S_ONLY_FOREST_TREES
const double PROB_FOREST_TREE           = 1.5;          // probability for forest trees (initialization)
#elif defined S_ONLY_SAVANNA_TREES
const double PROB_FOREST_TREE           = -1.5;          // probability for forest trees (initialization)
#else
const double PROB_FOREST_TREE           = 0.5;          // probability for forest trees (initialization)
#endif


#endif























