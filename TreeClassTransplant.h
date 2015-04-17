#ifndef TreeClass_h___
#define TreeClass_h___


#include <iostream>

#include "GlobalsTransplant.h"
#include "Fxn.h"


using namespace std;


class clTree
{
	private:
		int		Dpos_; 								// dormancy counter (days where nCGT>0)
		int		Dneg_; 								// standby counter (days where nCGT=0)
		int		Dormant_; 							// dormant yes =1, dormant no=0
		int		competitor_;						// index of tree for lightcompetition
		int		number_;							// number of tree
		int		age_;								// age of tree in days
		int		reprod_strategy_;					// reproduction strategy, 0=seeds, 1=root sucker
		
		double	Bl_;    							// leaf biomass, kg per plant
		double	Br_;    							// root biomass, kg per plant
		double	Bs_;    							// stem biomass, kg per plant
		
        double	Bls_;								// dead leaf biomass at tree
        double	Bld_;								// dead leaf biomass on ground, kg per plant
        double	Brd_;								// dead root biomass, kg per plant
        double	Bsd_;								// dead stem biomass, kg per plant

        double  gpp_;
        double  Rma_;
        double  Rgr_;

		double	Qsum_;
		double	Qi_;
		double	Ci_;
		double	gc_;    							// canopy stomatal conductance (mu mol/m2/s)   molar
		double	gb_;								// canopy boundary layer conductanced (m/s)    non-molar
		double	death_prob_;						// probability, that this tree will die
		double	Droot_;   							// rooting depth m
		double	nCGT_;
		double	Gw_;
		double	alloc_denom_;
		double	Aindex_;
		double	height_;							// height m
		double	canopy_radius_;
		double	canopy_area_;						// canopy area m^2
		double	LAI_;	 							// leaf area index of this tree
		double	site_param_;						// "quality" of site where plant growth
		double	stem_area_;							// in m^2
		
		
	private:
		void calDroot();
		void calPlantHeight();
		void calCanopyRadius();
		void calCanopyArea();
		void calStemArea();
		void calPlantLAI();
		void calSoilMoisture( arryLayer G_theta );
		void calGc( double A, double resp, double Ca, double P, double hs, double T );
		void calGb( double wind );
		void AllocateCorbon();
		void UpdateAllometry();
		void MassTurnover();
		void Respirate( double rRoot, double rStem, double rLeaf );
		double AllocLeaf();
		double AllocStem();
		double AllocRoot();
		double TopKillProb( double intensity );
		
		
	public:
		clTree( ) {}
		clTree( double init_mass, double pCanopy, int popsize, int number );
		~clTree();
		int  WillIDie(int frost);
        void setStateAfterElephants();
		void setStateAfterFire( double intensity );
		void setDeadBiomass();
		void RunPhysiology( double A0, double RmL, double mmsTOkgd, double T, double wind,
							arryLayer G_theta, double Ca, double P, double rh,
							double comp_height, double grass_height, int day,
							double T_fac, int frost );
		double getEt( double T, double P, double radnet_day, double s12, double SP_HEAT,
					  double gama12, double rho12, double VPD12 );
		int getSeedProduction();
		int getRootSucker();
		
		int    getCompetitor()			{ return competitor_; }
		int    getNumber()				{ return number_; }
		int    getDormant()				{ return Dormant_; }
		int    getReprodStrategy()		{ return reprod_strategy_; }
		double getBl()					{ return Bl_; }
		double getBr()					{ return Br_; }
		double getBs()					{ return Bs_; }
        double getBls()                 { return Bls_; }
        double getBld()					{ double tmp=Bld_; Bld_=0; return tmp; }
        double getBrd()					{ double tmp=Brd_; Brd_=0; return tmp; }
        double getBsd()					{ double tmp=Bsd_; Bsd_=0; return tmp; }
        double getGPP()                 { return gpp_; }
        double getRma()                 { return Rma_; }
        double getRgr()                 { return Rgr_; }
		double getQsum()				{ return Qsum_; }
		double getQi()					{ return Qi_; }
		double getCi()					{ return Ci_; }
		double getGw()					{ return Gw_; }
		double getGc()					{ return gc_; }
		double getGb()					{ return gb_; }
		double getDeathProb()			{ return death_prob_; }
		double getDroot()				{ return Droot_; }
		double getnCGT()				{ return nCGT_; }
		double getAindex()				{ return Aindex_; }
		double getHeight()				{ return height_; }
		double getCanopyRadius()		{ return canopy_radius_; }
		double getCanopyArea()			{ return canopy_area_; }
		double getStemArea()			{ return stem_area_; }
		double getLAI()					{ return LAI_; }
		
        void   setBs( double val )		{ Bs_  = val; }
        void   setBl( double val )		{ Bl_  = val; }
        void   setBls( double val )		{ Bls_ = val; }

};



// -------------------------------------------------------------------------------------------------------------------
// This constructor generates a tree of random size between 0 and max_Bl, whereas Bl is the leaf biomass.
// Bs and Br are in a constant ratio to Bl. 

clTree::clTree( double init_mass, double pCanopy, int popsize, int number )
{
	number_				=  number;
	age_				=  365*(int)init_mass;
	Dpos_				=  0;
	Dneg_				=  0;
	Dormant_			=  1;
	competitor_			= -1;
	
	Bl_					= init_mass*0.1;
	Br_					= init_mass*0.45;
	Bs_					= init_mass*0.45;
	Bld_				= 0.;
	Bls_                = 0.;
	Brd_				= 0.;
	Bsd_				= 0.;
	gpp_                = 0.;
	Rma_                = 0.;
	Rgr_                = 0.;
	Qsum_				= 1.;
	Qi_					= 1.;
	Ci_					= 1.;
	gc_					= 0;
	gb_					= 0;
	death_prob_			= 0;
	nCGT_				= 0.;
	Gw_					= 1.;
	Aindex_				= .2;
	canopy_area_		= MIN_CANOPY_AREA;
	stem_area_			= 0.00001;
	site_param_			= 0;
	
	reprod_strategy_	= 0;
    if ( drand48()<PROB_ROOT_SUCKER ) reprod_strategy_ = 1;
	
	
	// compute some allometric measures, functions in Fxn
	UpdateAllometry();
	
// we determine a neighbour tree for light competition.
	double has_comp = drand48();
	if ( has_comp > COMP_PAR_1-COMP_PAR_2*pCanopy-popsize*COMP_PAR_3 )
		competitor_ = (int) floor(popsize*drand48());
}



// -------------------------------------------------------------------------------------------------------------------

clTree::~clTree()
{
}



void clTree::calDroot()
{
	Droot_ = MyMin( MAX_ROOT_DEP_TREE, CylinderDepth( Br_, MIN_ROOT_RAD_TREE, ROOT_DENSITY_TREE, MAX_ROOT_DEP_TREE) );
	return;
}


// -------------------------------------------------------------------------------------------------------------------
// Standard power function, based on fits to Niklas&Enquist data, estimates the plant height from the stem biomass.
void clTree::calPlantHeight()
{
#   ifdef TRANS_ALLOMETRY
	height_ = exp( (log(Bs_)+3.5413)/3.5337 );               // Williams 2005
#	else
	height_ = HEIGHT_C1_TREE*pow(Bs_+Bl_,HEIGHT_C2_TREE);
#   endif
	return;
}


// -------------------------------------------------------------------------------------------------------------------
// Compute radius of the canopy
void clTree::calCanopyRadius()
{
#   ifdef TRANS_ALLOMETRY
	canopy_radius_ = sqrt( 3.0623641*pow(Bs_, 0.6616535)/3.141583 )/2.;
#	else
	canopy_radius_ = height_*GAMMA_CANOPY_TREE;
#   endif
	return;
}

// -------------------------------------------------------------------------------------------------------------------
// Compute the canopy area of a tree,
// assumes canopy radius is constant proportion of height
void clTree::calCanopyArea()
{
	canopy_area_ = MIN_CANOPY_AREA+M_PI*canopy_radius_*canopy_radius_;
	
	return;
}


// -------------------------------------------------------------------------------------------------------------------


void clTree::calStemArea()
{
#	ifdef TRANS_ALLOMETRY
	stem_area_ = exp( (log( Bs_ )+2.2111)/2.4831 )/100./2.;   // stem radius in m, Williams 2005
	stem_area_ = stem_area_*stem_area_*3.141593;              // stem area in m^2, Williams 2005
#	else
	stem_area_ = STEM_AREA_C1*height_/STEM_AREA_C2;    // from Higgins 2007 Ecology; 2000 to get radius in m
	stem_area_ = stem_area_*stem_area_*3.141593;    //  from Higgins 2007 Ecology;
#	endif
	return;
}

// -------------------------------------------------------------------------------------------------------------------
// compute leaf area index
void clTree::calPlantLAI()
{
	LAI_ = Bl_*SLA_TREE/MyMax(0.1,canopy_area_)+MIN_LAI;
	return;
}

// -------------------------------------------------------------------------------------------------------------------
// compute the allometric values of the tree
void clTree::UpdateAllometry()
{
	calDroot();
	calPlantHeight();
	calCanopyRadius();
	calCanopyArea();
	calStemArea();
	calPlantLAI();
	
	return;
}


// ---------------------------------------------------------------------------------------------------------------
int clTree::WillIDie(int frost)
{
	int will_I       = 0;
	double rnum = drand48();
	
	death_prob_ = 0;
	
	if ( frost==1 )                         death_prob_ += DEATH_PROB_FROST;
	if ( Dormant_==0 && nCGT_<Rma_ )        death_prob_ += DEATH_PROB_CARBON;
	if ( competitor_>0 )                    death_prob_ += DEATH_PROB_COMP;
	
	if ( rnum < death_prob_ ) will_I = 1;
	
#	ifdef S_NOMORT
	will_I = 0;
#	endif
	
	if ( number_==0 && will_I==1 )
	{
		Bl_ = 0.1;
		Bs_ = 0.1;
		Br_ = 0.1;
	}
	
	return will_I;
}


// -------------------------------------------------------------------------------------------------------------------
// allocation function of net carbon gain to the roots
double clTree::AllocRoot()
{
	return  (1.+A0_ROOT_TREE-Gw_)*alloc_denom_;
}

// -------------------------------------------------------------------------------------------------------------------
// allocation function of net carbon gain to the stem
double clTree::AllocStem()
{
	return  (1.+A0_STEM_TREE-Qi_)*alloc_denom_;
}

// -------------------------------------------------------------------------------------------------------------------
// allocation function of net carbon gain to the leaf
double clTree::AllocLeaf()
{
	return  (1.-Ci_)*alloc_denom_;
}

// -------------------------------------------------------------------------------------------------------------------
// call the alloaction functions for leaf, stem and root, and add it to the tree biomass
void clTree::AllocateCorbon()
{
	double al=0;
	double as=0;
	double ar=0;
	
	if ( nCGT_>0 && Dormant_==0 )
	{
		gpp_ = nCGT_;
		
		alloc_denom_ = 1./(3.+A0_ROOT_TREE+A0_STEM_TREE-Qi_-Gw_-Ci_);
		
		Bl_ += nCGT_*AllocLeaf();
		Bs_ += nCGT_*AllocStem();
		Br_ += nCGT_*AllocRoot();
	}
	
	
	return;
}


// -------------------------------------------------------------------------------------------------------------------
// run the tree physiology

void clTree::RunPhysiology( double A0, double RmL, double mmsTOkgd, double T, double wind,
							arryLayer G_theta, double Ca, double P, double rh,
							double comp_height, double grass_height, int day,
							double T_fac, int frost )
{
	double Ac;
	double Acs;			// canopy and canopy-stressed photosythesis
	double CG;
	double RmS;
	double RmR;			// maint respiration total, stem, root
	double RmLc;			// canopy leaf maint resp
	double RmLg;			// leaf maint resp in kg/day/plant
	
	gpp_                = 0.;
	Rma_                = 0.;
	Rgr_                = 0.;

	// water availability
	calSoilMoisture( G_theta );
	
	// light availability
	Qi_ = 1.;
	if ( grass_height > height_ ) { Qi_ *= LIGHT_COMP_GRASS_TREE_1/grass_height*height_+LIGHT_COMP_GRASS_TREE_2; }
	if ( comp_height  > height_ ) { Qi_ *= LIGHT_COMP_TREE_TREE_1/comp_height*height_+LIGHT_COMP_TREE_TREE_2; }
	Qsum_ = Qi_*FQsumTree( LAI_ );
	
	Ci_ = MyMin(1.,Bl_/(Bl_+Br_+Bs_)/A0_LEAF_TREE);
	
	// conopy photosynthesis
	Ac      = A0*Qsum_;		// light stresses canopy photosynthesis (mmol/m^2/s)
	Acs     = Ac*Gw_; 		// water (and light) stressed canopy photosynthesis (=A0*Qsum*Gw) (mmol/m^2/s)
	
	// carbon gain
	CG  = Acs*mmsTOkgd; //water stressed carbon gain kg/day/m^2 // NOTE
	CG  = MyMax( CG, CG*canopy_area_ );
	
	// components of maintenance respiration
	double bl_bs_ratio = Bl_/Bs_;
	RmLc = R_MAINT_RESP_C3*Acs; //main resp leaf is now function of water stress and canopy photosynthesis (mmol/m^2/s)
	RmLg = RmLc*mmsTOkgd; //leaf maint resp in kg/day/plant  NOTE
	RmLg = MyMax( RmLg, RmLg*canopy_area_ );
	RmS  = MaintRespFast( BETA_STEM_TREE*Bs_, BETA_N, UPSILON_STEM_TREE, T_fac );
	RmR  = MaintRespFast( BETA_ROOT_TREE*(Br_*(1.-bl_bs_ratio)), BETA_N, UPSILON_ROOT_TREE, T_fac )   // coarse roots
		                                + Br_*    bl_bs_ratio*RmLg;                                   // fine roots
	
	//growth resp kg/day
	Rgr_   = SIGMA_GROW_RESP_TREE*CG;
	
	//canopy conductance
	calGb( wind );  // m/s
	calGc( Acs, RmLc, Ca, P, rh, T );  // (mmol/m^2/s)
	
	double Ti;
	if      ( tmp_min_month > TI_CONST ) Ti = 0.;
	else                                 Ti = 1.5*(-1.+tmp_min_month/TI_CONST);
	
	Aindex_ = A0*( (WATER_IND_TREE*G_theta[3]+WATER_IND_TREE_1*Gw_) + Ti ) - RmL + site_param_;
	
	if ( Dormant_==0 )
	{
		nCGT_ = CG-Rgr_-RmR-RmS-RmLg;
		if ( nCGT_ < 0 )
		{
			nCGT_ = 0;
			Respirate( RmR, RmS, RmLg );
		}
	}
	else
	{
		nCGT_ = -0.01;
		Respirate( RmR, RmS, RmLg );
	}
	
	// if no carbon gain increment days without carbon gain
	if ( Dormant_==0 && Aindex_<=STRESS_INDEX_TREE ) Dneg_++;
	// if carbon gain then reset days without carbon gain counter
	if ( Dormant_==0 && Aindex_>STRESS_INDEX_TREE )  Dneg_ = 0; 
	// if dormant and carbon gain then increment dormancy break counter
	if ( Dormant_==1 && Aindex_>STRESS_INDEX_TREE  ) Dpos_++;
	// if dormant and carbon gain = 0 then reset dormancy break counter
	if ( Dormant_==1 && Aindex_<=STRESS_INDEX_TREE ) Dpos_ = 0; 
	// days with frost increase number of neg. days
	if ( frost==1 ) { Dneg_++; Dpos_=0; }
	// if in growth stage then dpos (dormancy break counter) remains zero.
	if ( Dormant_==0 ) Dpos_ = 0; 
	if ( Dormant_==1 ) Dneg_ = 0;
	
	
	//if days of negative nCGT > threshold then assign to dormant state
	if ( Dneg_>=D_NEG_TREE )
	{
		if (number_==0) phen_counter++;
		Dormant_   = 1;
		Dneg_      = 0;
        Bls_      += (1.-REM_BM_TREE)*Bl_;           // NOTE new version: permament turnover and litterfall
        Bl_        = REM_BM_TREE*Bl_;                // NOTE new version: permament turnover and litterfall
	}
    
	// if days of positive nCGT > threshold
	if ( Dpos_>=D_POS_TREE )
	{
		Dormant_ = 0; 
		Dpos_ = 0;
	}
	
    if ( Dormant_==0 )
    {
        Bls_ += HELPER_BLS_TREE * Bl_;
        Bl_  *= HELPER_BL_TREE;
    }
	
	// biomass turnover
	MassTurnover();
	
	// allocate carbon gain
	AllocateCorbon();
	
	UpdateAllometry();
	
	age_++;
	
	return;
}



// ---------------------------------------------------------------------------------------------------------------------------
// compute the water aviability of a tree. The value is between 0 and 1 and is influenced by the
// soil moisture and the rooting depth

void clTree::calSoilMoisture( arryLayer G_theta )
{
	int layers = 1;    // counts the number of layers, in which the plant has roots
	Gw_ = G_theta[0]*THICK[0];
	
	
	while ( Droot_ > DEPTH[layers-1] )
	{
		Gw_ += G_theta[layers]*THICK[layers];
		layers++;
	}
	
	Gw_ /= DEPTH[layers-1];
	
	return;
}


// --------------------------------------------------------------------------------------------------------------
void clTree::MassTurnover()
{
// 	new model for turnover: no stem turnover, fine root turnover equals leaf turnover
    
    
    double bl_bs_ratio = Bl_/Bs_;
    double f_root_turnover =     bl_bs_ratio *Br_*HELPER_BLD_TREE;
    double c_root_turnover = (1.-bl_bs_ratio)*Br_*HELPER_BRD_TREE;

    Bld_  = OMEGA_TRANS_TREE*Bls_;
    Bls_ -= OMEGA_TRANS_TREE*Bls_;

    Brd_  = (f_root_turnover+c_root_turnover);
    Br_  -= (f_root_turnover+c_root_turnover);

    return;
}


// --------------------------------------------------------------------------------------------------------------
void clTree::Respirate( double rRoot, double rStem, double rLeaf )
{
    double tmp_bm = Br_+Bs_+Bl_;

    Br_  = MyMax( 0.0001, Br_-rRoot );
    Bs_  = MyMax( 0.0001, Bs_-rStem );
    Bl_  = MyMax( 0.0001, Bl_-rLeaf );

    Rma_ = tmp_bm - (Br_+Bs_+Bl_);

    return;
}

// --------------------------------------------------------------------------------------------------------------


void clTree::setStateAfterFire( double intensity )
{
	double omega_top_kill;
	
#	ifdef TRANS_TOPKILL
	double p_tk_1 = 1.- ( 1. - 0.098*intensity/1000. + 0.037*intensity/1000.*log( 100.*sqrt(stem_area_/3.141593) ) );
	double p_tk_2 = 1.- ( 1. + 0.029*intensity/1000. - 0.010*intensity/1000.*log( 100.*sqrt(stem_area_/3.141593) ) );
	omega_top_kill = MyMax( p_tk_1, p_tk_2 );
#	else
	omega_top_kill = TopKillProb( intensity );
#	endif
	
	Bld_ = 0;
	Bsd_ = 0;
	
	if ( omega_top_kill > drand48() )
	{
		Bs_      = 0.01*Bs_;
		Bl_      = 0.001*Bl_;
        Bls_     = 0.001*Bls_;
		Dpos_    = 0;
		Dneg_    = 0;
		UpdateAllometry();
	}
	
	return;
}

// --------------------------------------------------------------------------------------------------------------


void clTree::setStateAfterElephants()
{
    
    Bld_     = 0.95*Bl_+0.95*Bls_;
    Bl_      = 0.05*Bl_;
    Bls_     = 0.05*Bls_;
    Bsd_     = 0.95*Bs_;
    Bs_      = 0.05*Bs_;
    Dpos_    = 0;
    Dneg_    = 0;
    UpdateAllometry();
    
    return;
}


// ---------------------------------------------------------------------------------------------------------------------------
// Function to compute the topkill probability for a tree
double clTree::TopKillProb( double intensity )
{
	double prob=1;
	
	if ( height_>0 )
		prob = (   exp( TOP_KILL_CONST - TOP_KILL_H*log( height_ ) + TOP_KILL_I*sqrt( intensity ) ) )/
		      ( 1.+exp( TOP_KILL_CONST - TOP_KILL_H*log( height_ ) + TOP_KILL_I*sqrt( intensity ) ) );
			  
	return prob;
}


void clTree::setDeadBiomass()
{
	Bld_ = 0;
// 	Bsd_ = 0;
	Brd_ = 0;
	return;
}


int clTree::getSeedProduction()
{
	int nSeeds = 0;
	
	if ( reprod_strategy_==0 && age_>3650 && nCGT_>0 )
	{
		nSeeds = (int)floor( nCGT_/(SEED_WEIGHT) );
	}
	
	return nSeeds;
}


int clTree::getRootSucker()
{
	int ret=0;
	
// 	if ( reprod_strategy_==1 && Br_>50 && drand48()>0.25 )
	if ( reprod_strategy_==1 && Br_>20 )
	{
		ret=1;
		Br_ -= 1.;
	}
	
	return ret;
}


// calculate gc_ in molar units
void clTree::calGc( double A, double resp, double Ca, double P, double hs, double T )
{
	double cs;
	double gb = gMolar( T, gb_, P );    // convert gb_ from non-molar to molar
	
	
	if ( gb>0 )
	{
		cs = Ca-CS_FACTOR*A*P/gb; 							//this line follows Arora
		gc_ = MyMax(0., M_TREE_C3*(A-resp)*hs/100.*P/cs+B_TREE_C3);	//ie for negative An we produce zero gs
	}
	else
		gc_ = 0.;
	
	return;
}


void clTree::calGb( double wind )
{
	gb_ = GetgbCANOPY( wind, height_ );
	
	return;
}



double clTree::getEt( double T, double P, double radnet_day, double s12, double SP_HEAT,
					  double gama12, double rho12, double VPD12 )
{
	double Et_tree = 0;
	double gc_nmolar  = gNonMolar( T, gc_, P );  //gn
	
	Et_tree = FOAPenMan( radnet_day, gb_, gc_nmolar, s12, SP_HEAT, gama12, rho12, VPD12 );  //mm/day
	
	if ( Dormant_==1 ) Et_tree=0;
	
	return Et_tree*canopy_area_;
}




#endif
















