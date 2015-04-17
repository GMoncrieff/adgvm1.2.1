#ifndef ElephantClass_h___
#define ElephantClass_h___


#include <iostream>

#include "MyMath.h"
#include "ElephantClassGlobals.h"

using namespace std;


class clElephants
{
    private:
        double visit_sequence_[365];            // sequence indicating when elephants visit site
        double daily_consumption_;

        double grass_frac_;

        double visit_ret_;    // visit frequency
        double util_inten_;   // factor for visit intensity
        double tr_mort_pr_;   // tree mortality after utilization (bark stripping)
        double max_tr_con_;   // maximum fraction (0 to 1) of leaf biomass consumed by el.
        double male_prob_;    // probability of bulls (for topkill)
        double wastage_;      // factor to describe wastage during biomass consumption
        double number_;       // elephant number
        double diet_part_;    // diet partitioning parameter
        double daily_requ_;   // daily biomass requirement per elephant (kg)
        double grazing_bm_;   // grazing rate
        double area_;         // area of park

    public:
        clElephants();
        clElephants( double visit_ret, double util_inten, double tr_mort_pr,
                     double max_tr_con, double male_prob, double wastage,
                     double number, double diet_part, double daily_requ,
                     double grazing_bm, double area );
        ~clElephants() {};
        void   GenerateVisitSequence( int day, int year );
        double getVisitSequence( int day )  { return visit_sequence_[day]; }
        double getTreeConsumption( int day );
        double getGrassConsumption( int day );
        void   getGrassFraction( double moisture );
        double getMaleProb()      { return male_prob_; }
        double getTreeMortality() { return tr_mort_pr_; }
        double getMaxTreeCons()   { return max_tr_con_; }
        double getGrazingBm()     { return grazing_bm_; }

};


// ------------------------------------------------------------
// standard constructor, initializes with 1 elephant per km^2
// ------------------------------------------------------------
clElephants::clElephants()
{
    for ( int i=0; i<365; i++ )
        visit_sequence_[i] = 0;
    grass_frac_         = 0;
    daily_consumption_  = 0;
    visit_ret_          = 0;
    util_inten_         = 0;
    tr_mort_pr_         = 0;
    max_tr_con_         = 0;
    male_prob_          = 0;
    wastage_            = 0;
    number_             = 0;
    diet_part_          = 0;
    daily_requ_         = 0;
    grazing_bm_         = 0;
    area_               = 0;
}


clElephants::clElephants( double visit_ret, double util_inten, double tr_mort_pr,
                          double max_tr_con, double male_prob, double wastage,
                          double number, double diet_part, double daily_requ,
                          double grazing_bm, double area )
{
    for ( int i=0; i<365; i++ )
        visit_sequence_[i] = 0;

    grass_frac_         = 0;
    daily_consumption_  = 0;
    visit_ret_          = visit_ret;
    util_inten_         = util_inten;
    tr_mort_pr_         = tr_mort_pr;
    max_tr_con_         = max_tr_con;
    male_prob_          = male_prob;
    wastage_            = wastage;
    number_             = number;
    diet_part_          = diet_part;
    daily_requ_         = daily_requ;
    grazing_bm_         = grazing_bm;
    area_               = area*100.;   // need area in ha
}


// ------------------------------------------------------------
// generates sequence of elephant visits in the grid cell
// ------------------------------------------------------------
// void clElephants::GenerateVisitSequence( int day, int year, double visit_ret,
//                                          double biomass_consumption,
//                                          double util_inten )
void clElephants::GenerateVisitSequence( int day, int year )
{
    if ( day==0 )
    {
        double biomass_consumption = number_*daily_requ_*util_inten_/area_;

        for ( int i=0; i<365; i++ )
            visit_sequence_[i] = 0;

        if ( visit_ret_ == 1. )
            visit_sequence_[(int)floor(365.*drand48())] = biomass_consumption;
        else if (visit_ret_ < 1. & visit_ret_ > 0. )    // less than once a year
        {
            int tmp_freq = MyRound( 1./visit_ret_ );  // visit all xxx years
            if ( year % tmp_freq==0 )
            {
//                 visit_sequence_[(int)floor(365.*drand48())] = (double)tmp_freq*biomass_consumption;
                visit_sequence_[(int)floor(365.*drand48())] = biomass_consumption;
            }
        }
        else if ( visit_ret_==0 );
        else      // more than once a year
        {
            int tmp_freq = MyRound( visit_ret_ );  // visit xxx times per year
            for ( int i=0; i<tmp_freq; i++ )
//                 visit_sequence_[(int)floor(365.*drand48())] += biomass_consumption/(double)tmp_freq;
                visit_sequence_[(int)floor(365.*drand48())] += biomass_consumption;
        }
//         for ( int i=0; i<365; i++ ) cout << visit_sequence_[i] << "  ";
//         cout << endl;
    }


    return;
}


double clElephants::getTreeConsumption( int day )
{
    return visit_sequence_[day]*(1.-grass_frac_)*wastage_;
}

double clElephants::getGrassConsumption( int day )
{
    return visit_sequence_[day]*grass_frac_;
}



void clElephants::getGrassFraction( double moisture )
{
    grass_frac_ = 1./( 1.+exp( (diet_part_-moisture)/ELFS_MOIST_FCT_B ) );

    return;
}


#endif













