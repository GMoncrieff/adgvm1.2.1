#ifndef ElephantClassGlobals_h___
#define ElephantClassGlobals_h___



// const double ELFS_MOIST_FCT_A       = 0.5;
// const double ELFS_MOIST_FCT_B       = 7.;
// const double ELFS_MOIST_FCT_A       = 0.5;
const double ELFS_MOIST_FCT_B       = 0.1;

// const double ELFS_FEM_TK_A          = 3.;
// const double ELFS_FEM_TK_B          = 10.;
// const double ELFS_FEM_TK_G          = 0.344840;

// const double ELFS_BUL_TK_A          = 3.;
// const double ELFS_BUL_TK_B          = 10.;
// const double ELFS_BUL_TK_G          = 0.145703;
// const double ELFS_BUL_TK_H          = 0.145703;
// const double ELFS_BUL_TK_S          = 10.;
// const double ELFS_BUL_EXP           = 0.2;

const double ELFS_TK_P = 3;
const double ELFS_TK_B = ELFS_TK_P/1.5;
const double ELFS_TK_M = 10;
const double ELFS_TK_V = 0.8;

const double ELFS_TK_HELPER_1 = pow(ELFS_TK_B,ELFS_TK_P)/MyFak(ELFS_TK_P-1.);
const double ELFS_TK_HELPER_2 = 1./ELFS_TK_V/pow( 2.*3.1415, 0.5 );

const double ELFS_TREE_FAC          = 1.33;  // tree biomass is "cheaper", 200/150

const double ELFS_NUM_USED_TREES    = 35.;   // trees utilized per elephant and day
const double ELFS_DAILY_CONS        = 150.;  // biomass consumption per day
const double ELFS_MAX_BM_PER_TREE   = ELFS_DAILY_CONS/ELFS_NUM_USED_TREES/4.; // 4: leaf, 2*stem, dead leaf

const double ELFS_KNP_AREA          = 1963300; // knp area in ha
// const double ELFS_KNP_AREA          = 1056600; // chobe
// const double ELFS_KNP_AREA          = 3000000; // serengeti


// const double ELFS_MORT_PROB         = 0.25;
// const double ELFS_TREE_CONS         = 0.2;      // maximum fraction of leaf biomass consumed by el.
// const double ELFS_MALE_PROB         = 0.2;      // probability of bulls (for topkill)
// const double ELFS_WASTAGE           = 1.2;      // factor to describe wastage during biomass consumption


// const double ELFS_MORT_PROB         = 0.25*XXX/10000.;
// const double ELFS_TREE_CONS         = 0.2 *YYY/10000.;      // maximum fraction of leaf biomass consumed by el.
// const double ELFS_MALE_PROB         = 0.2 *ZZZ/10000.;      // probability of bulls (for topkill)


#endif












