// This file defines the data for output in the file
// SysData.dat in the Output directory. Data are written
// every TIMESTEP_DAILY_DATA days, this constant is 
// defined in the Globals.h/GlobalsAustralia.h. These
// data can be analyzed with SysData.R
// Tree size data are written into the file SizeData.dat
// once per year. The file provides the number of trees
// in size classes with 0.5m steps. These data can be
// analyzed with SizeData.R.


if ( day_in_year%TIMESTEP_DAILY_DATA==0 )
{
	double total_combustion = // biomass, kg/m^2
		  MyGrassPop.getLeafLiveCombustion()
		+ MyGrassPop.getLeafDeadStCombustion()
		+ MyGrassPop.getLeafDeadLyCombustion()
		+ MyTreePop.getStemLiveCombustion()
		+ MyTreePop.getLeafLiveCombustion()
		+ MyTreePop.getStemDeadStCombustion()
		+ MyTreePop.getStemDeadLyCombustion()
		+ MyTreePop.getLeafDeadStCombustion()
		+ MyTreePop.getLeafDeadLyCombustion();

//	this is necessary as do not write every day so we cannot set these pools to zero after each fire
	MyTreePop.emptyCombustionPools();
	MyGrassPop.emptyCombustionPools();


// Biomasses given in kg biomass/m^2
//   -- multiplication scales to t biomass/ha
// Fluxes (gpp, rma, rgr, fire) given in kg biomass/m^2
//   -- multiplication with 0.44 scales to kg C/m^2
// Soil flux given in kg C/m^2
//   -- no multiplication required

	_SYSDATA_SysData <<
/* 1*/		setw(14) << count_years <<
/* 2*/		setw(14) << day_in_year <<
/* 3*/		setw(17) << MyGrassPop.getLeafBmLive() <<          // kg biomass/m^2
/* 4*/		setw(17) << MyGrassPop.getRootBmLive() <<          // kg biomass/m^2
/* 5*/		setw(17) << MyGrassPop.getLeafBmDeadSt() <<        // kg biomass/m^2
/* 6*/		setw(17) << MyGrassPop.getLeafBmDeadLy() <<        // kg biomass/m^2
/* 7*/		setw(17) << MyGrassPop.getRootBmDead() <<          // kg biomass/m^2
/* 8*/		setw(17) << MyGrassPop.getGPP()*0.44 <<            // kg C/m^2
/* 9*/		setw(17) << MyGrassPop.getRma()*0.44 <<            // kg C/m^2
/*10*/		setw(17) << MyGrassPop.getRgr()*0.44 <<            // kg C/m^2
/*11*/		setw(17) << MyTreePop.getpCanopySav() <<
/*12*/		setw(17) << MyTreePop.getpCanopyFor() <<
/*13*/		setw(17) << MyTreePop.getLeafBmLive() <<           // kg biomass/m^2
/*14*/		setw(17) << MyTreePop.getStemBmLive() <<           // kg biomass/m^2
/*15*/		setw(17) << MyTreePop.getRootBmLive() <<           // kg biomass/m^2
/*16*/		setw(17) << MyTreePop.getLeafBmDeadSt() <<         // kg biomass/m^2
/*17*/		setw(17) << MyTreePop.getLeafBmDeadLy() <<         // kg biomass/m^2
/*18*/		setw(17) << MyTreePop.getStemBmDeadSt() <<         // kg biomass/m^2
/*19*/		setw(17) << MyTreePop.getStemBmDeadLy() <<         // kg biomass/m^2
/*20*/		setw(17) << MyTreePop.getRootBmDead() <<           // kg biomass/m^2
/*21*/		setw(17) << MyTreePop.getGPP()*0.44 <<             // kg C/m^2
/*22*/		setw(17) << MyTreePop.getRma()*0.44 <<             // kg C/m^2
/*23*/		setw(17) << MyTreePop.getRgr()*0.44 <<             // kg C/m^2
/*24*/		setw(17) << MyTreePop.getMeanLai() <<
/*25*/		setw(10) << MyTreePop.getPopSize() <<
/*26*/		setw(17) << MyTreePop.getMeanHeight() <<
/*27*/		setw(17) << MyGrassPop.getC34Ratio() <<
/*28*/		setw(17) << tmp_mean <<
/*29*/		setw(17) << MyTreePop.getTallTrees() <<
/*30*/		setw(17) << MyTreePop.getBasalArea() <<
/*31*/		setw(17) << ca_par_preassure*10. <<//cs_preassure <<
/*32*/		setw(17) << Rain[day_in_year] <<
/*33*/		setw(17) << EtSite <<
/*34*/		setw(17) << mySoil.GetCarbonRelease() <<           // in C kg/m^2, already in C, no multipl. with 0.44
/*35*/		setw(17) << total_combustion*0.44 <<               // in C kg/m^2
		#ifdef S_MANIP_SEASONALITY
/*36*/		setw(17) << MAN_SEAS_value <<
		#endif
		#ifdef CLIM_SCEN
/*36*/		setw(17) << cs_mean_prec <<
		#endif
		#ifdef S_MANIP_FIRE
/*36*/		setw(17) << glob_fire_param <<
		#endif
		endl;
		
		
		if ( day_in_year==0 )  // write these data only once per year
        {
            for ( int str=0; str<NUM_SIZE_CLASSES; str++ )
                _SIZE_STRUC_size_struc[str] = 0;

            for ( int str=0; str<MyTreePop.getPopSize(); str++ )
                if( MyTreePop.getTreeHeight(str)<35. )
                    _SIZE_STRUC_size_struc[  (int)floor((MyTreePop.getTreeHeight(str)*2.))  ]++;

            for ( int str=0; str<NUM_SIZE_CLASSES; str++ )
                _SIZE_STRUC_SizeStrucFile << setw(7) << _SIZE_STRUC_size_struc[str];
            _SIZE_STRUC_SizeStrucFile << endl;
        }

}






