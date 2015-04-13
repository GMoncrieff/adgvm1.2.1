// This file defines fire output data, these
// data are written into the file FireData.dat
// in the ouput directory and can by analyzed with
// FireData.R. These data are written after each fire.


_FIREDATA_FireData 
	<< setw(17) << fire_num
	<< setw(17) << count_years
	<< setw(17) << day_in_year
	<< setw(17) << dead_fuel
	<< setw(17) << live_fuel
	<< setw(17) << dead_fuel_moisture
	<< setw(17) << live_fuel_moisture
	<< setw(17) << dead_fuel+live_fuel
	<< setw(17) << (live_fuel*live_fuel_moisture + dead_fuel*dead_fuel_moisture)/(live_fuel+dead_fuel)
	<< setw(17) << fire_intensity
	<< setw(17) << patchiness
	<< setw(17) << scorch
	<< setw(17) << cc_fine
	<< setw(17) << cc_coarse
	<< setw(17) << cc_heavy
	<< setw(17) << cc_tk_helper
	<< setw(17) << MyGrassPop.getLeafLiveCombustion()*10.  // t/ha
	<< setw(17) << MyGrassPop.getLeafDeadStCombustion()*10.  // t/ha
	<< setw(17) << MyGrassPop.getLeafDeadLyCombustion()*10.  // t/ha
	<< setw(17) << MyTreePop.getLeafLiveCombustion()*10.  // t/ha
	<< setw(17) << MyTreePop.getLeafDeadStCombustion()*10.  // t/ha
	<< setw(17) << MyTreePop.getLeafDeadLyCombustion()*10.  // t/ha
	<< setw(17) << MyTreePop.getStemLiveCombustion()*10.  // t/ha
	<< setw(17) << MyTreePop.getStemDeadStCombustion()*SFRAC_COARSE*10.  // t/ha
	<< setw(17) << MyTreePop.getStemDeadLyCombustion()*SFRAC_COARSE*10.  // t/ha
	<< setw(17) << MyTreePop.getStemDeadStCombustion()*SFRAC_HEAVY*10.  // t/ha
	<< setw(17) << MyTreePop.getStemDeadLyCombustion()*SFRAC_HEAVY*10.  // t/ha
	<< setw(17) << MyTreePop.getStemDeadStCombustion()*SFRAC_FINE*10.  // t/ha
	<< setw(17) << MyTreePop.getStemDeadLyCombustion()*SFRAC_FINE*10.  // t/ha
	<< setw(17) << N2Ograss(  MyGrassPop.getLeafLiveCombustion() +
							  MyGrassPop.getLeafDeadStCombustion() +
							  MyGrassPop.getLeafDeadLyCombustion() )*10.  // t/ha
	<< setw(17) << N2Oleaf(   MyTreePop.getLeafLiveCombustion() +
							  MyTreePop.getLeafDeadStCombustion() +
							  MyTreePop.getLeafDeadLyCombustion() )*10.  // t/ha
	<< setw(17) << N2Ocoarse( MyTreePop.getStemDeadStCombustion()*SFRAC_COARSE +
							  MyTreePop.getStemDeadLyCombustion()*SFRAC_COARSE )*10.  // t/ha
	<< setw(17) << N2Oheavy(  MyTreePop.getStemDeadStCombustion()*SFRAC_HEAVY +
							  MyTreePop.getStemDeadLyCombustion()*SFRAC_HEAVY )*10.  // t/ha
	<< setw(17) << N2Oshrub(  MyTreePop.getStemLiveCombustion() +
							  MyTreePop.getStemDeadStCombustion()*SFRAC_FINE +
							  MyTreePop.getStemDeadLyCombustion()*SFRAC_FINE )*10.  // t/ha
	<< setw(17) << CH4grass(  MyGrassPop.getLeafLiveCombustion() +
							  MyGrassPop.getLeafDeadStCombustion() +
							  MyGrassPop.getLeafDeadLyCombustion() )*10.  // t/ha
	<< setw(17) << CH4leaf(   MyTreePop.getLeafLiveCombustion() +
							  MyTreePop.getLeafDeadStCombustion() +
							  MyTreePop.getLeafDeadLyCombustion() )*10.  // t/ha
	<< setw(17) << CH4coarse( MyTreePop.getStemDeadStCombustion()*SFRAC_COARSE +
							  MyTreePop.getStemDeadLyCombustion()*SFRAC_COARSE )*10.  // t/ha
	<< setw(17) << CH4heavy(  MyTreePop.getStemDeadStCombustion()*SFRAC_HEAVY +
							  MyTreePop.getStemDeadLyCombustion()*SFRAC_HEAVY )*10.  // t/ha
	<< setw(17) << CH4shrub(  MyTreePop.getStemLiveCombustion() +
							  MyTreePop.getStemDeadStCombustion()*SFRAC_FINE +
							  MyTreePop.getStemDeadLyCombustion()*SFRAC_FINE )*10.  // t/ha
	<< setw(17) << ca_par_preassure*10.
	<< endl;














