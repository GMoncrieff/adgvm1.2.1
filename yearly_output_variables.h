// Definition of vector for annual data output, these data are written
// into YearlyData.***.dat in the output directory. If the length of
// this vector is changed, the constant YEARLY_DATA_LENGTH in the
// Globals.h/GlobalsAustralia.h must be changed.

/*  1 */        yearly_output_data[ 0] = longitude;
/*  2 */        yearly_output_data[ 1] = latitude;
/*  3 */        yearly_output_data[ 2] = count_years;
/*  4 */        yearly_output_data[ 3] = i_scen;
/*  5 */        yearly_output_data[ 4] = with_fire;
/*  6 */        yearly_output_data[ 5] = d_seed;
/*  7 */        yearly_output_data[ 6] = rain_sum,
/*  8 */        yearly_output_data[ 7] = evapo_sum;
/*  9 */        yearly_output_data[ 8] = evapo_grass_sum;
/* 10 */        yearly_output_data[ 9] = evapo_soil_sum;
/* 11 */        yearly_output_data[10] = MyGrassPop.getBlMaxYearC4()*10.;   // t/ha
/* 12 */        yearly_output_data[11] = MyGrassPop.getBrMaxYearC4()*10.;   // t/ha
/* 13 */        yearly_output_data[12] = MyGrassPop.getBlMaxYearC3()*10.;   // t/ha
/* 14 */        yearly_output_data[13] = MyGrassPop.getBrMaxYearC3()*10.;   // t/ha
/* 15 */        yearly_output_data[14] = MyGrassPop.getLeafBmDeadSt()*10. + MyGrassPop.getLeafBmDeadLy()*10.;   // t/ha
/* 16 */        yearly_output_data[15] = MyTreePop.getpCanopyYearMean(TR_SAV);
/* 17 */        yearly_output_data[16] = MyTreePop.getpCanopyYearMean(TR_FOR);
/* 18 */        yearly_output_data[17] = MyTreePop.getpCanopyYearMean(TR_BBS);
/* 19 */        yearly_output_data[18] = MyTreePop.getBlYearMax()*10.;       // t/ha
/* 20 */        yearly_output_data[19] = MyTreePop.getBsYearMax()*10.;       // t/ha
/* 21 */        yearly_output_data[20] = MyTreePop.getBrYearMax()*10.;       // t/ha
/* 22 */        yearly_output_data[21] = MyGrassPop.getC34RatioYearMean();
/* 23 */        yearly_output_data[22] = MyTreePop.getMeanHeightYearMean();
/* 24 */        yearly_output_data[23] = MyTreePop.getMaxHeightYearMean();
/* 25 */        yearly_output_data[24] = MyTreePop.getPopSize();
/* 26 */        yearly_output_data[25] = MyTreePop.getTreeNumType(TR_SAV);
/* 27 */        yearly_output_data[26] = MyTreePop.getTreeNumType(TR_FOR);
/* 28 */        yearly_output_data[27] = MyTreePop.getTreeNumType(TR_BBS);
/* 29 */        yearly_output_data[28] = MyTreePop.getMaxBasalAreaYearMean();
/* 30 */        yearly_output_data[29] = total_npp/(double)count_years;  // kg C/m^2
/* 31 */        yearly_output_data[30] = total_nee/(double)count_years;
/* 32 */        yearly_output_data[31] = mySoil.GetCarbonStored();
/* 33 */        yearly_output_data[32] = phen_counter;
/* 34 */        yearly_output_data[33] = fire_num;
/* 35 */        yearly_output_data[34] = total_fire_intensity/(double)fire_num;
/* 36 */        yearly_output_data[35] = IData.soil_N_[0];
/* 37 */        yearly_output_data[36] = IData.soil_C_[0];
/* 38 */        yearly_output_data[37] = tmp_mean;
/* 39 */        yearly_output_data[38] = ca_par_preassure*10.;
/* 40 */        yearly_output_data[39] = d_params[0];
/* 41 */        yearly_output_data[40] = d_params[1];
/* 42 */        yearly_output_data[41] = d_params[2];
/* 43 */        yearly_output_data[42] = d_params[3];
/* 44 */        yearly_output_data[43] = d_params[4];
/* 45 */        yearly_output_data[44] = d_params[5];
/* 46 */        yearly_output_data[45] = d_params[6];
/* 47 */        yearly_output_data[46] = d_params[7];
/* 48 */        yearly_output_data[47] = d_params[8];
/* 49 */        yearly_output_data[48] = d_params[9];
/* 50 */        yearly_output_data[49] = MyTreePop.getActiveDays();
/* 51 */        yearly_output_data[50] = MyGrassPop.getActiveDays();
/* 52 */        yearly_output_data[51] = evapo_ref_sum;
/* 53 */        yearly_output_data[52] = A0C3_mean;
/* 54 */        yearly_output_data[53] = A0C4_mean;
/* 55 */        yearly_output_data[54] = gs_C3_global;
/* 56 */        yearly_output_data[55] = gs_C4_global;
    






















