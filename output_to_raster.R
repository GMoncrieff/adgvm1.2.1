library(R.utils)
library(raster)

setwd("")

#how many replicates?

nreps<-5
# read in results
res <- read.table("/Users/glennmoncrieff/adgvm1/results/050515_master/result/OutputData/YearlyData.100.0.1.dat", quote="\"")
#res <- read.table("/Users/glennmoncrieff/adgvm1/results/DATE_BRANCH/result/OutputData/YearlyData.100.0.1.dat", quote="\"")


# /*  1 */        yearly_output_data[ 0] = longitude;
# /*  2 */        yearly_output_data[ 1] = latitude;
# /*  3 */        yearly_output_data[ 2] = count_years;
# /*  4 */        yearly_output_data[ 3] = i_scen;
# /*  5 */        yearly_output_data[ 4] = with_fire;
# /*  6 */        yearly_output_data[ 5] = d_seed;
# /*  7 */        yearly_output_data[ 6] = rain_sum,
# /*  8 */        yearly_output_data[ 7] = evapo_sum;
# /*  9 */        yearly_output_data[ 8] = evapo_grass_sum;
# /* 10 */        yearly_output_data[ 9] = evapo_soil_sum;
# /* 11 */        yearly_output_data[10] = MyGrassPop.getBlMaxYearC4()*10.;   // t/ha
# /* 12 */        yearly_output_data[11] = MyGrassPop.getBrMaxYearC4()*10.;   // t/ha
# /* 13 */        yearly_output_data[12] = MyGrassPop.getBlMaxYearC3()*10.;   // t/ha
# /* 14 */        yearly_output_data[13] = MyGrassPop.getBrMaxYearC3()*10.;   // t/ha
# /* 15 */        yearly_output_data[14] = MyGrassPop.getLeafBmDeadSt()*10. + MyGrassPop.getLeafBmDeadLy()*10.;   // t/ha
# /* 16 */        yearly_output_data[15] = MyTreePop.getpCanopySavYearMean();
# /* 17 */        yearly_output_data[16] = MyTreePop.getpCanopyForYearMean();
# /* 18 */        yearly_output_data[17] = MyTreePop.getBlYearMax()*10.;       // t/ha
# /* 19 */        yearly_output_data[18] = MyTreePop.getBsYearMax()*10.;       // t/ha
# /* 20 */        yearly_output_data[19] = MyTreePop.getBrYearMax()*10.;       // t/ha
# /* 21 */        yearly_output_data[20] = MyGrassPop.getC34RatioYearMean();
# /* 22 */        yearly_output_data[21] = MyTreePop.getMeanHeightYearMean();
# /* 23 */        yearly_output_data[22] = MyTreePop.getMaxHeightYearMean();
# /* 24 */        yearly_output_data[23] = MyTreePop.getPopSize();
# /* 25 */        yearly_output_data[24] = MyTreePop.getSavTreeNum();
# /* 26 */        yearly_output_data[25] = MyTreePop.getForTreeNum();
# /* 27 */        yearly_output_data[26] = MyTreePop.getMaxBasalAreaYearMean();
# /* 28 */        yearly_output_data[27] = total_npp/(double)count_years;  // t/ha
# /* 29 */        yearly_output_data[28] = total_nee/(double)count_years;
# /* 30 */        yearly_output_data[29] = mySoil.GetCarbonStored();
# /* 31 */        yearly_output_data[30] = phen_counter;
# /* 32 */        yearly_output_data[31] = fire_num;
# /* 33 */        yearly_output_data[32] = total_fire_intensity/(double)fire_num;
# /* 34 */        yearly_output_data[33] = IData.soil_N_[0];
# /* 35 */        yearly_output_data[34] = IData.soil_C_[0];
# /* 36 */        yearly_output_data[35] = tmp_mean;
# /* 37 */        yearly_output_data[36] = ca_par_preassure*10.;
# /* 38 */        yearly_output_data[37] = d_params[0];
# /* 39 */        yearly_output_data[38] = d_params[1];
# /* 40 */        yearly_output_data[39] = d_params[2];
# /* 41 */        yearly_output_data[40] = d_params[3];
# /* 42 */        yearly_output_data[41] = d_params[4];
# /* 43 */        yearly_output_data[42] = d_params[5];
# /* 44 */        yearly_output_data[43] = d_params[6];
# /* 45 */        yearly_output_data[44] = d_params[7];
# /* 46 */        yearly_output_data[45] = d_params[8];
# /* 47 */        yearly_output_data[46] = d_params[9];
# /* 48 */        yearly_output_data[47] = d_params[10];
# /* 49 */        yearly_output_data[48] = d_params[11];
# /* 50 */        yearly_output_data[49] = d_params[12];
# /* 51 */        yearly_output_data[50] = d_params[13];
# /* 52 */        yearly_output_data[51] = d_params[14];
# /* 53 */        yearly_output_data[52] = d_params[15];
# /* 54 */        yearly_output_data[53] = d_params[16];
# /* 55 */        yearly_output_data[54] = d_params[17];
# /* 56 */        yearly_output_data[55] = d_params[18];
# /* 57 */        yearly_output_data[56] = d_params[19];
# /* 58 */        yearly_output_data[57] = d_params[20];
# /* 59 */        yearly_output_data[58] = d_params[21];
# /* 60 */        yearly_output_data[59] = d_params[22];
# /* 61 */        yearly_output_data[60] = d_params[23];
# /* 62 */        yearly_output_data[61] = d_params[24];
# /* 63 */        yearly_output_data[62] = d_params[25];
# /* 64 */        yearly_output_data[63] = d_params[26];
# /* 65 */        yearly_output_data[64] = d_params[27];
# /* 66 */        yearly_output_data[65] = d_params[28];
# /* 67 */        yearly_output_data[66] = MyTreePop.getActiveDays();
# /* 68 */        yearly_output_data[67] = MyGrassPop.getActiveDays();
# /* 69 */        yearly_output_data[68] = evapo_ref_sum;
# /* 70 */        yearly_output_data[69] = A0C3_mean;
# /* 71 */        yearly_output_data[70] = A0C4_mean;

#variables of interest or VOI

coords  <-cbind(res[,1],res[,2])  #lat lon NB multiple coords per site, because multiple runs
rain    <-res[,7]/res[,3]      #MAP
NPP     <-res[,28]            #NPP
C4live<-rowSums(res[,11:12])      #C4 grass biomass live
C3live<-rowSums(res[,13:14])
C3C4<-res[,21] #C3 C4 ratio
treecover_sav<-res[,16]       #savanna tree canopy cover
treecover_for<-res[,17]       #forest tree canopy cover
tree_biomass<-rowSums(res[,18:20])  #above ground tree biomass
tree_height<-res[,23]         #tree height
phenology_grass<-res[,68]/365 # % days grass active
phenology_tree<-res[,67]/365  # % days tree active

sites<-unique(coords)
sites_string<-paste(as.character(sites[,1])," ",as.character(sites[,2]),sep="")
coord_string<-paste(as.character(coords[,1])," ",as.character(coords[,2]),sep="")

#create matrix from variable with rows as sites and columns as model runs
variable<-matrix(,nrow(sites),nreps)

for ( i in 1:nrow(sites)){  
  temp <- which(coord_string==sites_string[i])
  #replace VOI with name from list above
  ##########
  temp <- C4live[temp]
  #temp <- C3live[temp]
  #temp <- tree_biomass[temp]
  #temp <- treecover_for[temp]
  #temp <- treecover_sav[temp]
  ##########
  variable[i,]<-temp 
}

############################################## clean and plot

ncells<-nrow(variable)
valid_rep<-numeric(ncells)
valid_mean<-numeric(ncells)

for( i in 1:ncells){
  valid_rep[i]<-nreps-sum(is.na(variable[i,]))
}

for( i in 1:ncells){
  valid_mean[i]<-sum(variable[i,],na.rm=T)/valid_rep[i]
}

varXYZ<-cbind(sites,valid_mean)

# treebiomass_master<-rasterFromXYZ(varXYZ,digits=3)
# forcover_master<-rasterFromXYZ(varXYZ,digits=3)
 savcover_master<-rasterFromXYZ(varXYZ,digits=3)
# C4_master<-rasterFromXYZ(varXYZ,digits=3)
# C3_master<-rasterFromXYZ(varXYZ,digits=3)

################################################
#biome map
################################################
#rasters required for biome mapping:

#sav_can  -savnna tree cover
#for_can  -forest tree cover
#C34      -C3 C4 ratio
#grass3   -live c3 grass biomass
#grass4   -live c4 grass biomass
# 
# cancov <- sav_can+for_can
# grass <- grass3+grass4
# fcan <- for_can/cancov
# fcan[Which(cancov==0)]=0

#check all varialbe good and create base template
allinput<-stack(cancov,C34,grass,fcan)
base<-sum(allinput)*0

#create a raster for each biome
#woodland
woodland<-base
woodland[Which(cancov>0.1 & cancov<0.8 & fcan>0.5)]=1

#c4 savanna
savannaC4<-base
savannaC4[Which(cancov>0.1 & cancov<0.8 & fcan<=0.5 & C34 < 0.5)]=2

#c3 savanna
savannaC3<-base
savannaC3[Which(cancov>0.1 & cancov<0.8 & fcan<=0.5 & C34 >= 0.5)]=3


#c4 grassland
grasslandC4<-base
grasslandC4[Which(cancov<=0.1 & grass>0.5 & C34 < 0.5)]=4


#c3 grassland
grasslandC3<-base
grasslandC3[Which(cancov<=0.1 & grass>0.5 & C34 >= 0.5)]=5


#forest
forest<-base
forest[Which(cancov>=0.8)]=6


#desert
desert<-base
desert[Which(cancov<=0.1 & grass<=0.5)]=7

biome_map<-desert+forest+grasslandC3+grasslandC4+savannaC3+savannaC4+woodland
#biome_obs<-raster("/home/glenn/stability/vegmaps/pft_dom.asc")
#aggregate first
#biom_obs<-resample(real,biome_map,method="ngb")
       