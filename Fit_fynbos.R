source("MyDEoptim.R")

# compile code
# system("g++ -o Model_all Model.cpp -DW_SYSDATA -DCLIM_SCEN -DOPTIM_GLOBALS")
system("g++ -o Model_all Model.cpp -DCLIM_SCEN -DOPTIM_GLOBALS")


# v.des  1
# v.c4g  2
# v.c3g  3
# v.c4s  4
# v.wdl  5
# v.for  6
# v.c3s  7
# v.fyn  8

xcoo <- c(  19,  31.58,  18.5,  22,  31.266,  26,  28, 32.5 )
ycoo <- c( -34, -24.99, -32.5, -34, -25.166, -30, -28,  -27 )
vtpy <- c(   8,      4,     8,   8,       4,   2,   2,    6 )


nruns <- 3
nproc <- 12

getVegType <- function( d )
{
	# vectors for classification of cells
	v.c3s <- rep( 0, dim(d)[1] )
	v.c4s <- rep( 0, dim(d)[1] )
	v.for <- rep( 0, dim(d)[1] )
	v.c3g <- rep( 0, dim(d)[1] )
	v.c4g <- rep( 0, dim(d)[1] )
	v.wdl <- rep( 0, dim(d)[1] )
	v.des <- rep( 0, dim(d)[1] )
	v.fyn <- rep( 0, dim(d)[1] )
	
	# do vegetation classification for year 300
	for ( i in 1:dim(d)[1] )
	{
		if      ( d[i,16]+d[i,17]+d[i,18] < 0.1 & d[i,11]+d[i,13] < 1.5 ) v.des[i] <- v.des[i]+1    # desert
		else if ( d[i,16]+d[i,17]+d[i,18] < 0.1 & d[i,11]+d[i,13] > 1.5 )       # grasslands
		{
			if ( d[i,22]<0.5 )                            v.c4g[i] <- v.c4g[i]+1
# 			else                                          v.c3g[i] <- v.c3g[i]+1
			else                                          v.c4g[i] <- v.c4g[i]+1  ## do not distinguish between C3 and C4 grasslands
		}
		else if ( d[i,16]+d[i,17] > 0.8 )                 v.for[i] <- v.for[i]+1  # forest
# 		else if ( d[i,18] > 0.1 & d[i,16]+d[i,17] < 0.1 ) v.fyn[i] <- v.fyn[i]+1  # fynbos
		else if ( d[i,18] > 0.6 )                         v.fyn[i] <- v.fyn[i]+1  # fynbos
		else
		{
# 			if   ( d[i,32]> 0 & d[i,22]<0.5) v.c4s[i] <- v.c4s[i]+1
			if      ( d[i,16]>d[i,17] & d[i,22]<0.5 ) v.c4s[i] <- v.c4s[i]+1
			else if ( d[i,16]>d[i,17] & d[i,22]>0.5 ) v.c3s[i] <- v.c3s[i]+1
			else                                      v.for[i] <- v.for[i]+1  ## do not distinguish between forest and woodland
# 			else                                      v.wdl[i] <- v.wdl[i]+1
		}
	}
	vegtype <- cbind( v.des, v.c4g, v.c3g, v.c4s, v.wdl, v.for, v.c3s, v.fyn )
	vegtype_vec <- rep(0,dim(vegtype)[1])
	for ( i in 1:length(vegtype_vec) ) vegtype_vec[i] <- which(vegtype[i,]==1)

	d <- cbind( d[,1],d[,2], vegtype_vec )
	
	return( d )
}






descr <- c(
"DEATH_PROB_FROST0  ",
"DEATH_PROB_CARBON0 ",
"DEATH_PROB_COMP0   ",
"DEATH_PROB_FROST1  ",
"DEATH_PROB_CARBON1 ",
"DEATH_PROB_COMP1   ",
"DEATH_PROB_FROST2  ",
"DEATH_PROB_CARBON2 ",
"DEATH_PROB_COMP2   ",
"IGNITION_PROB      ",
"IGNITION_PAR_2     ",
"L_SAV_SAV          ",
"L_SAV_FOR          ",
"L_SAV_BBS          ",
"L_FOR_SAV          ",
"L_FOR_FOR          ",
"L_FOR_BBS          ",
"L_BBS_SAV          ",
"L_BBS_FOR          ",
"L_BBS_BBS          ",
"L_C4G_SAV          ",
"L_C4G_FOR          ",
"L_C4G_BBS          ",
"L_C3G_SAV          ",
"L_C3G_FOR          ",
"L_C3G_BBS          ",
"L_C4G_C4G          ",
"L_C3G_C3G          ",
"L_SAV_C4G          ",
"L_SAV_C3G          ",
"L_FOR_C4G          ",
"L_FOR_C3G          ",
"L_BBS_C4G          ",
"L_BBS_C3G          "
)

parset <- c(
0.001,
0.001,
0.001,
0.001,
0.001,
0.001,
0.001,
0.001,
0.001,
0.1,
0.1,
0.5,
0.5,
0.5,
0.5,
0.5,
0.5,
0.5,
0.5,
0.5,
0.5,
0.5,
0.5,
0.5,
0.5,
0.5,
0.5,
0.5,
0.5,
0.5,
0.5,
0.5,
0.5,
0.5
)


parset1 <- c(
0.001,
0.001,
0.001,
0.001,
0.001,
0.0001,
0.001,
0.001,
0.0001,
0.01,
0.1,
0.1,
0.1,
0.1,
0.3,
0.3,
0.3,
0.99,
0.99,
0.99,
0.1,
0.1,
0.1,
0.1,
0.1,
0.1,
0.2,
0.2,
0.2,
0.2,
0.5,
0.5,
0.999,
0.999
)




par_lower <- rep(0.0001,34)
par_upper <- rep(0.9999,34)



par_initial <- rbind(
  parset1,
c( 0.9999, 0.3644492, 1e-04, 0.2848697, 0.2074071, 0.7863725, 1e-04, 1e-04, 0.9999, 1e-04, 0.06113167, 0.835574, 1e-04, 0.603885, 0.6586143, 1e-04, 0.7283344, 0.3798453, 0.9999, 0.02925457, 0.9999, 0.5847419, 0.9999, 1e-04, 0.1586156, 1e-04, 0.2287690, 0.4015056, 0.1337562, 0.4083678, 0.1232729, 0.4748842, 0.694666, 1e-04 ),
c( 0.7576149, 1e-04, 1e-04, 0.1389601, 0.3459932, 0.912872, 0.4372959, 0.79994, 1e-04, 0.9999, 1e-04, 0.9999, 1e-04, 0.8509298, 0.9999, 1e-04, 0.791888, 1e-04, 0.4205591, 1e-04, 0.971159, 0.814358, 1e-04, 0.1733561, 0.510018, 0.9866322, 0.9999, 0.9999, 0.9149228, 0.5266742, 1e-04, 1e-04, 1e-04, 1e-04 ),
  parset
)


cat("\n\n")

# wrapper for model
aDGVM_wrapper <- function( parset )
{
# 	cat( "\n" )
	
# 	for ( i in 1:length(descr) )
# 	{
# 		cat( "double const ", descr[i], " = ", parset[i]/parfac[i], ";\n" )
# 	}
# 	cat( "\n" )

	# delete old output files
	system("rm OutputData/fynbos_fit/*.dat")
	
	parset_str <- NULL
	for ( i in 1:length(parset) ) parset_str <- paste(parset_str, parset[i])
	
	
	system("pgrep Model_all | wc > numproc")
	nn <- read.table("numproc")[1,1]

	t_start <- Sys.time()

	# run sinmulations
	for ( l in 1:length(xcoo) )
	{
		for ( k in 1:nruns )
		{
			system("pgrep Model_all | wc > numproc")
			nn <- read.table("numproc")[1,1]
# 			cat( nn, "\n")

			while ( nn>=nproc )
			{
# 				cat ("nn>nproc\n")
				Sys.sleep(3)
# 					nproc_ex <- 1
				system("pgrep Model_all | wc > numproc")
				nn <- read.table("numproc")[1,1]
# 				cat( nn, "\n")
			}
			
			system( paste( "./Model_all 100 251", xcoo[l], ycoo[l], "1", k, " 3 fynbos_fit", parset_str, " &" ) )
			#cat( paste( "./Model_all 100 251", xcoo[l], ycoo[l], "1", k, " 3 fynbos_fit", parset_str, " &" ), "\n" )
			nproc_ex <- 0
# 			cat ("nn<nproc, ", k, " ", l, "\n")
		}
	}
	
	
	output_filename <- "OutputData/fynbos_fit/YearlyData.250.3.1.dat"
	done         <- 0
	time_counter <- 0
	
	while ( done==0 )
	{
		if ( file.exists(output_filename ) )
		{
			yy <- read.table("OutputData/fynbos_fit/YearlyData.250.3.1.dat")
			if ( dim(yy)[1]== (length(xcoo)*nruns) ) done <- 1
		}
		Sys.sleep(3)
		time_counter <- time_counter+(3/60)
		if ( time_counter >=( 30 ) ) done <- 1
		cat(time_counter, " ")
	}
	
	cat("\n\nTime : ", time_counter, "\n", sep="" )
	
	vegtype <- 0
	
	# calculate sos
	sos <- 0
	if ( time_counter>=( 30 ) ) sos <- 99
	else
	{
		vegtype <- getVegType( yy )
		for ( i in 1:length(xcoo) )
		{
			dd <- subset( vegtype, vegtype[,1]==xcoo[i] )
			sos <- sos + sum(abs( sign(dd[,3]-rep(vtpy[i],dim(dd)[1])) ))
		}
	}
	
	
	cat("VAL ", parset, "\n" )
	cat("VTY ", vegtype, "\n" )


	cat("const double DEATH_PROB_FROST[NUM_TR_TYPES]  = { ", parset[1], ", ", parset[4], ", ", parset[7], " };\n", sep="" )
	cat("const double DEATH_PROB_CARBON[NUM_TR_TYPES] = { ", parset[2], ", ", parset[5], ", ", parset[8], " };\n", sep="" )
	cat("const double DEATH_PROB_COMP[NUM_TR_TYPES]   = { ", parset[3], ", ", parset[6], ", ", parset[9], " };\n", sep="" )
	cat("const double IGNITION_PROB  = ", parset[10], ";\n", sep="" )
	cat("const double IGNITION_PAR_2 = ", parset[11], ";\n", sep="" )
	cat("const double L_SAV_SAV = ",      parset[12], ";\n", sep="" )
	cat("const double L_SAV_FOR = ",      parset[13], ";\n", sep="" )
	cat("const double L_SAV_BBS = ",      parset[14], ";\n", sep="" )
	cat("const double L_FOR_SAV = ",      parset[15], ";\n", sep="" )
	cat("const double L_FOR_FOR = ",      parset[16], ";\n", sep="" )
	cat("const double L_FOR_BBS = ",      parset[17], ";\n", sep="" )
	cat("const double L_BBS_SAV = ",      parset[18], ";\n", sep="" )
	cat("const double L_BBS_FOR = ",      parset[19], ";\n", sep="" )
	cat("const double L_BBS_BBS = ",      parset[20], ";\n", sep="" )
	cat("const double L_C4G_SAV = ",      parset[21], ";\n", sep="" )
	cat("const double L_C4G_FOR = ",      parset[22], ";\n", sep="" )
	cat("const double L_C4G_BBS = ",      parset[23], ";\n", sep="" )
	cat("const double L_C3G_SAV = ",      parset[24], ";\n", sep="" )
	cat("const double L_C3G_FOR = ",      parset[25], ";\n", sep="" )
	cat("const double L_C3G_BBS = ",      parset[26], ";\n", sep="" )
	cat("const double L_C4G_C4G = ",      parset[27], ";\n", sep="" )
	cat("const double L_C3G_C3G = ",      parset[28], ";\n", sep="" )
	cat("const double L_SAV_C4G = ",      parset[29], ";\n", sep="" )
	cat("const double L_SAV_C3G = ",      parset[30], ";\n", sep="" )
	cat("const double L_FOR_C4G = ",      parset[31], ";\n", sep="" )
	cat("const double L_FOR_C3G = ",      parset[32], ";\n", sep="" )
	cat("const double L_BBS_C4G = ",      parset[33], ";\n", sep="" )
	cat("const double L_BBS_C3G = ",      parset[34], ";\n", sep="" )



	cat("OPT ------------------------------------------------------------------------------------- sos = ", round(sos,3), "\n\n" )
	
	
	return( sos )
}


# optim( parset, aDGVM_wrapper )
oo <- MyDEoptim( aDGVM_wrapper, lower=par_lower, upper=par_upper, initial=par_initial, control=list(itermax = 200, NP=50, refresh=1, digits=4) )



















