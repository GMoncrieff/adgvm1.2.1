

read_data <- 1

dimx <- 97
dimy <- 75

# coo <- read.table("south_africa_coo.dat")

xx <- read.table("za_x_coords.dat")[,1]
yy <- read.table("za_y_coords.dat")[,1]




# /* 11 */        yearly_output_data[10] = MyGrassPop.getBlMaxYearC4()*10.;   // t/ha
# /* 13 */        yearly_output_data[12] = MyGrassPop.getBlMaxYearC3()*10.;   // t/ha
# /* 16 */        yearly_output_data[15] = MyTreePop.getpCanopyYearMean(TR_SAV);
# /* 17 */        yearly_output_data[16] = MyTreePop.getpCanopyYearMean(TR_FOR);
# /* 18 */        yearly_output_data[17] = MyTreePop.getpCanopyYearMean(TR_BBS);
# /* 22 */        yearly_output_data[21] = MyGrassPop.getC34RatioYearMean();
# /* 26 */        yearly_output_data[25] = MyTreePop.getTreeNumType(TR_SAV);
# /* 27 */        yearly_output_data[26] = MyTreePop.getTreeNumType(TR_FOR);
# /* 28 */        yearly_output_data[27] = MyTreePop.getTreeNumType(TR_BBS);
# /* 34 */        yearly_output_data[33] = fire_num;
# /* 40 */        yearly_output_data[39] = d_params[0]; CO2
# /* 41 */        yearly_output_data[40] = d_params[1]; Temp
# /* 42 */        yearly_output_data[41] = d_params[2]; MAP
# /* 43 */        yearly_output_data[42] = d_params[3]; Seasonality
# /* 44 */        yearly_output_data[43] = d_params[4];
# /* 45 */        yearly_output_data[44] = d_params[5];
# /* 46 */        yearly_output_data[45] = d_params[6];
# /* 47 */        yearly_output_data[46] = d_params[7];
# /* 48 */        yearly_output_data[47] = d_params[8];
# /* 49 */        yearly_output_data[48] = d_params[9];
# /* 50 */        yearly_output_data[49] = MyTreePop.getActiveDays();
# /* 51 */        yearly_output_data[50] = MyGrassPop.getActiveDays();
# /* 52 */        yearly_output_data[51] = evapo_ref_sum;
# /* 53 */        yearly_output_data[52] = A0C3_mean;
# /* 54 */        yearly_output_data[53] = A0C4_mean;
# /* 55 */        yearly_output_data[54] = gs_C3_global;
# /* 56 */        yearly_output_data[55] = gs_C4_global;






# /* 11 */        yearly_output_data[10] = MyGrassPop.getBlMaxYearC4()*10.;   // t/ha
# /* 13 */        yearly_output_data[12] = MyGrassPop.getBlMaxYearC3()*10.;   // t/ha
# /* 16 */        yearly_output_data[15] = MyTreePop.getpCanopySavYearMean();
# /* 17 */        yearly_output_data[16] = MyTreePop.getpCanopyForYearMean();
# /* 21 */        yearly_output_data[20] = MyGrassPop.getC34RatioYearMean();
# /* 32 */        yearly_output_data[31] = fire_num;
# /* 38 */        yearly_output_data[37] = d_params[0]; CO2
# /* 39 */        yearly_output_data[38] = d_params[1]; Temp
# /* 40 */        yearly_output_data[39] = d_params[2]; MAP
# /* 41 */        yearly_output_data[40] = d_params[3]; Seasonality

library(sna)
library(sp)

bo <- read.table("za_border.dat")



getWater <- function( d )
{
	vm <- subset( d, point.in.polygon( d[,1],d[,2],bo[,1],bo[,2])==1 )
	
	xres <- abs(vm[1:(dim(vm)[1]-1),1]-vm[2:(dim(vm)[1]),1])
	xres <- min( xres[xres>0] )
	yres <- abs(vm[1:(dim(vm)[1]-1),2]-vm[2:(dim(vm)[1]),2])
	yres <- min( yres[yres>0] )
	
# 	mx <- matrix( -10, (max(vm[,1])-min(vm[,1])+xres)/xres,(max(vm[,2])-min(vm[,2])+yres)/yres )
	mx <- matrix( -1e15, dimx, dimy )
	
	for ( i in 1:dim(vm)[1] )
	{
		x_index <- (vm[i,1]-min(xx)+xres)/xres
		y_index <- (vm[i,2]-min(yy)+yres)/yres
		mx[x_index, y_index] <- (vm[i,7]-vm[i,52])/vm[i,7]/vm[i,3]
	}
	
	return( mx )
}




getC34R <- function( d )
{
	vm <- subset( d, point.in.polygon( d[,1],d[,2],bo[,1],bo[,2])==1 )
	
	xres <- abs(vm[1:(dim(vm)[1]-1),1]-vm[2:(dim(vm)[1]),1])
	xres <- min( xres[xres>0] )
	yres <- abs(vm[1:(dim(vm)[1]-1),2]-vm[2:(dim(vm)[1]),2])
	yres <- min( yres[yres>0] )
	
# 	mx <- matrix( -10, (max(vm[,1])-min(vm[,1])+xres)/xres,(max(vm[,2])-min(vm[,2])+yres)/yres )
	mx <- matrix( -10, dimx, dimy )
	
	for ( i in 1:dim(vm)[1] )
	{
		x_index <- (vm[i,1]-min(xx)+xres)/xres
		y_index <- (vm[i,2]-min(yy)+yres)/yres
		mx[x_index, y_index] <- vm[i,22]
	}
	
	return( mx )
}



getCSav <- function( d )
{
	vm <- subset( d, point.in.polygon( d[,1],d[,2],bo[,1],bo[,2])==1 )
	
	xres <- abs(vm[1:(dim(vm)[1]-1),1]-vm[2:(dim(vm)[1]),1])
	xres <- min( xres[xres>0] )
	yres <- abs(vm[1:(dim(vm)[1]-1),2]-vm[2:(dim(vm)[1]),2])
	yres <- min( yres[yres>0] )
	
# 	mx <- matrix( -10, (max(vm[,1])-min(vm[,1])+xres)/xres,(max(vm[,2])-min(vm[,2])+yres)/yres )
	mx <- matrix( -10, dimx, dimy )
	
	for ( i in 1:dim(vm)[1] )
	{
		x_index <- (vm[i,1]-min(xx)+xres)/xres
		y_index <- (vm[i,2]-min(yy)+yres)/yres
		mx[x_index, y_index] <- vm[i,16]
	}
	
	return( mx )
}

getCFor <- function( d )
{
	vm <- subset( d, point.in.polygon( d[,1],d[,2],bo[,1],bo[,2])==1 )
	
	xres <- abs(vm[1:(dim(vm)[1]-1),1]-vm[2:(dim(vm)[1]),1])
	xres <- min( xres[xres>0] )
	yres <- abs(vm[1:(dim(vm)[1]-1),2]-vm[2:(dim(vm)[1]),2])
	yres <- min( yres[yres>0] )
	
# 	mx <- matrix( -10, (max(vm[,1])-min(vm[,1])+xres)/xres,(max(vm[,2])-min(vm[,2])+yres)/yres )
	mx <- matrix( -10, dimx, dimy )
	
	for ( i in 1:dim(vm)[1] )
	{
		x_index <- (vm[i,1]-min(xx)+xres)/xres
		y_index <- (vm[i,2]-min(yy)+yres)/yres
		mx[x_index, y_index] <- vm[i,17]
	}
	
	return( mx )
}

getCBbs <- function( d )
{
	vm <- subset( d, point.in.polygon( d[,1],d[,2],bo[,1],bo[,2])==1 )
	
	xres <- abs(vm[1:(dim(vm)[1]-1),1]-vm[2:(dim(vm)[1]),1])
	xres <- min( xres[xres>0] )
	yres <- abs(vm[1:(dim(vm)[1]-1),2]-vm[2:(dim(vm)[1]),2])
	yres <- min( yres[yres>0] )
	
# 	mx <- matrix( -10, (max(vm[,1])-min(vm[,1])+xres)/xres,(max(vm[,2])-min(vm[,2])+yres)/yres )
	mx <- matrix( -10, dimx, dimy )
	
	for ( i in 1:dim(vm)[1] )
	{
		x_index <- (vm[i,1]-min(xx)+xres)/xres
		y_index <- (vm[i,2]-min(yy)+yres)/yres
		mx[x_index, y_index] <- vm[i,18]
	}
	
	return( mx )
}







getNPP <- function( d )
{
	vm <- subset( d, point.in.polygon( d[,1],d[,2],bo[,1],bo[,2])==1 )
	
	xres <- abs(vm[1:(dim(vm)[1]-1),1]-vm[2:(dim(vm)[1]),1])
	xres <- min( xres[xres>0] )
	yres <- abs(vm[1:(dim(vm)[1]-1),2]-vm[2:(dim(vm)[1]),2])
	yres <- min( yres[yres>0] )
	
# 	mx <- matrix( -10, (max(vm[,1])-min(vm[,1])+xres)/xres,(max(vm[,2])-min(vm[,2])+yres)/yres )
	mx <- matrix( -10, dimx, dimy )
	
	for ( i in 1:dim(vm)[1] )
	{
		x_index <- (vm[i,1]-min(xx)+xres)/xres
		y_index <- (vm[i,2]-min(yy)+yres)/yres
		mx[x_index, y_index] <- vm[i,30]
	}
	
	return( mx )
}


getFint <- function( d )
{
	vm <- subset( d, point.in.polygon( d[,1],d[,2],bo[,1],bo[,2])==1 )
	
	xres <- abs(vm[1:(dim(vm)[1]-1),1]-vm[2:(dim(vm)[1]),1])
	xres <- min( xres[xres>0] )
	yres <- abs(vm[1:(dim(vm)[1]-1),2]-vm[2:(dim(vm)[1]),2])
	yres <- min( yres[yres>0] )
	
# 	mx <- matrix( -10, (max(vm[,1])-min(vm[,1])+xres)/xres,(max(vm[,2])-min(vm[,2])+yres)/yres )
	mx <- matrix( -10, dimx, dimy )
	
	for ( i in 1:dim(vm)[1] )
	{
		x_index <- (vm[i,1]-min(xx)+xres)/xres
		y_index <- (vm[i,2]-min(yy)+yres)/yres
		mx[x_index, y_index] <- vm[i,35]
		if ( vm[i,35]=="Nan" ) mx[x_index, y_index] <- 0
	}
	
	return( mx )
}


getFnum <- function( d )
{
	vm <- subset( d, point.in.polygon( d[,1],d[,2],bo[,1],bo[,2])==1 )
	
	xres <- abs(vm[1:(dim(vm)[1]-1),1]-vm[2:(dim(vm)[1]),1])
	xres <- min( xres[xres>0] )
	yres <- abs(vm[1:(dim(vm)[1]-1),2]-vm[2:(dim(vm)[1]),2])
	yres <- min( yres[yres>0] )
	
# 	mx <- matrix( -10, (max(vm[,1])-min(vm[,1])+xres)/xres,(max(vm[,2])-min(vm[,2])+yres)/yres )
	mx <- matrix( -10, dimx, dimy )
	
	for ( i in 1:dim(vm)[1] )
	{
		x_index <- (vm[i,1]-min(xx)+xres)/xres
		y_index <- (vm[i,2]-min(yy)+yres)/yres
		mx[x_index, y_index] <- vm[i,34]/(vm[i,3]-30)
	}
	
	return( mx )
}


getXCoo <- function( d )
{
# 	vm <- d
	vm <- subset( d, point.in.polygon( d[,1],d[,2],bo[,1],bo[,2])==1 )
	
	xres <- abs(vm[1:(dim(vm)[1]-1),1]-vm[2:(dim(vm)[1]),1])
	xres <- min( xres[xres>0] )
	yres <- abs(vm[1:(dim(vm)[1]-1),2]-vm[2:(dim(vm)[1]),2])
	yres <- min( yres[yres>0] )
	
# 	mx <- matrix( -10, (max(vm[,1])-min(vm[,1])+xres)/xres,(max(vm[,2])-min(vm[,2])+yres)/yres )
	mx <- matrix( -10, dimx, dimy )
	
	for ( i in 1:dim(vm)[1] )
	{
		x_index <- (vm[i,1]-min(xx)+xres)/xres
		y_index <- (vm[i,2]-min(yy)+yres)/yres
		mx[x_index, y_index] <- vm[i,1]
	}
	
	return( mx )
}


getYCoo <- function( d )
{
	vm <- subset( d, point.in.polygon( d[,1],d[,2],bo[,1],bo[,2])==1 )
	
	xres <- abs(vm[1:(dim(vm)[1]-1),1]-vm[2:(dim(vm)[1]),1])
	xres <- min( xres[xres>0] )
	yres <- abs(vm[1:(dim(vm)[1]-1),2]-vm[2:(dim(vm)[1]),2])
	yres <- min( yres[yres>0] )
	
# 	mx <- matrix( -10, (max(vm[,1])-min(vm[,1])+xres)/xres,(max(vm[,2])-min(vm[,2])+yres)/yres )
	mx <- matrix( -10, dimx, dimy )
	
	for ( i in 1:dim(vm)[1] )
	{
		x_index <- (vm[i,1]-min(xx)+xres)/xres
		y_index <- (vm[i,2]-min(yy)+yres)/yres
		mx[x_index, y_index] <- vm[i,2]
	}
	
	return( mx )
}


getVegType <- function( d )
{
	d <- subset( d, point.in.polygon( d[,1],d[,2],bo[,1],bo[,2])==1 )
	
	# vectors for classification of cells
	v.c3s <- rep( 0, dim(d)[1] )
	v.c4s <- rep( 0, dim(d)[1] )
	v.for <- rep( 0, dim(d)[1] )
	v.c3g <- rep( 0, dim(d)[1] )
	v.c4g <- rep( 0, dim(d)[1] )
	v.wdl <- rep( 0, dim(d)[1] )
	v.des <- rep( 0, dim(d)[1] )
	v.fyn <- rep( 0, dim(d)[1] )
	
	d <- d[,1:47]    ## added columns to YearlyData files, this truncates to old size
	
	# do vegetation classification for year 300
	for ( i in 1:dim(d)[1] )
	{
		if      ( d[i,16]+d[i,17]+d[i,18] < 0.1 & d[i,11]+d[i,13] < 1.5 ) v.des[i] <- v.des[i]+1    # desert
		else if ( d[i,16]+d[i,17]+d[i,18] < 0.1 & d[i,11]+d[i,13] > 1.5 )       # grasslands
		{
			if ( d[i,22]<0.5 )                            v.c4g[i] <- v.c4g[i]+1
			else                                          v.c3g[i] <- v.c3g[i]+1
		}
		else if ( d[i,16]+d[i,17] > 0.8 )                 v.for[i] <- v.for[i]+1  # forest
# 		else if ( d[i,18] > 0.1 & d[i,16]+d[i,17] < 0.1 ) v.fyn[i] <- v.fyn[i]+1  # fynbos
# 		else if ( d[i,18] > 0.4 | d[i,28]>2000 )          v.fyn[i] <- v.fyn[i]+1  # fynbos
		else if ( d[i,28]>2500 )          v.fyn[i] <- v.fyn[i]+1  # fynbos
		else
		{
# 			if   ( d[i,32]> 0 & d[i,22]<0.5) v.c4s[i] <- v.c4s[i]+1
			if      ( d[i,16]>d[i,17] & d[i,22]<0.5 ) v.c4s[i] <- v.c4s[i]+1
			else if ( d[i,16]>d[i,17] & d[i,22]>0.5 ) v.c3s[i] <- v.c3s[i]+1
			else                                      v.wdl[i] <- v.wdl[i]+1
		}
	}
	vm <- cbind( d, v.des, v.c4g, v.c3g, v.c4s, v.wdl, v.for, v.c3s, v.fyn )
	
	xres <- abs(vm[1:(dim(vm)[1]-1),1]-vm[2:(dim(vm)[1]),1])
	xres <- min( xres[xres>0] )
	yres <- abs(vm[1:(dim(vm)[1]-1),2]-vm[2:(dim(vm)[1]),2])
	yres <- min( yres[yres>0] )
	
# 	mx <- matrix( 0, (max(vm[,1])-min(vm[,1])+xres)/xres,(max(vm[,2])-min(vm[,2])+yres)/yres )
	mx <- matrix( 0, dimx, dimy )
	
	for ( i in 1:dim(vm)[1] )
	{
		x_index <- (vm[i,1]-min(xx)+xres)/xres
		y_index <- (vm[i,2]-min(yy)+yres)/yres
# 		x_index <- (vm[i,1]-min(vm[,1])+xres)/xres
# 		y_index <- (vm[i,2]-min(vm[,2])+yres)/yres
		mx[x_index, y_index] <- which(vm[i,48:55]==1)
# 		print(c(x_index,y_index))
	}
	
	return( mx )
}


# ---------------------------------------------------------------------------------
# --- Read input data -------------------------------------------------------------
# ---------------------------------------------------------------------------------


# cbind(1:339,(1:339+1761)) # get years

if ( read_data==1 )
{
	vmpst.0 <- getVegType( read.table("YearlyData.139.3.0.dat") )
	vmcur.0 <- getVegType( read.table("YearlyData.251.3.0.dat") )
	vmfut.0 <- getVegType( read.table("YearlyData.339.3.0.dat") )

	vmpst.1 <- getVegType( read.table("YearlyData.139.3.1.dat") )
	vmcur.1 <- getVegType( read.table("YearlyData.251.3.1.dat") )
	vmfut.1 <- getVegType( read.table("YearlyData.339.3.1.dat") )

	vm_30.1 <- getVegType( read.table("YearlyData.269.3.1.dat") )
	vm_50.1 <- getVegType( read.table("YearlyData.289.3.1.dat") )


	vmsav.1 <- getCSav( read.table("YearlyData.251.3.1.dat") )
	vmfor.1 <- getCFor( read.table("YearlyData.251.3.1.dat") )
	vmbbs.1 <- getCBbs( read.table("YearlyData.251.3.1.dat") )

	vmsav.0 <- getCSav( read.table("YearlyData.251.3.0.dat") )
	vmfor.0 <- getCFor( read.table("YearlyData.251.3.0.dat") )
	vmbbs.0 <- getCBbs( read.table("YearlyData.251.3.0.dat") )

	grpst.1 <- getC34R( read.table("YearlyData.139.3.1.dat") )
	grcur.1 <- getC34R( read.table("YearlyData.251.3.1.dat") )
	grfut.1 <- getC34R( read.table("YearlyData.339.3.1.dat") )

	grpst.0 <- getC34R( read.table("YearlyData.139.3.0.dat") )
	grcur.0 <- getC34R( read.table("YearlyData.251.3.0.dat") )
	grfut.0 <- getC34R( read.table("YearlyData.339.3.0.dat") )


# 	xcoo <- getXCoo( read.table("YearlyData.139.3.1.dat") )
# 	ycoo <- getYCoo( read.table("YearlyData.139.3.1.dat") )

# 	xx <- unique(sort(xcoo))
# 	xx <- xx[2:length(xx)]
# 	yy <- unique(sort(ycoo))
# 	yy <- yy[1:(length(yy)-1)]


	fipst.1 <- getFint( read.table("YearlyData.139.3.1.dat") )
	ficur.1 <- getFint( read.table("YearlyData.251.3.1.dat") )
	fifut.1 <- getFint( read.table("YearlyData.339.3.1.dat") )

	fnpst.1 <- getFnum( read.table("YearlyData.139.3.1.dat") )
	fncur.1 <- getFnum( read.table("YearlyData.251.3.1.dat") )
	fnfut.1 <- getFnum( read.table("YearlyData.339.3.1.dat") )

	pppst.0 <- getNPP( read.table("YearlyData.139.3.0.dat") )
	ppcur.0 <- getNPP( read.table("YearlyData.251.3.0.dat") )
	ppfut.0 <- getNPP( read.table("YearlyData.339.3.0.dat") )

	pppst.1 <- getNPP( read.table("YearlyData.139.3.1.dat") )
	ppcur.1 <- getNPP( read.table("YearlyData.251.3.1.dat") )
	ppfut.1 <- getNPP( read.table("YearlyData.339.3.1.dat") )


	wapst.0 <- getWater( read.table("YearlyData.139.3.0.dat") )
	wacur.0 <- getWater( read.table("YearlyData.251.3.0.dat") )
	wafut.0 <- getWater( read.table("YearlyData.339.3.0.dat") )

	wapst.1 <- getWater( read.table("YearlyData.139.3.1.dat") )
	wacur.1 <- getWater( read.table("YearlyData.251.3.1.dat") )
	wafut.1 <- getWater( read.table("YearlyData.339.3.1.dat") )


}


# ---------------------------------------------------------------------------------
# --- Plot vegetation maps --------------------------------------------------------
# ---------------------------------------------------------------------------------

ccc <- terrain.colors(6)
vt_col <- c( "white", "grey", rgb(186/255, 161/255, 7/255), rgb(163/255, 71/255, 0/255), rgb(173/255, 21/255, 20/255), rgb(1/255, 84/255, 160/255), rgb(73/255, 14/255, 82/255), "black" )
vt_brk <- c( -100, 0,  1,       2,        3,          4,             5,           6, 7, 8 )
colors12<-matrix(c(1.000,0.750,0.500,1.000,0.500,0.000 ,1.000,1.000,0.600 ,1.000,1.000,0.200 ,0.700,1.000,0.550, 
                 0.200,1.000,0.000,0.650,0.930,1.000 ,0.100,0.700,1.000 ,0.800,0.750,1.000 ,0.400,0.300,1.000, 
                 1.000,0.600,0.750,0.900,0.100,0.200 ),ncol=3,nrow=12,byrow=T)
col12<-rgb(colors12)

#          ocean    v.c4g,    v.sav,    v.wdl,   v.for,    v.c3g,   v.des
# 	vm <- cbind( d,   v.des,    v.c4g,   v.c3g,   v.c4s,    v.wdl,   v.for,   v.c3s,   v.fyn )
# vt_col<-c("white",grey(0.85),col12[1],col12[3],col12[11],col12[7],col12[9],col12[5],"black")#col12[10])
vt_col<-c("white",grey(0.85),col12[1],col12[3],col12[11],col12[7],col12[9],col12[5], "brown")#col12[10])


minx <- min(xx)-0.3
maxx <- max(xx)+0.3
miny <- min(yy)-0.3
maxy <- max(yy)+0.3


graphics.off()
pdf( file="za_fynbos_maps.pdf", height=4.8, width=8.3 )
par( mfcol=c(2,3), mar=c(1,1,1,1), oma=c(0,2,2,0), bty="n", xaxt="n", yaxt="n" )


image( xx,yy,vmpst.0,   breaks=vt_brk, col=vt_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
image( xx,yy,vmpst.1,   breaks=vt_brk, col=vt_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)

image( xx,yy,vmcur.0,   breaks=vt_brk, col=vt_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
image( xx,yy,vmcur.1,   breaks=vt_brk, col=vt_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)

image( xx,yy,vmfut.0,   breaks=vt_brk, col=vt_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
image( xx,yy,vmfut.1,   breaks=vt_brk, col=vt_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)

# 
# image( vmpst.0,   breaks=vt_brk, col=vt_col )
# image( vmpst.1,   breaks=vt_brk, col=vt_col )
# 
# image( vmcur.0,   breaks=vt_brk, col=vt_col )
# image( vmcur.1,   breaks=vt_brk, col=vt_col )
# 
# image( vmfut.0,   breaks=vt_brk, col=vt_col )
# image( vmfut.1,   breaks=vt_brk, col=vt_col )


mtext("With fire ", 2, line=0, outer=T, at=0.25 )
mtext("Fire suppression", 2, line=0, outer=T, at=0.75 )

mtext("1900", 3, line=0, outer=T, at=0.16 )
mtext("2012", 3, line=0, outer=T, at=0.50 )
mtext("2100", 3, line=0, outer=T, at=0.83 )

graphics.off()



graphics.off()
pdf( file="za_fynbos_series.pdf", height=4.8, width=5.5 )
par( mfcol=c(2,2), mar=c(1,1,1,1), oma=c(0,0,0,0), bty="n", xaxt="n", yaxt="n" )

image( xx,yy,vmcur.1,   breaks=vt_brk, col=vt_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
text(18,-23.,"2012",cex=1.3)
image( xx,yy,vm_50.1,   breaks=vt_brk, col=vt_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
text(18,-23.,"2050",cex=1.3)

image( xx,yy,vm_30.1,   breaks=vt_brk, col=vt_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
text(18,-23.,"2030",cex=1.3)
image( xx,yy,vmfut.1,   breaks=vt_brk, col=vt_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
text(18,-23.,"2100",cex=1.3)

graphics.off()



# ---------------------------------------------------------------------------------
# --- Plot C3:C4 ratio ------------------------------------------------------------
# ---------------------------------------------------------------------------------


# fi_brk <- c( -100, 1/c(30, 20, 15, 10, 6, 3, 2), 100)
gr_brk <- seq(0,1,by=0.1)
gr_col <- c(colorRampPalette(c("red","blue"))(length(gr_brk)-1) )


graphics.off()
pdf( file="za_fynbos_c34.pdf", height=4.8, width=8.3 )
par( mfcol=c(2,3), mar=c(1,1,1,1), oma=c(0,2,2,0), bty="n", xaxt="n", yaxt="n" )

image( xx,yy,grpst.0,   breaks=gr_brk, col=gr_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
image( xx,yy,grpst.1,   breaks=gr_brk, col=gr_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)

image( xx,yy,grcur.0,   breaks=gr_brk, col=gr_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
image( xx,yy,grcur.1,   breaks=gr_brk, col=gr_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)

image( xx,yy,grfut.0,   breaks=gr_brk, col=gr_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
image( xx,yy,grfut.1,   breaks=gr_brk, col=gr_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)


mtext("With fire", 2, line=0, outer=T, at=0.25 )
mtext("No fire",   2, line=0, outer=T, at=0.75 )

mtext("1900", 3, line=0, outer=T, at=0.16 )
mtext("2012", 3, line=0, outer=T, at=0.50 )
mtext("2100", 3, line=0, outer=T, at=0.83 )


graphics.off()
pdf( file="gr_legend.pdf", width=7.5, height=.23 )
par( mfcol=c(1,1), mar=c(0,0,0,0), bty="n", xaxt="n", yaxt="n" )
plot( -100, -100, xlim=c(0,500), ylim=c(0,10) )

for ( i in 2:length(gr_col) )
{
	x_step <- 60
	x_val  <- (i-2)*x_step+1
	rect( x_val, 1, x_val+13, 9, col=gr_col[i], border=gr_col[i])
	if (i<length(gr_col)) text( x_val+9,4, pos=4, cex=1, paste( round(gr_brk[i]*100,1), "-", round(gr_brk[i+1]*100,1), sep="" ))
}
text( x_val+9,4, pos=4, cex=1, paste( ">",round(gr_brk[i]*100,1), "", sep="" ))


graphics.off()





# ---------------------------------------------------------------------------------
# --- Plot tree cover maps --------------------------------------------------------
# ---------------------------------------------------------------------------------


# fi_brk <- c( -100, 1/c(30, 20, 15, 10, 6, 3, 2), 100)
vm_brk <- seq(0,1,0.01)
vm_col <- c(colorRampPalette(c("grey","red"))(length(vm_brk)-1) )


graphics.off()
pdf( file="za_fynbos_cover.pdf", height=4.8, width=8.3 )
par( mfcol=c(2,3), mar=c(1,1,1,1), oma=c(0,2,2,0), bty="n", xaxt="n", yaxt="n" )

image( xx,yy,vmsav.0,   breaks=vm_brk, col=vm_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
image( xx,yy,vmsav.1,   breaks=vm_brk, col=vm_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)

image( xx,yy,vmfor.0,   breaks=vm_brk, col=vm_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
image( xx,yy,vmfor.1,   breaks=vm_brk, col=vm_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)

image( xx,yy,vmbbs.0,   breaks=vm_brk, col=vm_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
image( xx,yy,vmbbs.1,   breaks=vm_brk, col=vm_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)


mtext("With fire", 2, line=0, outer=T, at=0.25 )
mtext("No fire",   2, line=0, outer=T, at=0.75 )

mtext("Savanna", 3, line=0, outer=T, at=0.16 )
mtext("Forest", 3, line=0, outer=T, at=0.50 )
mtext("Burning bush", 3, line=0, outer=T, at=0.83 )


graphics.off()








# ---------------------------------------------------------------------------------
# --- Plot water maps -------------------------------------------------------------
# ---------------------------------------------------------------------------------


# fi_brk <- c( -100, 1/c(30, 20, 15, 10, 6, 3, 2), 100)
# wa_brk <- seq(-1550, 1550, 100)
wa_brk <- c( -1000, -1, -0.2, -0.1, seq(-.1, .1, length=50), 0.1, 0.2, 1, 1000 )
wa_col <- c(colorRampPalette(c("red","grey","blue"))(length(wa_brk)-1) )


graphics.off()
pdf( file="za_fynbos_water.pdf", height=4.8, width=8.3 )
par( mfcol=c(2,3), mar=c(1,1,1,1), oma=c(0,2,2,0), bty="n", xaxt="n", yaxt="n" )

image( xx,yy,wapst.0,   breaks=wa_brk, col=wa_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
image( xx,yy,wapst.1,   breaks=wa_brk, col=wa_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)

image( xx,yy,wacur.0,   breaks=wa_brk, col=wa_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
image( xx,yy,wacur.1,   breaks=wa_brk, col=wa_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)

image( xx,yy,wafut.0,   breaks=wa_brk, col=wa_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
image( xx,yy,wafut.1,   breaks=wa_brk, col=wa_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)


mtext("With fire ", 2, line=0, outer=T, at=0.25 )
mtext("Fire suppression", 2, line=0, outer=T, at=0.75 )

mtext("1900", 3, line=0, outer=T, at=0.16 )
mtext("2012", 3, line=0, outer=T, at=0.50 )
mtext("2100", 3, line=0, outer=T, at=0.83 )

graphics.off()















# ---------------------------------------------------------------------------------
# --- Plot fire maps --------------------------------------------------------------
# ---------------------------------------------------------------------------------


# fi_brk <- c( -100, 1/c(30, 20, 15, 10, 6, 3, 2), 100)
fi_brk <- c( -0.1, 10, 200, 350, 500, 750, 1000, 2000, 100000 )
fi_col <- c(colorRampPalette(c(grey(0.1),rgb(214/255,96/255,77/255),rgb(254/255,224/255,139/255)))(length(fi_brk)-1) )

fn_brk <- c( -1, c(0.0000001, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45), 100)
fn_brk <- c( -1, 1/c(100000, 10, 5, 4, 3, 2.5, 2), 100)
fn_col <- c(colorRampPalette(c(grey(0.1),rgb(214/255,96/255,77/255),rgb(254/255,224/255,139/255)))(length(fn_brk)-1) )

graphics.off()
pdf( file="za_fynbos_fire.pdf", height=4.8, width=8.3 )
par( mfcol=c(2,3), mar=c(1,1,1,1), oma=c(0,2,2,0), bty="n", xaxt="n", yaxt="n" )

image( xx,yy,fipst.1,   breaks=fi_brk, col=fi_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
image( xx,yy,fnpst.1,   breaks=fn_brk, col=fn_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)

image( xx,yy,ficur.1,   breaks=fi_brk, col=fi_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
image( xx,yy,fncur.1,   breaks=fn_brk, col=fn_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)

image( xx,yy,fifut.1,   breaks=fi_brk, col=fi_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
image( xx,yy,fnfut.1,   breaks=fn_brk, col=fn_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)


mtext("Fire frequency", 2, line=0, outer=T, at=0.25 )
mtext("Fire intensity", 2, line=0, outer=T, at=0.75 )

mtext("1900", 3, line=0, outer=T, at=0.16 )
mtext("2012", 3, line=0, outer=T, at=0.50 )
mtext("2100", 3, line=0, outer=T, at=0.83 )


graphics.off()
pdf( file="fi_legend.pdf", width=7.5, height=.23 )
par( mfcol=c(1,1), mar=c(0,0,0,0), bty="n", xaxt="n", yaxt="n" )
plot( -100, -100, xlim=c(0,500), ylim=c(0,10) )

fi_brk <- fi_brk/1000
for ( i in 2:length(fi_col) )
{
	x_step <- 60
	x_val  <- (i-2)*x_step+1
	rect( x_val, 1, x_val+13, 9, col=fi_col[i], border=fi_col[i])
	if (i<length(fi_col)) text( x_val+9,4, pos=4, cex=1, paste( round(fi_brk[i],1), "-", round(fi_brk[i+1],1), sep="" ))
}
text( x_val+9,4, pos=4, cex=1, paste( ">",round(fi_brk[i],1), "MJ/s/m", sep="" ))


graphics.off()


graphics.off()
pdf( file="fn_legend.pdf", width=7.5, height=.23 )
par( mfcol=c(1,1), mar=c(0,0,0,0), bty="n", xaxt="n", yaxt="n" )
plot( -100, -100, xlim=c(0,500), ylim=c(0,10) )

for ( i in 2:length(fn_col) )
{
	x_step <- 60
	x_val  <- (i-2)*x_step+1
	rect( x_val, 1, x_val+13, 9, col=fn_col[i], border=fn_col[i])
	if (i<length(fn_col)) text( x_val+9,4, pos=4, cex=1, paste( round(fn_brk[i],2), "-", round(fn_brk[i+1],2), sep="" ))
}
text( x_val+9,4, pos=4, cex=1, paste( ">",round(fn_brk[i],2), " fires/year", sep="" ))

graphics.off()


# ---------------------------------------------------------------------------------
# --- Plot npp maps ---------------------------------------------------------------
# ---------------------------------------------------------------------------------


pp_brk <- c( -1, seq(0,0.8,0.1))
# pp_col <- c(colorRampPalette(c(grey(0.1),rgb(214/255,96/255,77/255),rgb(254/255,224/255,139/255)))(length(fi_brk)-1) )
pp_col <- terrain.colors(length(fi_brk))

graphics.off()
pdf( file="za_fynbos_npp.pdf", height=4.8, width=8.3 )
par( mfcol=c(2,3), mar=c(1,1,1,1), oma=c(0,2,2,0), bty="n", xaxt="n", yaxt="n" )

image( xx,yy,pppst.0,   breaks=pp_brk, col=pp_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
image( xx,yy,pppst.1,   breaks=pp_brk, col=pp_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)

image( xx,yy,ppcur.0,   breaks=pp_brk, col=pp_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
image( xx,yy,ppcur.1,   breaks=pp_brk, col=pp_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)

image( xx,yy,ppfut.0,   breaks=pp_brk, col=pp_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)
image( xx,yy,ppfut.1,   breaks=pp_brk, col=pp_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)


mtext("With fire ", 2, line=0, outer=T, at=0.25 )
mtext("Fire suppression", 2, line=0, outer=T, at=0.75 )

mtext("1900", 3, line=0, outer=T, at=0.16 )
mtext("2012", 3, line=0, outer=T, at=0.50 )
mtext("2100", 3, line=0, outer=T, at=0.83 )

graphics.off()


graphics.off()
pdf( file="npp_legend.pdf", width=7.5, height=.23 )
par( mfcol=c(1,1), mar=c(0,0,0,0), bty="n", xaxt="n", yaxt="n" )
plot( -100, -100, xlim=c(0,500), ylim=c(0,10) )

for ( i in 2:length(pp_col) )
{
	x_step <- 60
	x_val  <- (i-2)*x_step+1
	rect( x_val, 1, x_val+13, 9, col=pp_col[i], border=pp_col[i])
	if (i<length(pp_col)) text( x_val+9,4, pos=4, cex=1, paste( round(pp_brk[i],1), "-", round(pp_brk[i+1],1), sep="" ))
}
text( x_val+9,4, pos=4, cex=1, paste( ">",round(pp_brk[i],1), "kg C/m²", sep="" ))

graphics.off()




# --------------------------------------------------------------------------------------
# --- Get transition matrices ----------------------------------------------------------
# --------------------------------------------------------------------------------------


getTransMat <- function( m1, m2 )
{
  mr <- matrix( 0, 8, 8 )  # matrix that stores transitions between the 7 biome types

  for ( i in 1:dim(m1)[1] )
  {
    for ( j in 1:dim(m1)[2] )
    {
      t1 <- m1[i,j]
      t2 <- m2[i,j]
      mr[t1,t2] <- mr[t1,t2]+1
    }
  }

  cat("\n")
  cat("from, to  des   c4g   c3g   c4s   wdl   for   c3s   fyn\n")

  cat("des    ", format( mr[1,1:8], width=5), "\n")
  cat("c4g    ", format( mr[2,1:8], width=5), "\n")
  cat("c3g    ", format( mr[3,1:8], width=5), "\n")
  cat("c4s    ", format( mr[4,1:8], width=5), "\n")
  cat("wdl    ", format( mr[5,1:8], width=5), "\n")
  cat("for    ", format( mr[6,1:8], width=5), "\n")
  cat("c3s    ", format( mr[7,1:8], width=5), "\n")
  cat("fyn    ", format( mr[8,1:8], width=5), "\n")
  cat("\n")

  return(mr)
}


m.p_xx <- getTransMat( vmpst.1, vmpst.1 )
m.p_fs <- getTransMat( vmpst.1, vmpst.0 )

m.c_xx <- getTransMat( vmpst.1, vmcur.1 )
m.c_fs <- getTransMat( vmcur.1, vmcur.0 )

m.f_xx <- getTransMat( vmcur.1, vmfut.1 )
m.f_fs <- getTransMat( vmfut.1, vmfut.0 )



# --------------------------------------------------------------------------------------
# --- Plot transition graphs -----------------------------------------------------------
# --------------------------------------------------------------------------------------



gplot.layout.myla_var <- function( d, layout.par ) return( matrix(c(2,0,  3.8,1,  0.2,1,  3.8,3,  2.,2,  2,4,  0.2,3 )/9+0.06,7,2, byrow=T))



plotGraphVar <- function(m)
{
# 	biomes <- c("Desert", "C4-grl", "C3-grl", "C4-Sav", "Woodl", "Forest", "C3-Sav", "Fynbos")
	labs <- rep(0,7)
	pvec <- c(8,2:7)
	for ( i in 1:length(pvec) ) labs[i] <- paste(round(100*sum(m[,pvec[i]])/sum(m), 1))
print(labs)
	
	gplot((m[c(8,2:7),c(8,2:7)]), diag=F, vertex.cex=7.5, vertex.sides=40, vertex.rot=45, displaylabels=TRUE, label=labs, label.pos=5, arrowhead.cex=1.5, loop.cex=1., vertex.col=vt_col[c(8,2:7)+1], edge.lwd=m[c(8,2:7),c(8,2:7)]/20, thresh=15, mode="myla_var", jitter=F, label.cex=1.6, new=F )#, label.col=c("black", "black", "white", "white", "white", "white", "white"))
}

# plot(-100,-100, xlim=c(0,maxx), ylim=c(0,maxy))
# plotGraphVar( m.p_fs )



graphics.off()
pdf( file="za_fynbos_graphs.pdf", height=5.6, width=8.3 )
par( mfcol=c(2,3), mar=c(1,1,1,1), oma=c(0,2,2,0), bty="n", xaxt="n", yaxt="n" )

maxx <- 0.55
maxy <- 0.55

plot(-100,-100, xlim=c(0,maxx), ylim=c(0,maxy))
plotGraphVar( m.p_fs )
plot(-100,-100, xlim=c(0,maxx), ylim=c(0,maxy))
plotGraphVar( m.p_xx )

plot(-100,-100, xlim=c(0,maxx), ylim=c(0,maxy))
plotGraphVar( m.c_fs )
plot(-100,-100, xlim=c(0,maxx), ylim=c(0,maxy))
plotGraphVar( m.c_xx )

plot(-100,-100, xlim=c(0,maxx), ylim=c(0,maxy))
plotGraphVar( m.f_fs )
plot(-100,-100, xlim=c(0,maxx), ylim=c(0,maxy))
plotGraphVar( m.f_xx )

mtext("With fire ", 2, line=0, outer=T, at=0.25 )
mtext("Fire suppression", 2, line=0, outer=T, at=0.75 )

mtext("1900", 3, line=0, outer=T, at=0.16 )
mtext("2012", 3, line=0, outer=T, at=0.50 )
mtext("2100", 3, line=0, outer=T, at=0.83 )

graphics.off()




# system("pdflatex figures.tex")


