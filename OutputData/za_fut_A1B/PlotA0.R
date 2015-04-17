

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




getA0C3 <- function( d )
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
		x_index <- (vm[i,1]-min(vm[,1])+xres)/xres
		y_index <- (vm[i,2]-min(vm[,2])+yres)/yres
		mx[x_index, y_index] <- vm[i,53]
	}
	
	return( mx )
}


getA0C4 <- function( d )
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
		x_index <- (vm[i,1]-min(vm[,1])+xres)/xres
		y_index <- (vm[i,2]-min(vm[,2])+yres)/yres
		mx[x_index, y_index] <- vm[i,54]
	}
	
	return( mx )
}


# ---------------------------------------------------------------------------------
# --- Read input data -------------------------------------------------------------
# ---------------------------------------------------------------------------------


# cbind(1:339,(1:339+1761)) # get years

if ( read_data==1 )
{
	c3 <- getA0C3( read.table("YearlyData.1.3.0.dat") )
	c4 <- getA0C4( read.table("YearlyData.1.3.0.dat") )
}


# ---------------------------------------------------------------------------------
# --- Plot vegetation maps --------------------------------------------------------
# ---------------------------------------------------------------------------------


minx <- min(xx)-0.3
maxx <- max(xx)+0.3
miny <- min(yy)-0.3
maxy <- max(yy)+0.3


# ---------------------------------------------------------------------------------
# --- Plot C3:C4 ratio ------------------------------------------------------------
# ---------------------------------------------------------------------------------


# fi_brk <- c( -100, 1/c(30, 20, 15, 10, 6, 3, 2), 100)
gr_brk <- seq(-4,4,by=0.25)
gr_col <- c(colorRampPalette(c("red","white","blue"))(length(gr_brk)-1) )

a0_brk <- seq(0,14,by=0.5)
a0_col <- c(colorRampPalette(c("grey","blue"))(length(a0_brk)-1) )

graphics.off()
pdf( file="za_fynbos_a0.pdf", height=2.4, width=8.3 )
par( mfcol=c(1,3), mar=c(1,1,1,1), oma=c(0,2,2,0), bty="n", xaxt="n", yaxt="n" )

image( xx,yy,c3,   breaks=a0_brk, col=a0_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)

image( xx,yy,c4,   breaks=a0_brk, col=a0_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)

image( xx,yy,c3-c4,   breaks=gr_brk, col=gr_col, xlim=c(minx,maxx), ylim=c(miny,maxy) )
lines(bo[,1],bo[,2],type="l",lwd=0.5)



# mtext("With fire", 2, line=0, outer=T, at=0.25 )
# mtext("No fire",   2, line=0, outer=T, at=0.75 )

mtext("C3",    3, line=0, outer=T, at=0.16 )
mtext("C4",    3, line=0, outer=T, at=0.50 )
mtext("C3-C4", 3, line=0, outer=T, at=0.83 )


graphics.off()









