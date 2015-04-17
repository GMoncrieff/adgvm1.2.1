library(sna)
library(sp)

bo <- read.table("za_border.dat")

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
	
# 	d <- d[,1:47]    ## added columns to YearlyData files, this truncates to old size
	
	vv <- rep( 0, nrow(d) )

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
		vv[i] <- which( c(v.des[i], v.c4g[i], v.c3g[i], v.c4s[i], v.wdl[i], v.for[i], v.c3s[i], v.fyn[i])==1 )

	}
	
	vm <- cbind( d, vv )
	
	return( vm )
}



v <- getVegType( read.table("YearlyData.251.3.1.dat") )
v[,7] <- v[,7]/v[,3]

colors12<-matrix(c(1.000,0.750,0.500,1.000,0.500,0.000 ,1.000,1.000,0.600 ,1.000,1.000,0.200 ,0.700,1.000,0.550, 
                 0.200,1.000,0.000,0.650,0.930,1.000 ,0.100,0.700,1.000 ,0.800,0.750,1.000 ,0.400,0.300,1.000, 
                 1.000,0.600,0.750,0.900,0.100,0.200 ),ncol=3,nrow=12,byrow=T)
col12<-rgb(colors12)



#         v.des      v.c4g    v.c3g    v.c4s     v.wdl    v.for    v.c3s     v.fyn
vt_col<-c(grey(0.85),col12[1],col12[3],col12[11],col12[7],col12[9],col12[5], "brown" )


ind <- 7

graphics.off()
pdf( file="treenumbers.pdf", height=5, width=14 )

par( mar=c(4,4,1,1), mfcol=c(1,3) )

plot( -100000, -100000, xlim=c( min(v[,ind]), max(v[,ind]) ), ylim=c(0, max(v[,26:28]) ) )
for ( i in 1:8 ) 
{
	g <- subset( v, v[,57]==i )
	points( g[,ind], g[,26], col=vt_col[i], pch=18 )
}

plot( -100000, -100000, xlim=c( min(v[,ind]), max(v[,ind]) ), ylim=c(0, max(v[,26:28]) ) )
for ( i in 1:8 ) 
{
	g <- subset( v, v[,57]==i )
	points( g[,ind], g[,27], col=vt_col[i], pch=18 )
}

plot( -100000, -100000, xlim=c( min(v[,ind]), max(v[,ind]) ), ylim=c(0, max(v[,26:28]) ) )
for ( i in 1:8 ) 
{
	g <- subset( v, v[,57]==i )
	points( g[,ind], g[,28], col=vt_col[i], pch=18 )
}

graphics.off()




graphics.off()
pdf( file="biomes.pdf", height=5, width=5 )

par( mar=c(4,4,1,1), mfcol=c(1,1) )

plot( -100000, -100000, xlim=c( min(v[,ind]), max(v[,ind]) ), ylim=c(1, 8 ) )
for ( i in 1:8 ) 
{
	g <- subset( v, v[,57]==i )
	points( g[,ind], g[,57], col=vt_col[i], pch=18 )
}


graphics.off()












