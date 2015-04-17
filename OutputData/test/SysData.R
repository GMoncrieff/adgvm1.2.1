


SysData <- function(ta, te)
{
    filelist   <- dir(pattern="SysData_.*")
    filedata   <- file.info(filelist)
    filenewest <- sort(filedata$mtime, decreasing=T)[1]
    fileind    <- which(filedata$mtime==filenewest)
    filename   <- filelist[fileind]
    cat ( "\nRead ", filename, "\n\n" )

    t <- read.table(filename)
#     t <- read.table("SysData_0.dat")
    
    len <- dim(t)[1]
    x <- t[,1]+t[,2]/365
    
    cat( "Mean annual rain     ", mean(t[,33])*365, "\n" ) 
    cat( "Mean annual evapotr. ", mean(t[,34])*365, "\n" )
    
    cat( "\n" )

    graphics.off()
    
#     x11(width=20,height=11.5)
	pdf(file="SysData.pdf", width=10, height=12)
    par( mfcol=c(4,1), mar=c(2,4,1,1) )
    

# -------------------------------------------------------------------
# ------------- plot 1 ----------------------------------------------
# -------------------------------------------------------------------

    m  <- max( t[,3:6] )*10
    plot( -100, -100, ylab="Biomass (t/h)", xlab="Years", ylim=c(0,m), xlim=c(ta,te) )
    lines( x, t[,3]*10, col="green", lty=1, lwd=2 )
    lines( x, t[,4]*10, col="red",   lty=1, lwd=2 )
    lines( x, t[,5]*10, col="darkgreen", lty=1 )
    lines( x, t[,6]*10, col="darkgreen", lty=2 )
    legend( ta, m, c("grass leaf", "grass root", "standing dead grass", "lying dead grass"), col=c("green", "red", "darkgreen", "darkgreen"), lwd=c(2,2,1,1), lty=c(1,1,1,2), merge=T, cex=1.3, bty="n" )
    
# -------------------------------------------------------------------
# ------------- plot 2 ----------------------------------------------
# -------------------------------------------------------------------

    m  <- max( t[,14:20]*10 )
    plot( -100, -100, ylab="Biomass (t/ha)", xlab="Years", ylim=c(0,m), xlim=c(ta,te) )
    lines( x, t[,14]*10, col="green", lty=1, lwd=2 )
    lines( x, t[,15]*10, col="blue",  lty=1, lwd=2 )
    lines( x, t[,16]*10, col="red",   lty=1, lwd=2 )

    lines( x, t[,17]*10, col="darkgreen", lty=1 )
    lines( x, t[,18]*10, col="darkgreen", lty=2 )

    lines( x, t[,19]*10, col="darkblue", lty=1 )
    lines( x, t[,20]*10, col="darkblue", lty=2 )

    legend( ta, m, c("tree leaf", "tree stem", "tree root", "hanging dead leaf", "lying dead leaf", "standing dead stems", "lying dead stems"), lwd=c(2,2,2,1,1,1,1), col=c("green", "blue", "red", "darkgreen", "darkgreen", "darkblue", "darkblue"), lty=c(1,1,1,1,2,1,2), merge=T, cex=1.3, bty="n" )
    
# -------------------------------------------------------------------
# ------------- plot 3 ----------------------------------------------
# -------------------------------------------------------------------

    plot( -100, -100, ylab="Cover", xlab="Years", ylim=c(0,1), xlim=c(ta,te) )
    lines( x, t[,11], col="green", lty=1, lwd=2 )
    lines( x, t[,12], col="blue",  lty=1, lwd=2 )
    lines( x, t[,13], col="red",   lty=1, lwd=2 )
    lines( x, t[,28], col="grey",  lty=1, lwd=2 )


    legend( ta, 1, c("Savanna tree cover", "Forest tree cover", "Burning bush cover", "C3:C4 ratio"), lwd=c(2,2,2,2), col=c("green", "blue", "red", "grey"), lty=c(1,1,1,1), merge=T, cex=1.3, bty="n" )



# -------------------------------------------------------------------
# ------------- plot 2 ----------------------------------------------
# -------------------------------------------------------------------

# NPP plots

    npp_gr <- t[, 8]-t[, 9]-t[,10]
    npp_tr <- t[,22]-t[,23]-t[,24]
    ann_nee <- rep(0,max(t[,1]))
    ann_npp <- rep(0,max(t[,1]))

    for ( i in 1:length(ann_nee) )
    {
        tmp <- subset( t, t[,1]==i )
#         ann_nee[i] <- sum(tmp[,8]-tmp[,9]-tmp[,10]+tmp[,22]-tmp[,23]-tmp[,24]-tmp[,35]-tmp[,36])
        ann_nee[i] <- sum(tmp[,8]-tmp[,9]-tmp[,10]+tmp[,22]-tmp[,23]-tmp[,24]-tmp[,35])
        ann_npp[i] <- sum(tmp[,8]-tmp[,9]-tmp[,10]+tmp[,22]-tmp[,23]-tmp[,24])
    }
    ann_nee <- ann_nee*t[2,2]  # data are not written every day
    ann_npp <- ann_npp*t[2,2]  #  --> multiplication with timestep

    cat( "Mean annual npp      ", mean(ann_npp), " kg C/m^2/year\n" )
    cat( "Mean annual nee      ", mean(ann_nee), " kg C/m^2/year\n" )

    ann_nee <- ann_nee/100   # 100 is scaling factor for plots
    ann_npp <- ann_npp/100   # 100 is scaling factor for plots

    m  <- max( max(npp_gr+npp_tr-t[,35]), max(ann_nee), max(ann_npp) )
    n  <- min( min(npp_gr+npp_tr-t[,35]), min(ann_nee), max(ann_npp) )

    plot( -100, -100, ylab="Fluxes (kg C/m^2/day or kg C/m^2/year/100)", xlab="Years", ylim=c(n,m), xlim=c(ta,te) )
    lines( x, npp_gr, col="darkred", lty=1, lwd=2 )
    lines( x, npp_tr, col="darkgreen", lty=1, lwd=2 )
    lines( x, -t[,35], col="darkblue", lty=1, lwd=2 )  # soil
# 	lines( x, npp_gr+npp_tr-t[,34], col="black", lwd=1 )
# 	points( x, -t[,35]/100, pch=20, lwd=5, col="magenta" )  # fire 
	points( (1:length(ann_nee))+1., ann_nee, pch=20, lwd=5, col="cyan" )
	lines(  (1:length(ann_nee))+1., ann_nee, col="cyan" )
	points( (1:length(ann_npp))+1., ann_npp, pch=20, lwd=5, col="magenta" )
	lines(  (1:length(ann_npp))+1., ann_npp, col="magenta" )
    abline( h=0 )

    legend( ta, m, c("NPP grass", "NPP tree", "Soil emmisions", "NPP annual/100", "NEE annual/100"), col=c("darkred", "darkgreen", "darkblue", "magenta", "cyan"), lwd=c(2,2,2,5,5), lty=c(1,1,1,3,3), merge=T, cex=1.3, bty="n" )






    graphics.off()


}











