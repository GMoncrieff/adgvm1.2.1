

# IPCC SRES future simulations

# system("g++ -O3 -o Model_za_fut Model.cpp -DCLIM_SCEN")
# g++ -O3 -o Model_za_fut Model.cpp -DCLIM_SCEN -DS_NOTREES

coords <- read.table("south_africa_coo.dat")

# unlink("job*")

# run_vec <- 1:4
run_vec <- 0

fir_vec <- c( 1, 0 )
# fir_vec <- c( 1 )
cnt    <- 0

# use odd number to have with/without fire simulations in one file
# n_proc <- 15   ## number of processors to be used
n_proc <- 11   ## number of processors to be used


cseq <- sample(dim(coords)[1])

# for ( x in 1:dim(coords)[1] )
for ( x in cseq )
{
  for ( fire in fir_vec )
  {
    for ( run in run_vec )
    {
#       cline <- paste( "./Model_za_fut 100 339", coords[x,1], coords[x,2], fire, run, " 2 za_fut_A2" )
      cline <- paste( "./Model_za_fut 1000 339", coords[x,1], coords[x,2], fire, run, " 3 za_fut_A1B" )
#        cline <- paste( "./Model_za_fut 100 339", coords[x,1], coords[x,2], fire, run, " 4 za_fut_B1" )

#       cline <- paste( "./Model_za_fut 3 1", coords[x,1], coords[x,2], fire, run, " 0 ttt" )


      write.table(cline, file=paste("job_", cnt, sep=""), append=T, col.names=F, row.names=F, quote=F )
      cnt <- cnt+1
              
      if ( cnt==n_proc ) cnt <- 0
    }
  }
}

for ( i in 0:(n_proc-1) )
  write.table(paste("./job_", i, " &", sep=""), file="jobs.run", append=T, col.names=F, row.names=F, quote=F )

system("chmod +x job*")




















