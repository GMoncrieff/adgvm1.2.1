library(doMC)
registerDoMC(detectCores())

Outfolder <- "OutputFolderName"
Modelhome <- "/Users/glennmoncrieff/adgvm1/R/"
Datahome <- "/Users/glennmoncrieff/adgvm/data/OutputData"
Currentmodel <-"LOCATIONXXX/DATEXXX"
setwd(paste(Modelhome,Currentmodel,sep="")
      
#parameters
modelname <- "adgvm"  #name of compiled model
ntrees <- 100         #number of tress at model start
fire <- 1             # 0 = off, 1 = on
years <- 100          # number of years to run
climscen <- 0         # 0 = ambient, 3 = IPCC A1B

#read in list of sites
all_ba<-read.table("sites.csv",sep=",") #read data with coords and ba vals

#create a vector to store model calls. one call for each site
c_line<-numeric(nrow(all_ba))

#loop through sites in parallel
foreach (i=1:nrow(all_ba), .combine=c) %dopar% { 

      #random number seed
      randseed_ <- Sys.time()
      randseed <- as.integer(randseed_)
      
      #create call
      c_line[i] <- paste( "./", modelname, ntrees, years, all_ba[i,1], all_ba[i,2], fire ,randseed, climscen, Outfolder, sep = " ")
      
      # call adgvm
      system(c_line[i])
      
      numnow<- i
      numnow
      
  }


setwd(paste(Datahome,Outfolder,sep="")

# read in results
res    <- read.table("YearlyData."years,climscen,fire,".dat")

save(res,file=paste(Datahome,Outfolder,"adgvm_data.Rdata",sep=""))

