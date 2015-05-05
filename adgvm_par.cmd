#! /bin/bash

#PBS -N Glenn_aDGVM1.2.1

#PBS -l nodes=5:ppn=24
#Memory required per processor (specify how many megabytes)
#PBS -l pmem=1000mb
# You must specify Wall Clock time (hh:mm:ss) [Maximum allowed 30 days = 720:00:00]
#PBS -l walltime=360:00:00
#PBS -o myjob.output
# Send job stdout to file "myjob.output"
#PBS -e myjob.err


homefolder="/home/gmoncrieff/adgvm/adgvm1/DATE_BRANCH/result"
current="/home/gmoncrieff/adgvm/adgvm1/DATE_BRANCH/"

mv -f Globals_cluster.h Globals.h

g++ -O3 -o adgvm Model.cpp

mpic++ -g -O0  -c call_adgvm_par.cpp -o call_adgvm_par.o
mpic++ -g  call_adgvm_par.o -o call_adgvm_par

mkdir $homefolder

cd $current

cp adgvm call_adgvm_par sites.dat $homefolder

cd $homefolder

mkdir OutputData

chmod +x adgvm call_adgvm_par

mpirun ./call_adgvm_par

echo DONE

