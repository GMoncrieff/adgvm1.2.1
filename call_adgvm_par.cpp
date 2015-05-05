#include <mpi.h>
#include <iostream>
#include <string>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <sstream>
#include <vector>
#include <stdio.h>

using namespace std;

//code for replacing # in site file
inline string convert(string str_in)
{
     for(int i = 0; i < str_in.length(); i++)
     {
          switch(str_in[i])
          {
               case '#':
               str_in[i] = ' ';
               break;
          }
     }

     return str_in;
}

int main(int argc, char *argv[])
{
	cout << "start" << endl;
	int rank, size;
	  //processor name
	int namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int mystart, myend; // loop boundaries
	 int x;


// adgvm pars
	string EXE1 = "./adgvm";
//	string EXE2 = "./adgvm_for";
	string SIM_SITE;
	string NTREE = "100";
	string NYEAR = "100";
	string FIRE = "1";
	string CLIM = "0";
	string DIR1 = "OutputData";
//	string DIR2 = "OutputData2";
/*	string CLD = "20";
	string CAN_EX1 = "50";
	string CAN_EX2 = "40";
	string IG1 = "10";
	string IG2 = "10";
	string USTEM = "15";
	string RESP_T = "35";  //35
	string RESP_G = "35";  //35
	string SUCK1 = "0";
	string SUCK2 = "0";
	string BETA = "25";
	string GERM1 = "25";  //25
	string GERM2 = "25";  //25
	string DFROST1 = "10";
	string DFROST2 = "10";
	string DCARB1 = "10";  //10
	string DCARB2 = "10";  //10
	string DCOMP1 = "10";
	string DCOMP2 = "50";
	string C_SS = "50";
	string C_S4 = "50";
	string C_S3 = "15";
	string C_FS = "50";
	string TK_C1 = in1; //43
	string TK_C2 = "53";
	string TK_H1 = in2; //5003
	string TK_H2 = "4003";
	string TK_I1 = in3;	//4408
	string TK_I2 = "5408";
	string D1_1 = in4;
	string D1_2 = "382";
	string D2_1 = in5;
	string D2_2 = "5375"; */


// read in list of sites
	ifstream sitefile;
	sitefile.open("sites.dat");
	vector<string> sites;
	string asite;

	while (sitefile >> asite) {

	  sites.push_back(asite);

	}
	sitefile.close();

	float RESMAT [sites.size()][2];

	int nsites = sites.size();

	//initialize parralelization
	  MPI_Init(&argc, &argv);
	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	  MPI_Comm_size(MPI_COMM_WORLD, &size);
	  MPI_Get_processor_name(processor_name,&namelen);


	  mystart = (nsites / size) * rank;
	    if (nsites % size > rank){
	      mystart += rank;
	      myend = mystart + (nsites / size) + 1;
	    }else{
	      mystart += nsites % size;
	      myend = mystart + (nsites / size);
	    }

	  printf("CPU%d %d ~ %d\n",rank,mystart,myend);


	// loop through sites
	  for(x=mystart; x<myend; x++){
	   // cout << sites.at( x ) << " ";
	   SIM_SITE = sites.at( x );
	   SIM_SITE = convert(SIM_SITE);

// generate random seed, first digit is site number so that list can be sorted later.
	   string NRAN;
	   stringstream out;
	   out << x;
	   NRAN = out.str();

	   srand(time(0));
	   int RAN1 = (rand() % 10);
	   srand(RAN1);

	   string RAN11;
	   stringstream out1;
	   out1 << RAN1;
	   RAN11 = out1.str();


	   string IRAN = NRAN+RAN11;

	    //cout << IRAN << endl;

	   string CALL1 = EXE1+" "+NTREE+" "+NYEAR+" "+SIM_SITE+" "+FIRE+" "+IRAN+" "+CLIM+" "+DIR1;
//	   string CALL2 = EXE2+" "+NTREE+" "+NYEAR+" "+SIM_SITE+" "+FIRE+" "+IRAN+" "+CLIM+" "+DIR2;
	/*	+" "+CLD+" "+CAN_EX1+" "+CAN_EX2+" "+IG1+" "+IG2+" "+USTEM+" "+RESP_T+" "+RESP_G+" "+SUCK1+" "+SUCK2+" "+BETA+" "+GERM1+" "+GERM2+" "
	   +DFROST1+" "+DFROST2+" "+DCARB1+" "+DCARB2+" "+DCOMP1+" "+DCOMP2+" "+C_SS+" "+C_S4+" "+C_S3+" "+C_FS+" "+TK_C1+" "+TK_C2+" "+TK_H1+" "+TK_H2+" "+TK_I1+" "+TK_I2+" "+D1_1+" "+D1_2+" "+D2_1+" "+D2_2;
	*/
	   system(CALL1.c_str());
//	   system(CALL2.c_str());
	}

	MPI_Barrier(MPI_COMM_WORLD);

	//control passes to master
	if(rank == 0){
    cout << "DONE" << endl;
	}
	MPI_Finalize();
        return 0;
}
