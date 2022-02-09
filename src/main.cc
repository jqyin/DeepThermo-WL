//WL ground state searching

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mpi.h"
#include "parameter.hpp"
#include "rand.hpp"
#include "alloy.hpp"
#include "pt.hpp"
#include "wanglandau.hpp"
#include "allreduce.h"

int nprocs;
int myrank;
//#define NDEBUG

int main(int argc, char *argv[])
{
 
	int i,j,k,m;
	int STARTED;  //Shows whether or not the restart has been applied
	FILE *ofp_run;
	char s[512];
	double tmp_flat=0.0, volume;
	int nT;
	time_t t1,t2;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	//ini random number;
	srand(atoi(argv[4])*myrank);
	shelltimeseed(rand()+19*myrank+19273);

        //start tensorflow session;
        std::string model_dir;
        SessionOptions options;
        options.config.mutable_gpu_options()->set_visible_device_list(std::to_string(myrank%6));
        model_dir = "./models/exported/fp32-opt/medium/MoNbTaVW/model_MoNbTaVW";
        model.LoadModel(model_dir, options);
        

	//Error message if the number of arguments is incorrect
	if (argc != 3 && myrank == 0) ErrorMsg(0, "");

	//Reads Input Parameters
	ReadInput(argv[1]);
  	//Initializes the System
	ini_T(MTi,MTf,nprocs);
	ini_sys();

        if(myrank == 0)
            ini_alloy(0);
        MPI_Bcast(Atom, N_3, MPI_SHORT, 0, MPI_COMM_WORLD);
        MPI_Bcast(NT, NE, MPI_INT, 0, MPI_COMM_WORLD);
        ini_apos();
        currEtot = Etot();
	printf("Etot = %g\n", currEtot);

	//warm up
	parallel_tempering(nprocs,MDROP,MSAMPS,MSEP,0,0);

        if(myrank == 0){
                ofp_run=fopen("run.dat","a");
                fprintf(ofp_run,"#seeds: %d,%d,%d\n",314159265,362436069,atoi(argv[4]));
        }
	//Initialize the Wang-Landau sampling parameters
	initWL();
        numf=1;
	STARTED = 0;
	size_t msgSize = D1BINS*sizeof(double);
	for(lnwlf=1.0;lnwlf>ModFactorFinal;lnwlf=lnwlf/IterationFactor)
	{
		IterSweeps=0;
		resetWL();
		tmp_flat=0.0;
		MPI_Barrier(MPI_COMM_WORLD);
#ifdef  GLOBAL_UPDATE
		if(myrank == 0)
			global_update();
		MPI_Bcast(wlH, D1BINS, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
		MPI_Bcast(wllng, D1BINS, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
#endif		
		while(tmp_flat <= Flatness)
		{		
			sweepWL(N*N*N/nprocs, 1);
			IterSweeps+=1;
			TotalSweeps+=1; 
			if( (TotalSweeps % 10) == 0)
			{
				RingAllreduce(wllngi,msgSize,&wllngd);
				RingAllreduce(wlHi,msgSize,&wlHd);
				//MPI_Allreduce(wllngi, wllngd, D1BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				//MPI_Allreduce(wlHi, wlHd, D1BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			        for(i=0;i<D1BINS;i++){
					wlH[i] += wlHd[i];
					wllng[i] += wllngd[i];
				}	
				memset(wllngi, 0, D1BINS*sizeof(double));
				memset(wlHi, 0, D1BINS*sizeof(double));
				tmp_flat=flatWL();
				//if( (IterSweeps %10)==0 &&  myrank == 0)
				if(myrank == 0)
					 write_DOS_H();
			};
		}
		if(myrank == 0)
			write_DOS_H();
		if(myrank == 0){
                        fprintf(ofp_run,"%g\t%d\t%d\t%g\t%g\n",lnwlf,TotalSweeps*D1BINS,IterSweeps*D1BINS,tmp_flat,numbelow_flat);
                        fflush(ofp_run);
                }
		numf=numf+1;
	}
	if(myrank == 0)
		fclose(ofp_run);

        /*resetWL();              
        int count=0;
	MPI_Barrier(MPI_COMM_WORLD);

        while(count < ProductionBinSamps){ // (lnwlf > 1e-6){
		sweepWL( N*N*N/nprocs, 1);
		std::cout << "myrank: " << myrank << " E: " << currEtot << std::endl;
		if( (IterSweeps % 1 ) == 0)
		{
                	MPI_Allreduce(wllngi, wllngd, D1BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                        MPI_Allreduce(wlHi, wlHd, D1BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                        for(i=0;i<D1BINS;i++){
                        	wlH[i] += wlHd[i];
                                wllng[i] += wllngd[i];
                        }
                        memset(wllngi, 0, D1BINS*sizeof(double));
                        memset(wlHi, 0, D1BINS*sizeof(double));
                }
		count++;
                IterSweeps+=1;
                TotalSweeps+=1;

                lnwlf = 1.0/TotalSweeps;
        }*/
        if(myrank == 0){
                write_DOS_H();
                thermoqs();
        }
	
	freeWL();
	freePT();
	MPI_Finalize();
}



