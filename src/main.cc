//WL ground state searching

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <chrono>
#include "mpi.h"
#include "parameter.hpp"
#include "rand.hpp"
#include "alloy.hpp"
#include "pt.hpp"
#include "wanglandau.hpp"
#include "allreduce.h"
#include "client.hpp"

int nprocs;
int myrank;
//#define NDEBUG

int main(int argc, char *argv[])
{
 
	int i,j,k,m;
	int STARTED;  //Shows whether or not the restart has been applied
	FILE *ofp_run, *ofp_comm;
	char s[512];
	double tmp_flat=0.0, volume;
	double ptEmin, ptEmax;
	int nT, nsweeps;
	time_t start, end;
	std::chrono::high_resolution_clock::time_point tik,tok;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	//ini random number;
	srand(atoi(argv[4])*myrank);
	shelltimeseed(rand()+19*myrank+19273);

        std::string model_dir, vae_dir;
        //model_dir = "./models/exported/mixed-opt/dense/MoNbTaW";
        model_dir = "./models/exported/mixed-opt/conv/MoNbTaW";
        //model_dir = "./models/exported/fp32-opt/small/MoNbTaW";
        //model_dir = "./models/exported/fp32-opt/large/MoNbTaW/model_MoNbTaW";
        //model_dir = "./models/exported/fp16-opt/large/MoNbTaW/model_MoNbTaW";

#ifdef DL_MODEL

#ifdef BACKEND_TF
        //start tensorflow session;
        SessionOptions options;
        options.config.mutable_gpu_options()->set_visible_device_list(std::to_string(myrank%GPUperNode));
        model.LoadModel(model_dir, options);
#else
	LoadClientModel(model_dir);	

#endif

#endif
        vae_dir = "./models/vae/models";
	LoadClientModel(vae_dir);	
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

	//warm up with parallel tempering 
	parallel_tempering(nprocs,MDROP,MSAMPS,MSEP,0,metropolis);
	MPI_Allreduce(&currEtot, &ptEmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&currEtot, &ptEmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	if(ptEmin/N_3 < WLD1min)
		WLD1min = ptEmin/N_3;
	if(ptEmax/N_3 > WLD1max)
		WLD1max = ptEmax/N_3;
	printf("ptEmin = %g  ptEmax = %g\n", ptEmin, ptEmax);
	
        if(myrank == 0){
                ofp_run=fopen("run.dat","a");
                ofp_comm=fopen("comm.dat","w");
		start = time(NULL);
                fprintf(ofp_run,"#seeds: %d,%d,%d\n",314159265,362436069,atoi(argv[4]));
        }


	//Initialize the Wang-Landau sampling parameters
	initWL();
        numf=1;
	STARTED = 0;
	size_t msgSize = D1BINS*sizeof(double);
        nsweeps = N*N*N/nprocs > 0 ? N*N*N/nprocs : 1 ;
	for(lnwlf=1.0;lnwlf>ModFactorFinal;lnwlf=lnwlf/IterationFactor)
	{
		IterSweeps=0;
		resetWL();
		tmp_flat=0.0;
		MPI_Barrier(MPI_COMM_WORLD);
#ifdef  GLOBAL_UPDATE
		global_update(N_3, Flatness*lnwlf);
#endif		
		while(tmp_flat <= Flatness)
		{		
			if(TotalSweeps%int(pow(10, numf))==0)
                                vae_update(WLdos);

			sweepWL(N_3, WLdos);
			IterSweeps+=1;
			TotalSweeps+=1; 
			if( (TotalSweeps % 1) == 0)
			{
				MPI_Barrier(MPI_COMM_WORLD);
				if(IterSweeps == 10 && myrank == 0)
					tik = std::chrono::high_resolution_clock::now();
#ifdef COMM_RING
				RingAllreduce(wllngi,msgSize,&wllngd,myrank,nprocs);
				RingAllreduce(wlHi,msgSize,&wlHd,myrank,nprocs);
				MPI_Barrier(MPI_COMM_WORLD);
#else
				//MPI_Allreduce(wllngi, wllngd, D1BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				//MPI_Allreduce(wlHi, wlHd, D1BINS, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
				MPI_Allreduce(wlHi, wlHd, D1BINS, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);
#endif
				if(IterSweeps == 10 && myrank == 0)
					tok = std::chrono::high_resolution_clock::now();
			        for(i=0;i<D1BINS;i++){
					wlH[i] += wlHd[i];
					wllng[i] += wlHd[i]*lnwlf;    //wllngd[i];
				}	
				memset(wllngi, 0, D1BINS*sizeof(double));
				//memset(wlHi, 0, D1BINS*sizeof(int));
				memset(wlHi, 0, D1BINS*sizeof(unsigned short));
				tmp_flat=flatWL();
				if( (IterSweeps %100)==0 &&  myrank == 0)
					write_DOS_H();
			};
		}
		if(myrank == 0)
			write_DOS_H();
		if(myrank == 0){
                        fprintf(ofp_run,"%g\t%d\t%d\t%g\t%g\n",lnwlf,TotalSweeps*D1BINS,IterSweeps*D1BINS,tmp_flat,numbelow_flat);
			std::chrono::duration<double, std::milli> ms_double = tok - tik;
                        fprintf(ofp_comm, "comm time: %f (ms)\nbandwidth: %f (GB/s)\n", 1.0*ms_double.count(), 2.0*msgSize/1e+6/ms_double.count());
                        fflush(ofp_run);
                        fflush(ofp_comm);

                }
		numf=numf+1;
	}
	if(myrank == 0){
		end = time(NULL);
		fprintf(ofp_run, "wl time: %ld (s) \n", end-start); 
		fclose(ofp_run);
	}

	// production runs
        resetWL();              
        int count=0;
	MPI_Barrier(MPI_COMM_WORLD);
        while(count < ProductionBinSamps){ // (lnwlf > 1e-6){
		sweepWL(N_3, WLproduction);
		//std::cout << "myrank: " << myrank << " E: " << currEtot << std::endl;
		if( (IterSweeps % 1 ) == 0)
		{
                        MPI_Allreduce(wlHi, wlHd, D1BINS, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
			//MPI_Allreduce(wlHi, wlHd, D1BINS, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                        for(i=0;i<D1BINS;i++){
                        	wlH[i] += wlHd[i];
                                wllng[i] += wlHd[i]*lnwlf; 
                        }
			//memset(wlHi, 0, D1BINS*sizeof(int));
			memset(wlHi, 0, D1BINS*sizeof(unsigned short));
                }
		count++;
                IterSweeps+=1;
                TotalSweeps+=1;

                lnwlf = 1.0/TotalSweeps;
        }
	memset(wllngi, 0, D1BINS*sizeof(double));
	memset(wllngd, 0, D1BINS*sizeof(double));
	MPI_Reduce(op, wllngi, D1BINS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(op2, wllngd, D1BINS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if(myrank == 0){
                write_DOS_H();
                thermoqs();
        }
	
	freeWL();
	freePT();
	MPI_Finalize();
}



