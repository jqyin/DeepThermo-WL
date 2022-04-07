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

#ifndef TF_BACKEND
#include "client.hpp"
#endif

int main(int argc, char *argv[])
{
 
	int i,j,k,m;
	FILE *ofp_run, *ofp_comm;
	char s[512];
	double tmp_flat=0.0, volume;
	double ptEmin, ptEmax;
	int nT, nsweeps;
	time_t start, end;
	std::chrono::high_resolution_clock::time_point tik,tok;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpiState.nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpiState.myrank);
	//ini random number;
	srand(atoi(argv[4])*mpiState.myrank);
	shelltimeseed(rand()+19*mpiState.myrank+19273);

        std::string model_dir;
        model_dir = "./models/";

#ifdef TF_BACKEND
        //start tensorflow session;
        SessionOptions options;
        options.config.mutable_gpu_options()->set_visible_device_list(std::to_string(mpiState.myrank%GPUperNode));
        model[0].LoadModel(model_dir+"/encoder", options);
        model[1].LoadModel(model_dir+"/decoder", options);
#else
	LoadClientModel(model_dir);	
#endif
	//Error message if the number of arguments is incorrect
	if (argc != 3 && mpiState.myrank == 0) ErrorMsg(0, "");

	//Reads Input Parameters
	ReadInput(argv[1]);
  	//Initializes the System
	ini_T(ptState.MTi,ptState.MTf,mpiState.nprocs);
	ini_sys();

        if(mpiState.myrank == 0)
            ini_alloy(0);
        MPI_Bcast(alloyState.Atom, alloyState.N_3, MPI_SHORT, 0, MPI_COMM_WORLD);
        MPI_Bcast(alloyState.NT, NE, MPI_INT, 0, MPI_COMM_WORLD);
        ini_apos();
        alloyState.currEtot = Etot();
	printf("Etot = %g\n", alloyState.currEtot);

	//warm up with parallel tempering 
	parallel_tempering(mpiState.nprocs,ptState.MDROP,ptState.MSAMPS,ptState.MSEP,0,metropolis);
	MPI_Allreduce(&alloyState.currEtot, &ptEmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&alloyState.currEtot, &ptEmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	if(ptEmin/alloyState.N_3 < wlState.WLD1min)
		wlState.WLD1min = ptEmin/alloyState.N_3;
	if(ptEmax/alloyState.N_3 > wlState.WLD1max)
		wlState.WLD1max = ptEmax/alloyState.N_3;
	printf("ptEmin = %g  ptEmax = %g\n", ptEmin, ptEmax);
	
        if(mpiState.myrank == 0){
                ofp_run=fopen("run.dat","a");
                ofp_comm=fopen("comm.dat","w");
		start = time(NULL);
                fprintf(ofp_run,"#seeds: %d,%d,%d\n",314159265,362436069,atoi(argv[4]));
        }


	//Initialize the Wang-Landau sampling parameters
        wlState.numf=1;
	size_t msgSize = wlState.D1BINS*sizeof(double);
        nsweeps = alloyState.N*alloyState.N*alloyState.N/mpiState.nprocs > 0 ? alloyState.N*alloyState.N*alloyState.N/mpiState.nprocs : 1 ;
	initWL();
	if(mpiState.myrank==0)
		printf("Ecut: %f\n", (alloyState.E_0/wlState.invdWLD1+wlState.WLD1min)*alloyState.N_3);
	for(wlState.lnwlf=wlState.ModFactorInit;wlState.lnwlf>wlState.ModFactorFinal;wlState.lnwlf=wlState.lnwlf/wlState.IterationFactor)
	{
		wlState.IterSweeps=0;
		resetWL();
		tmp_flat=0.0;
		MPI_Barrier(MPI_COMM_WORLD);
#ifdef  GLOBAL_UPDATE
		global_update(alloyState.N_3, wlState.Flatness*wlState.lnwlf);
#endif		
		while(tmp_flat <= wlState.Flatness)
		{		
			if(wlState.TotalSweeps%int(pow(2, wlState.numf))==0 && alloyState.currEtot < (alloyState.E_0/wlState.invdWLD1+wlState.WLD1min)*alloyState.N_3)
                                vae_update(WLdos);

			sweepWL(alloyState.N_3, WLdos);
			wlState.IterSweeps+=1;
			wlState.TotalSweeps+=1; 
			if( (wlState.TotalSweeps % 1) == 0)
			{
				MPI_Barrier(MPI_COMM_WORLD);
				if(wlState.IterSweeps == 10 && mpiState.myrank == 0)
					tik = std::chrono::high_resolution_clock::now();
#ifdef COMM_RING
				RingAllreduce(wlState.wllngi,msgSize,&wlState.wllngd,mpiState.myrank,mpiState.nprocs);
				RingAllreduce(wlState.wlHi,msgSize,&wlState.wlHd,mpiState.myrank,mpiState.nprocs);
				MPI_Barrier(MPI_COMM_WORLD);
#else
				MPI_Allreduce(wlState.wlHi, wlState.wlHd, wlState.D1BINS, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);
#endif
				if(wlState.IterSweeps == 10 && mpiState.myrank == 0)
					tok = std::chrono::high_resolution_clock::now();
			        for(i=0;i<wlState.D1BINS;i++){
					wlState.wlH[i] += wlState.wlHd[i];
					wlState.wllng[i] += wlState.wlHd[i]*wlState.lnwlf;    
				}	
				memset(wlState.wllngi, 0, wlState.D1BINS*sizeof(double));
				memset(wlState.wlHi, 0, wlState.D1BINS*sizeof(unsigned short));
				tmp_flat=flatWL(WLdos);
				if( (wlState.IterSweeps %100)==0 &&  mpiState.myrank == 0)
					write_DOS_H();
			};
		}
		if(mpiState.myrank == 0)
			write_DOS_H();
		if(mpiState.myrank == 0){
                        fprintf(ofp_run,"%g\t%d\t%d\t%g\t%g\n",wlState.lnwlf,wlState.TotalSweeps*wlState.D1BINS,wlState.IterSweeps*wlState.D1BINS,tmp_flat,wlState.numbelow_flat);
			std::chrono::duration<double, std::milli> ms_double = tok - tik;
                        fprintf(ofp_comm, "comm time: %f (ms)\nbandwidth: %f (GB/s)\n", 1.0*ms_double.count(), 2.0*msgSize/1e+6/ms_double.count());
                        fflush(ofp_run);
                        fflush(ofp_comm);

                }
		wlState.numf=wlState.numf+1;
		alloyState.E_0 = alloyState.E_0 - int(0.01*wlState.D1BINS);  
	}
	if(mpiState.myrank == 0){
		end = time(NULL);
		fprintf(ofp_run, "wl time: %ld (s) \n", end-start); 
		fclose(ofp_run);
	}

	// production runs
        resetWL();              
        int count=0;
	MPI_Barrier(MPI_COMM_WORLD);
        while(count < wlState.ProductionBinSamps){ 
		sweepWL(alloyState.N_3, WLproduction);
		if( (wlState.IterSweeps % 1 ) == 0)
		{
                        MPI_Allreduce(wlState.wlHi, wlState.wlHd, wlState.D1BINS, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
                        for(i=0;i<wlState.D1BINS;i++){
                        	wlState.wlH[i] += wlState.wlHd[i];
                                wlState.wllng[i] += wlState.wlHd[i]*wlState.lnwlf; 
                        }
			memset(wlState.wlHi, 0, wlState.D1BINS*sizeof(unsigned short));
                }
		count++;
                wlState.IterSweeps+=1;
                wlState.TotalSweeps+=1;

        }
	memset(wlState.wllngi, 0, wlState.D1BINS*sizeof(double));
	memset(wlState.wllngd, 0, wlState.D1BINS*sizeof(double));
	MPI_Reduce(alloyState.op, wlState.wllngi, wlState.D1BINS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(alloyState.op2, wlState.wllngd, wlState.D1BINS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if(mpiState.myrank == 0){
                write_DOS_H();
                thermoqs();
        }
	
	freeWL();
	freePT();
	MPI_Finalize();
}



