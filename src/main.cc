//Classical Hisenberg Model

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mpi.h"
#include "parameter.hpp"
#include "rand.hpp"
#include "alloy.hpp"
#include "wanglandau.hpp"

#define Nprof 32 
#define SyncT 1

int nprocs;
int myrank;

//#define NDEBUG
double fit(double* x,double* y,int ndata);
void update_prof();
void recalibrate();

//#define NDEBUG

int main(int argc, char *argv[])
{
 
	int i,j,k,m;
	int STARTED;  //Shows whether or not the restart has been applied
	FILE *ofp_run, *fp_err;
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

	//Error message if the number of arguments is incorrect
	if (argc !=5 && argc!=6) ErrorMsg(0, "");

	//Reads Input Parameters
	ReadInput(argv[1]);

	//Get Random Number Seed
	//shelltimeseed(atoi(argv[4]));
	//gettimeseed();
  
  	//Initializes the System
	if(argc==5 )
		initialize( NULL);
	else
		initialize(argv[5]);

	//Initialize the Wang-Landau sampling parameters

	initWL();
	nummoves=0;  //used in wlhybrid() to keep track of total number of moves

	ProductionRun = atoi(argv[3]);
	//starts the production run if it is turned on in the input file
	if( ProductionRun == 1)
	{
		production_run();
		exit(1);
	}
	else if( (ProductionRun != 0) && (ProductionRun != 1) )
	{
		fprintf(stderr,"ProductionRun is not equal to 1 or 0! Is set to: %d \n\n",ProductionRun);
		exit(1);
	};
	
	if(myrank == 0){
		ofp_run=fopen("run.dat","a");
		fprintf(ofp_run,"#seeds: %d,%d,%d\n",314159265,362436069,atoi(argv[4]));
		fp_err = fopen("err.dat", "w");
	}
	//Initialization
	numf=1;  //Keeps track of the number of interations  
	STARTED = 0;
	//Main Wang-Landau sampling loop
	for(lnwlf=1.0;lnwlf>ModFactorFinal;lnwlf=lnwlf/IterationFactor)
    {
		IterSweeps=0;	//Keeps track of the number of sweeps per iteration
		resetWL();		//Resets the WL histgram      
		tmp_flat=0.0;	//Resets the flatness for the loop below
		MPI_Barrier(MPI_COMM_WORLD);
#ifdef  GLOBAL_UPDATE
		if(myrank == 0)
			global_update();
		MPI_Bcast(wlH, D1BINS, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
		MPI_Bcast(wllng, D1BINS, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
#endif		
	/*	volume = 0.0;
		for(i=0;i<N;i++)	
			for(j=0;j<N;j++)	
				for(k=0;k<N;k++)
					volume += vol[0][i][j][k] + volb[0][i][j][k];	
		if(myrank == 0){
	  		 fprintf(ofp_run, "trackE = %18.10e, trueE= %18.10e, Volume= %18.10e\n", currEtot, Etot(), volume); 
			 fflush(ofp_run);
		}*/

  
		//time(&t1);
		//This while loop is executed until the histogram is sampled
		while(  tmp_flat <= Flatness  )
		{
		
			//Performs a number () of WL sweeps
			sweepWL( N*N );
			IterSweeps+=1;		
			TotalSweeps+=1;
	 		
			if( (IterSweeps % SyncT ) == 0)
			{// after a long run, still not flat;
			//	time(&t2);
			//	printf("Time to do 10 MCS= %d second",(int)(t2-t1));		
				MPI_Allreduce(wllngi, wllngt, D1BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				MPI_Allreduce(wlHi, wlHt, D1BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			        for(i=0;i<D1BINS;i++){
					wlH[i] += wlHt[i];
					wllng[i] += wllngt[i];
				}	
				memset(wllngi, 0, D1BINS*sizeof(double));
				memset(wlHi, 0, D1BINS*sizeof(double));
				tmp_flat=flatWL();
				if( (IterSweeps %100)==0 &&  myrank == 0)
					write_DOS_H();
				//	feedback_dt();

			};

		};
		if(myrank == 0)
			write_DOS_H();

		//thermoqs();
		if(myrank == 0){
			fprintf(ofp_run,"%g\t%d\t%d\t%g\t%g\n",lnwlf,TotalSweeps*D1BINS,IterSweeps*D1BINS,tmp_flat,numbelow_flat);
			fflush(ofp_run);
		}
		numf=numf+1;  //Used in write_DOS_H() to lable each DOS
	};
	if(myrank == 0)
		fclose(ofp_run); 	

#ifdef PROF
	MPI_Allreduce(wlHui, wlHu, D1BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(wlHdi, wlHd, D1BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	update_prof();
#endif 	

	resetWL();		//Resets the WL histgram      
	int count=0;
	while(count < ProductionBinSamps){ // (lnwlf > 1e-6){
		
			sweepWL( N*N/nprocs );
			if( (IterSweeps % SyncT ) == 0)
			{// after a long run, still not flat;
			//	time(&t2);
			//	printf("Time to do 10 MCS= %d second",(int)(t2-t1));		
				MPI_Allreduce(wllngi, wllngt, D1BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				MPI_Allreduce(wlHi, wlHt, D1BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			        for(i=0;i<D1BINS;i++){
					wlH[i] += wlHt[i];
					wllng[i] += wllngt[i];
				}	
				memset(wllngi, 0, D1BINS*sizeof(double));
				memset(wlHi, 0, D1BINS*sizeof(double));
			}
			if(TotalSweeps % 10000 == 0 && myrank == 0){
				count_nn(TotalSweeps, fp_err);
				write_DOS_H();
			}
			count++;
			IterSweeps+=1;		
			TotalSweeps+=1;

			lnwlf = 1.0/TotalSweeps;
#ifdef PROF
			if(TotalSweeps % 10000 == 0){
				MPI_Allreduce(wlHui, wlHu, D1BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				MPI_Allreduce(wlHdi, wlHd, D1BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				update_prof();			
			}

#endif 	
	}

	//Prints out the final normalized density of states, thermodynamics, and order parameters
	//write_DOS_H();
	//write_restart();
	//write_normDOS();
	if(myrank == 0){
		write_DOS_H();
		recalibrate();
		thermoqs();
	}
	freeWL();
	freeS();
	MPI_Finalize();
}
void recalibrate(){
	for(int i=0; i<D1BINS; i++)
		wllng[i] -= log(prof[i]);
}

double fit(double* x,double* y,int ndata)
{
	int i;
	double t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss;
	double  b=0.0;

	for (i=0;i<ndata;i++) {
		sx += x[i];
		sy += y[i];
	}
	ss=ndata;

	sxoss=sx/ss;

	for (i=0;i<ndata;i++) {
		t=x[i]-sxoss;
		st2 += t*t;
		b += t*y[i];
	}

	b /= st2;
	
	return b;

}



void update_prof(){
	int i,j, inv, cnt ,idx;
	double avgH, avgD;
	double* f = (double*) malloc(D1BINS*sizeof(double));

	for(i=0;i<D1BINS;i++){
		if(mask[i] == 1)
			f[i] = wlHu[i]/(wlHd[i]+wlHu[i]);
	}
	inv = (D1BINS-1)/Nprof;
	double* d = (double*) malloc(inv*sizeof(double));
	double df, maxD =-1;
	double x[Nprof], y[Nprof];
	idx = 0;
	for(i = 0; i< inv; i++){
		avgH = 0.0;
		cnt=0;
		while(cnt < Nprof){
			if(mask[idx] == 1){
				x[cnt] = -2*N*N+idx*4;
				y[cnt] = f[idx];
				avgH += (wlHu[idx] + wlHd[idx]);
				cnt++;
			}
			idx++;
		}
		avgH /= Nprof;
		df = fit(x, y, Nprof);
		d[i] = 1.0/sqrt(avgH*fabs(df));
		if(d[i] > maxD)
			maxD = d[i];
	}


/*	for(int iti = 0; iti< D1BINS/10; iti++){

		avgH = 0.0;
		if(iti >= 2 && iti <= (D1BINS/10-3) ){
				df = -0.2*f[(iti-2)*10] - 0.1*f[(iti-1)*10] + 0.1*f[(iti+1)*10] + 0.2*f[(iti+2)*10] ;
				for(i=(iti-2)*10; i < (iti+2)*10; i++)
					avgH += (wlHu[i] + wlHd[i]);
				avgH /= 40;

		}else if( iti == 0){
				df = -0.77143*f[iti*10] + 0.18571*f[(iti+1)*10] + 0.57143*f[(iti+2)*10] + 0.38571*f[(iti+3)*10] - 0.37143*f[(iti+4)*10];	
				for(i=iti*10; i < (iti+4)*10; i++)
					avgH += (wlHu[i] + wlHd[i]);
				avgH /= 40;

		}else if(iti == 1){
				df = -0.48571*f[(iti-1)*10] +  0.04286*f[iti*10] + 0.28571*f[(iti+1)*10] + 0.24286*f[(iti+2)*10] - 0.08571*f[(iti+3)*10];	
				for(i=(iti-1)*10; i < (iti+3)*10; i++)
					avgH += (wlHu[i] + wlHd[i]);
				avgH /= 40;
		}else if(iti == (D1BINS/10-2) ){
				df = 0.08571*f[(iti-3)*10] - 0.24286*f[(iti-2)*10] - 0.28571*f[(iti-1)*10] - 0.04286 *f[iti*10] + 0.48571*f[(iti+1)*10];	
				for(i=(iti-3)*10; i < (iti+1)*10; i++)
					avgH += (wlHu[i] + wlHd[i]);
				avgH /= 40;
		}else if(iti == (D1BINS/10-1) ){
				df = 0.37143*f[(iti-4)*10] - 0.38571*f[(iti-3)*10] - 0.57143*f[(iti-2)*10] - 0.18571 *f[(iti-1)*10] + 0.77143*f[iti*10];	
				for(i=(iti-4)*10; i < iti*10; i++)
					avgH += (wlHu[i] + wlHd[i]);
				avgH /= 40;
		}else{
			exit(77);
		}	

		df = fabs(df)*invdWLD1*10;
		prof[iti] = 1/sqrt(avgH*df);
		if(prof[iti] > maxD)
			maxD = prof[iti];
	}*/
		char s[512];
		FILE* fp;
	if(myrank == 0){

		sprintf(s,"diff%03d.dat",numf-1);
		fp = fopen(s, "w");
		for(i =0;i < inv; i++){
			d[i] /= maxD;
			fprintf(fp,"%g\t%g\n", (i+0.5)*Nprof*4-2*N*N, d[i]);
		}
		fflush(fp);
		fclose(fp);
	}

	idx = 0;
	for(i = 0; i< inv; i++){
		cnt=0;
		while(cnt < Nprof){
			if(mask[idx] == 1){
				prof[idx] = d[i];
				cnt++;
			}
			idx++;
		}
			
	}


	if(myrank == 0){
		sprintf(s, "current%03d.dat", numf-1);
		fp = fopen(s, "w");
		for( i =0;i <D1BINS;i++)
			if(mask[i] == 1)
				fprintf(fp,"%g\t%g\t%g\t%g\n", -2.0*N*N+i*4, f[i], prof[i], (wlHu[i]+wlHd[i]));
		fflush(fp);
		fclose(fp);
	}

	free(f);
	free(d);
}





