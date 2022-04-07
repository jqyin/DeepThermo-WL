#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<memory.h>
#include<assert.h>
#include <chrono>
#include "mpi.h"

#include "rand.hpp"
#include "wanglandau.hpp"
#include "alloy.hpp"
#include "pt.hpp"

double flatWL(SamplingMode mode)
{
	int i,k;
	double Hi,avg,min;
  
	k=0;  //number of sampled bins, used to calculate average
	min=1.0e300;  //the minimum sampled bin in the histogram wlH[][]
	avg=0.0;  //average height of the histogram wlH[][]
	wlState.numbelow_flat=0.0;  //number of bins below the flatness criteria
  
	//This loop finds the flatness and the number of unsampled bins
	for(i=0;i<wlState.D1BINS;++i){
		Hi = wlState.wlH[i];
		if(wlState.mask[i]==1)
		{			
			//minimum of the histogram
			if( (Hi<min) )
			{
				min=Hi;
			};
			//average of the histogram
			if( (Hi>0.0) )
			{
				avg+=Hi;
				k+=1;
			};	
		};
  	}
	avg/=1.0*k;  //Calculates the average height of the histogram (for  wlH > 0.0)
	if(k==0)
		avg=1.0;  //This keeps from dividing by zero
  
	wlState.numbelow_flat = 0.0;  
	//returns the current flatenss of the histogram		
	return min/avg;
  
}


//reset accumulated histogram
void resetWL()
{
	int i;
	for(i=0;i<wlState.D1BINS;++i){
		wlState.wlH[i]=0;
	}
}

//reads in a mask.dat file
void readmask(void)
{
  int i,k,tmp;
  FILE *ifp;
  
  ifp=fopen("mask.dat","r");
  
  for(k=0;k<wlState.D1BINS;++k)
    {
      fscanf(ifp,"%d\t%d\n",&i,&tmp);
      wlState.mask[i]=tmp;
    };
}

//reads in a mask.dat file
void writemask(void)
{
  int i,k,tmp;
  FILE *ofp;
  
  ofp=fopen("mask.dat","w");
  
  for(i=0;i<wlState.D1BINS;++i)
      fprintf(ofp,"%g\t%d\n",i/wlState.invdWLD1+wlState.WLD1min,wlState.mask[i]);
}

//reads in the g.dat file
void readg(void)
{
	int i;
	FILE *ifp;
	char s1[512];
	double tmp;
  
	ifp=fopen("g.dat","r");
	if(ifp != NULL){
		//Reading in the DOS and histogram from the restart function
		fscanf(ifp,"#%d\t%lg\t%d\t%d\t%lg\n",&wlState.numf,&wlState.ModFactorInit,&wlState.TotalSweeps,&wlState.IterSweeps, &tmp);
		printf("restart from: numf %d  lnwlf %g \n", wlState.numf, wlState.ModFactorInit);
		for(i=0;i<wlState.D1BINS;++i)
    		{
			fscanf(ifp,"%lg\t%lg\t%lg\t%lg \n",&tmp,&(wlState.wllng[i]),&tmp, &tmp );
    		};
		fclose(ifp);
	}
  
}

//writes out the DOS and the histogram
void write_DOS_H(void)
{
	int i;
	FILE *ofp;
	char s1[512];
	double maxg=-1.0e300;
	double ming=1.0e300;
  
	//uncomment the following line to output in log_10
	sprintf(s1,"DOS_H_iter%03d.dat",wlState.numf);
	ofp=fopen(s1,"w");
  
	//Find the minimum of the DOS
	for(i=0;i<wlState.D1BINS;++i)
	{
		if((wlState.wllng[i]>maxg) && (wlState.mask[i]==1))
			maxg=wlState.wllng[i];
		if((wlState.wllng[i]<ming) && (wlState.mask[i]==1))
			ming=wlState.wllng[i];
	};
  
	//Label each iteration with the mod. factor, number of sweeps, and flatness
	fprintf(ofp,"#%d\t%g\t%d\t%d\t%g\n",wlState.numf,wlState.lnwlf,wlState.TotalSweeps,wlState.IterSweeps,flatWL(WLdos)); //1.0*acceptdiff/attemptdiff, 1.0*acceptrot/attemptrot  );
	for(i=0;i<wlState.D1BINS;++i)
	{
		//Printing Out 1D Data
		fprintf(ofp,"%.8f\t%18.10e\t%18.10e\t%g\n",i/wlState.invdWLD1+wlState.WLD1min,(wlState.wllng[i]-maxg),1.0*wlState.wlH[i], 1.0*wlState.acceptrot[i]/wlState.attemptrot[i]);
	};
  
	fflush(ofp);
	fclose(ofp);
}

//Write out the normalized DOS //Currently only works for smaller chain sizes, this is because of underflow
void initWL(void)
{
	int i,j,k;
	int n;
	double Ei;


	wlState.TotalSweeps=1;		//Total number of sweeps
	wlState.lnwlf=wlState.ModFactorInit;

	//Primary Binning Direction for WL Simulation
	wlState.D1BINS = (int)((wlState.WLD1max - wlState.WLD1min)/(wlState.dWLD1*alloyState.invN));
	//Inverse Bin Width
	wlState.invdWLD1=1.0/(wlState.dWLD1*alloyState.invN);
	//correct the int cut off;
     	wlState.WLD1max = wlState.WLD1min +  wlState.D1BINS/wlState.invdWLD1;   

	wlState.wlH = (double*)malloc( wlState.D1BINS * sizeof(double) );
	wlState.wlHd = (unsigned short*)malloc(wlState.D1BINS * sizeof(unsigned short) );
	wlState.wlHi = (unsigned short*)malloc( wlState.D1BINS * sizeof(unsigned short) );
	wlState.wllng = (double*)malloc( wlState.D1BINS * sizeof(double ) );
	wlState.wllng_prior = (double*)malloc( wlState.D1BINS * sizeof(double ) );
	wlState.wllngd = (double*)malloc( wlState.D1BINS * sizeof(double ) );
	wlState.wllngi = (double*)malloc( wlState.D1BINS * sizeof(double ) );
	wlState.mask = (int*)malloc( wlState.D1BINS * sizeof(int ) );

	alloyState.op = (double*)malloc( wlState.D1BINS * sizeof(double ) );
	alloyState.op2 = (double*)malloc( wlState.D1BINS * sizeof(double ) );

	wlState.attemptrot = (int *)malloc(wlState.D1BINS*sizeof(int));
	wlState.acceptrot = (int *)malloc(wlState.D1BINS*sizeof(int));

	//Check to see if memory was allocated properly
	if ( (wlState.wlH == NULL) || (wlState.wllng == NULL) || (wlState.mask == NULL) )
    	{
        	fprintf(stderr,"\nFailure to allocate memory for 'Wang-Landau 2D Arrays'.  See 'initWL()'.\n");
       	 	exit(1);
    	};

	//Initialize Arrays
	for(i=0;i<wlState.D1BINS;++i)
	{
		wlState.wllng[i]=0.0;
		wlState.wllng_prior[i]=0.0;
		wlState.wllngd[i]=0.0;
		wlState.wllngi[i]=0.0;
		wlState.wlH[i]=0;
		wlState.wlHd[i]=0;
		wlState.wlHi[i]=0;
		wlState.mask[i]=1;
		wlState.attemptrot[i]=0;
		wlState.acceptrot[i]=0;
		alloyState.op[i] = alloyState.op2[i]= 0;

	};
	for(i = 0 ; i<1000;i++)
		wlState.list[i] = false;

	//Run the standard MC routine until the configuration has energy within the WL simulation energy range
	ptState.pT=100;
	while( (alloyState.currEtot*alloyState.invN> wlState.WLD1max)  )   // || (currEtot*invN<WLD1min) )
    	{
		mchybrid(metropolis);
		fprintf(stderr,"rank%d Relaxing:  %g of %g \n", mpiState.myrank, alloyState.currEtot*alloyState.invN,wlState.WLD1max);
    	};

	if(alloyState.currEtot*alloyState.invN<wlState.WLD1min){ // extend ground energy;
		fprintf(stderr,"currEtot/N: %g  \n",alloyState.currEtot*alloyState.invN);
	};
	wlState.LOWESTE = 0; 
	
        alloyState.E_0 = int(0.2*wlState.D1BINS);
	//checkpoint
	readg();
	alloyState.E_0 = alloyState.E_0 - int((wlState.numf-1)*0.01*wlState.D1BINS);
	MPI_Barrier(MPI_COMM_WORLD);
}

void freeWL(){
	free(wlState.wlH);
	free(wlState.wllng);
	free(wlState.wllng_prior);
	free(wlState.wlHi);
	free(wlState.wllngi);
	free(wlState.wlHd);
	free(wlState.wllngd);
	free(wlState.mask);
	free(wlState.attemptrot);
	free(wlState.acceptrot);

}
//attempts to run one WL hybrid move per bin
void sweepWL(int nsweeps, SamplingMode mode)
{
	for(int i=0;i<nsweeps;++i)
   	{
		BondSwap(mode);
	}
}


void global_update(int nsweeps, double Flat){
	int i, cnt=0;
	double duration = 0.0;
	double tmp_flat = 0.0;
	double maxg=-1.0e300;
	double ming=1.0e300;
	FILE * fp, *fp_time;
	std::chrono::high_resolution_clock::time_point start,end;
 	if(wlState.numf == 2 && mpiState.myrank == 0){
		fp_time = fopen("infer.dat", "w");	
	}
	while(tmp_flat <= Flat)
	{		
		sweepWL(nsweeps, WLprior);
		if(wlState.numf == 2)
			start = std::chrono::high_resolution_clock::now();
		vae_update(WLprior);
		if(wlState.numf == 2){
			end = std::chrono::high_resolution_clock::now();	
			std::chrono::duration<double, std::milli> ms_double = end - start;
			duration += ms_double.count(); 
			cnt++;
		}
		
		wlState.IterSweeps+=1;
		wlState.TotalSweeps+=1; 
		if( (wlState.TotalSweeps % 1) == 0)
		{
			//MPI_Allreduce(wlHi, wlHd, D1BINS, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(wlState.wlHi, wlState.wlHd, wlState.D1BINS, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
			for(i=0;i<wlState.D1BINS;i++){
				wlState.wlH[i] += wlState.wlHd[i];
				wlState.wllng_prior[i] += wlState.wlHd[i]*wlState.lnwlf;    //wllngd[i];
			}	
			memset(wlState.wllngi, 0, wlState.D1BINS*sizeof(double));
			memset(wlState.wlHi, 0, wlState.D1BINS*sizeof(unsigned short));
			tmp_flat=flatWL(WLprior);
			if( (wlState.IterSweeps %100)==0 &&  mpiState.myrank == 0)
				write_DOS_H();
		}
	}
 	if(wlState.numf == 2 && mpiState.myrank == 0){
		fprintf(fp_time, "infers/ms: %f\n", 2.0*cnt/duration);
		fflush(fp_time);	
	}

	if(mpiState.myrank == 0)
		fp = fopen("global-update-prior.dat", "w");
	for(i=0;i<wlState.D1BINS;++i){
                wlState.wllng_prior[i] = wlState.wllng_prior[i]*wlState.lnwlf;
		wlState.wllng[i] += wlState.wllng_prior[i];	
		if(mpiState.myrank == 0)
			fprintf(fp,"%g\t%18.10e\n",i/wlState.invdWLD1+wlState.WLD1min, wlState.wllng_prior[i]);
	}
	if(mpiState.myrank == 0)
		fclose(fp);
        
	resetWL();
	
}


//routine to encapsulate the WangLandau algorithm, returns 1 if accepted and 0 if rejected
int WangLandau(double Ei, double Ef, SamplingMode mode) //, double Enbi, double Enbf)
{
	int iti;//indices of initial config
	int fti;//indices of final config
	double lngi,lngf;  //values for the density of states, intial and final
	double R;

	//Primary direction index
	iti=(int) ((Ei*alloyState.invN-wlState.WLD1min)*wlState.invdWLD1);
	fti=(int) ((Ef*alloyState.invN-wlState.WLD1min)*wlState.invdWLD1);

	//This statement simply prints out the lowest 20 energy configurations;
	if(fti < wlState.LOWESTE+10 && fti >= wlState.LOWESTE)
	{	
		if(! wlState.list[fti - wlState.LOWESTE]){
			//write_xyz(fti - LOWESTE);
			wlState.list[fti - wlState.LOWESTE] = true;
		}
	};

	wlState.attemptrot[iti] += 1;
    	assert ( iti>=0 && iti < wlState.D1BINS);   

	if((fti<0)||(fti>=wlState.D1BINS)||(wlState.mask[fti]==0))
    	{
			wlState.wlHi[iti]+=1;
			wlState.wllngi[iti]+= wlState.lnwlf;	

			return 0;
    	}
	else
    	{
		//inside of bounds
		if(mode == WLdos || mode == WLproduction){
			lngi=wlState.wllngi[iti] + wlState.wllng[iti];
			lngf=wlState.wllngi[fti] + wlState.wllng[fti];
		}
		else if(mode == WLprior){//prior
			lngi=wlState.wllngi[iti] + wlState.wllng[iti] + wlState.wllng_prior[iti];
			lngf=wlState.wllngi[fti] + wlState.wllng[fti] + wlState.wllng_prior[fti];
		}
      		R=exp(lngi-lngf);
		if(randd1()<R)
		{
			//accept
			if(mode == WLprior || mode == WLdos){
				wlState.acceptrot[iti] += 1;
	    			wlState.wlHi[fti]+=1;
				wlState.wllngi[fti]+= wlState.lnwlf;
			}else if(mode == WLproduction){
				OrderParameter(fti);
	    			wlState.wlHi[fti]+=1;
				wlState.wllngi[fti]+= wlState.lnwlf;	
			}	
			return 1;
		}
		else
		{
			//reject
			if(mode == WLprior || mode == WLdos){
				wlState.wlHi[iti]+=1;
				wlState.wllngi[iti]+= wlState.lnwlf;
			}else if(mode == WLproduction){
				OrderParameter(iti);
				wlState.wlHi[iti]+=1;
				wlState.wllngi[iti]+= wlState.lnwlf;	
			}	
			return 0;
		};
    	};

}


