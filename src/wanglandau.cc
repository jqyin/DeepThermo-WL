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

extern int myrank, nprocs;

double flatWL(SamplingMode mode)
{
	int i,k;
	double Hi,avg,min;
  
	k=0;  //number of sampled bins, used to calculate average
	min=1.0e300;  //the minimum sampled bin in the histogram wlH[][]
	avg=0.0;  //average height of the histogram wlH[][]
	numbelow_flat=0.0;  //number of bins below the flatness criteria
  
	//This loop finds the flatness and the number of unsampled bins
	for(i=0;i<D1BINS;++i){
		Hi = wlH[i];
		if(mask[i]==1)
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
  
	numbelow_flat = 0.0;  
	//returns the current flatenss of the histogram		
	return min/avg;
  
}


//reset accumulated histogram
void resetWL()
{
	int i;
	for(i=0;i<D1BINS;++i){
		wlH[i]=0;
	}
}

//reads in a mask.dat file
void readmask(void)
{
  int i,k,tmp;
  FILE *ifp;
  
  ifp=fopen("mask.dat","r");
  
  for(k=0;k<D1BINS;++k)
    {
      fscanf(ifp,"%d\t%d\n",&i,&tmp);
      mask[i]=tmp;
    };
}

//reads in a mask.dat file
void writemask(void)
{
  int i,k,tmp;
  FILE *ofp;
  
  ofp=fopen("mask.dat","w");
  
  for(i=0;i<D1BINS;++i)
      fprintf(ofp,"%g\t%d\n",i/invdWLD1+WLD1min,mask[i]);
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
		fscanf(ifp,"#%d\t%lg\t%d\t%d\t%lg\n",&numf,&ModFactorInit,&TotalSweeps,&IterSweeps, &tmp);
		printf("restart from: numf %d  lnwlf %g \n", numf, ModFactorInit);
		for(i=0;i<D1BINS;++i)
    		{
			fscanf(ifp,"%lg\t%lg\t%lg\t%lg \n",&tmp,&(wllng[i]),&tmp, &tmp );
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
	sprintf(s1,"DOS_H_iter%03d.dat",numf);
	ofp=fopen(s1,"w");
  
	//Find the minimum of the DOS
	for(i=0;i<D1BINS;++i)
	{
		if((wllng[i]>maxg) && (mask[i]==1))
			maxg=wllng[i];
		if((wllng[i]<ming) && (mask[i]==1))
			ming=wllng[i];
	};
  
	//Label each iteration with the mod. factor, number of sweeps, and flatness
	fprintf(ofp,"#%d\t%g\t%d\t%d\t%g\n",numf,lnwlf,TotalSweeps,IterSweeps,flatWL(WLdos)); //1.0*acceptdiff/attemptdiff, 1.0*acceptrot/attemptrot  );
	for(i=0;i<D1BINS;++i)
	{
		//Printing Out 1D Data
		fprintf(ofp,"%.8f\t%18.10e\t%18.10e\t%g\n",i/invdWLD1+WLD1min,(wllng[i]-maxg),1.0*wlH[i], 1.0*acceptrot[i]/attemptrot[i]);
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


	TotalSweeps=1;		//Total number of sweeps
	lnwlf=ModFactorInit;
  

	//Primary Binning Direction for WL Simulation
	D1BINS = (int)((WLD1max - WLD1min)/(dWLD1*invN));
	//Inverse Bin Width
	invdWLD1=1.0/(dWLD1*invN);
	std::cout << "D1BINS: " << D1BINS << std::endl;
	//correct the int cut off;
     	WLD1max = WLD1min +  D1BINS/invdWLD1;   

	wlH = (double*)malloc( D1BINS * sizeof(double) );
	wlHd = (unsigned short*)malloc( D1BINS * sizeof(unsigned short) );
	wlHi = (unsigned short*)malloc( D1BINS * sizeof(unsigned short) );
	wllng = (double*)malloc( D1BINS * sizeof(double ) );
	wllng_prior = (double*)malloc( D1BINS * sizeof(double ) );
	wllngd = (double*)malloc( D1BINS * sizeof(double ) );
	wllngi = (double*)malloc( D1BINS * sizeof(double ) );
	mask = (int*)malloc( D1BINS * sizeof(int ) );

	op = (double*)malloc( D1BINS * sizeof(double ) );
	op2 = (double*)malloc( D1BINS * sizeof(double ) );


	attemptrot = (int *)malloc(D1BINS*sizeof(int));
	acceptrot = (int *)malloc(D1BINS*sizeof(int));

	//Check to see if memory was allocated properly
	if ( (wlH == NULL) || (wllng == NULL) || (mask == NULL) )
    	{
        	fprintf(stderr,"\nFailure to allocate memory for 'Wang-Landau 2D Arrays'.  See 'initWL()'.\n");
       	 	exit(1);
    	};

	//Initialize Arrays
	for(i=0;i<D1BINS;++i)
	{
		wllng[i]=0.0;
		wllng_prior[i]=0.0;
		wllngd[i]=0.0;
		wllngi[i]=0.0;
		wlH[i]=0;
		wlHd[i]=0;
		wlHi[i]=0;
		mask[i]=1;
		attemptrot[i]=0;
		acceptrot[i]=0;
		op[i] = op2[i]= 0;

	};
	for(i = 0 ; i<1000;i++)
		list[i] = false;

	//Run the standard MC routine until the configuration has energy within the WL simulation energy range
	pT=100;
	while( (currEtot*invN> WLD1max)  )   // || (currEtot*invN<WLD1min) )
    	{
		mchybrid(metropolis);
		fprintf(stderr,"rank%d Relaxing:  %g of %g \n", myrank, currEtot*invN,WLD1max);
    	};

	if(currEtot*invN<WLD1min){ // extend ground energy;
		fprintf(stderr,"currEtot/N: %g  \n",currEtot*invN);
	};
	LOWESTE = 0; 
	
        E_0 = int(0.2*D1BINS);
	//checkpoint
	readg();
	E_0 = E_0 - int((numf-1)*0.01*D1BINS);
	MPI_Barrier(MPI_COMM_WORLD);
}

void freeWL(){
	int i;
	free(wlH);
	free(wllng);
	free(wllng_prior);
	free(wlHi);
	free(wllngi);
	free(wlHd);
	free(wllngd);
	free(mask);
	free(attemptrot);
	free(acceptrot);



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
 	if(numf == 2 && myrank == 0){
		fp_time = fopen("infer.dat", "w");	
	}
	while(tmp_flat <= Flat)
	{		
		sweepWL(nsweeps, WLprior);
		if(numf == 2)
			start = std::chrono::high_resolution_clock::now();
		vae_update(WLprior);
		if(numf == 2){
			end = std::chrono::high_resolution_clock::now();	
			std::chrono::duration<double, std::milli> ms_double = end - start;
			duration += ms_double.count(); 
			cnt++;
		}
		
		IterSweeps+=1;
		TotalSweeps+=1; 
		if( (TotalSweeps % 1) == 0)
		{
			//MPI_Allreduce(wlHi, wlHd, D1BINS, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(wlHi, wlHd, D1BINS, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
			for(i=0;i<D1BINS;i++){
				wlH[i] += wlHd[i];
				wllng_prior[i] += wlHd[i]*lnwlf;    //wllngd[i];
			}	
			memset(wllngi, 0, D1BINS*sizeof(double));
			memset(wlHi, 0, D1BINS*sizeof(unsigned short));
			tmp_flat=flatWL(WLprior);
			if( (IterSweeps %100)==0 &&  myrank == 0)
				write_DOS_H();
		}
	}
 	if(numf == 2 && myrank == 0){
		fprintf(fp_time, "infers/ms: %f\n", 2.0*cnt/duration);
		fflush(fp_time);	
	}

	for(i=0;i<D1BINS;++i)
	{
		if((wllng_prior[i]>maxg) && (mask[i]==1))
			maxg=wllng_prior[i];
		if((wllng[i]<ming) && (mask[i]==1))
			ming=wllng[i];
	}
	if(myrank == 0)
		fp = fopen("global-update-prior.dat", "w");
	for(i=0;i<D1BINS;++i){
                wllng_prior[i] = wllng_prior[i]*lnwlf;
		wllng[i] += wllng_prior[i];	
		if(myrank == 0)
			fprintf(fp,"%g\t%18.10e\n",i/invdWLD1+WLD1min, wllng_prior[i]);
	}
	if(myrank == 0)
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
	iti=(int) ((Ei*invN-WLD1min)*invdWLD1);
	fti=(int) ((Ef*invN-WLD1min)*invdWLD1);

	//This statement simply prints out the lowest 20 energy configurations;
	if(fti < LOWESTE+10 && fti >= LOWESTE)
	{	
		if(! list[fti - LOWESTE]){
			//write_xyz(fti - LOWESTE);
			list[fti - LOWESTE] = true;
		}
	};

	attemptrot[iti] += 1;
    	assert ( iti>=0 && iti < D1BINS);   

	if((fti<0)||(fti>=D1BINS)||(mask[fti]==0))
    	{
			wlHi[iti]+=1;
			wllngi[iti]+= lnwlf;	

			return 0;
    	}
	else
    	{
		//inside of bounds
		if(mode == WLdos || mode == WLproduction){
			lngi=wllngi[iti] + wllng[iti];
			lngf=wllngi[fti] + wllng[fti];
		}
		else if(mode == WLprior){//prior
			lngi=wllngi[iti] + wllng[iti] + wllng_prior[iti];
			lngf=wllngi[fti] + wllng[fti] + wllng_prior[fti];
		}
      		R=exp(lngi-lngf);
		if(randd1()<R)
		{
			//accept
			if(mode == WLprior || mode == WLdos){
				acceptrot[iti] += 1;
	    			wlHi[fti]+=1;
				wllngi[fti]+= lnwlf;
			}else if(mode == WLproduction){
				OrderParameter(fti);
	    			wlHi[fti]+=1;
				wllngi[fti]+= lnwlf;	
			}	
			return 1;
		}
		else
		{
			//reject
			if(mode == WLprior || mode == WLdos){
				wlHi[iti]+=1;
				wllngi[iti]+= lnwlf;
			}else if(mode == WLproduction){
				OrderParameter(iti);
				wlHi[iti]+=1;
				wllngi[iti]+= lnwlf;	
			}	
			return 0;
		};
    	};

}


