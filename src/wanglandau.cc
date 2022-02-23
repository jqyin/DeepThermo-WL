#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<memory.h>
#include<assert.h>

#include "rand.hpp"
#include "wanglandau.hpp"
#include "alloy.hpp"
#include "pt.hpp"

extern int myrank, nprocs;

double flatWL(void)
{
	int i,k;
	double Hi,avg,min;
  
	k=0;  //number of sampled bins, used to calculate average
	min=1.0e300;  //the minimum sampled bin in the histogram wlH[][]
	avg=0.0;  //average height of the histogram wlH[][]
	numbelow_flat=0.0;  //number of bins below the flatness criteria
  
	//This loop finds the flatness and the number of unsampled bins
	for(i=0;i<D1BINS;++i){
		Hi = wlH[i]*((i < E_0)? k_f : 1.0);
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
  
	//This loop counts the number of bins below the flatness criteria
  /*
    for(i=0;i<D1BINS;++i)
    for(j=0;j<D2BINS;++j)
    if(mask[i]==1)
    {
    if( wlH[i]<Flatness*avg )
    numbelow_flat+=1.0;
    };
    //Stores the percentage of states below the flatness criteria	
    numbelow_flat = (numbelow_flat)/(1.0*D1BINS*D2BINS);
  */
	numbelow_flat = 0.0;  //Turned off for now
	//returns the current flatenss of the histogram		
	return min/avg;
  
}


//reset accumulated histogram
void resetWL()
{
	int i;
 //   double maxg = -1e300;
	for(i=0;i<D1BINS;++i){
		wlH[i]=0;
//		if(wllng[i] > maxg)
//			maxg = wllng[i];
	}
//	for(i=0;i<D1BINS;++i)
//		wllng[i] -= maxg;
  
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
  
//	sprintf(s1,"DOS_H_iter%03d.dat",numf);
	ifp=fopen("g.dat","r");
	//Reading in the DOS and histogram from the restart function
	fscanf(ifp,"#%d\t%lg\t%d\t%d\t%lg\t%lg\n",&numf,&lnwlf,&TotalSweeps,&i, &tmp, &tmp);
	//fprintf(stderr,"#%d\t%lg\t%d\t%d\t%g\n",numf,lnwlf,TotalSweeps,IterSweeps,tmp);
	for(i=0;i<D1BINS;++i)
    {
		fscanf(ifp,"%lg\t%lg\t%lg\t%lg \n",&tmp,&(wllng[i]),&(wlH[i]),&tmp );
		//fprintf(stderr,"%g\t%18.10e\t%g\n",tmp,wllng[i],wlH[i]);
    };
	fclose(ifp);
  
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
	//tmp=log10(exp(1.0));
	sprintf(s1,"DOS_H_iter%03d.dat",numf);
	ofp=fopen(s1,"w");
	//ofp=fopen("g.dat","w");
  
	//Find the minimum of the DOS
	for(i=0;i<D1BINS;++i)
	{
		if((wllng[i]>maxg) && (mask[i]==1))
			maxg=wllng[i];
		if((wllng[i]<ming) && (mask[i]==1))
			ming=wllng[i];
	};
  
	//Label each iteration with the mod. factor, number of sweeps, and flatness
	fprintf(ofp,"#%d\t%g\t%d\t%d\t%g\n",numf,lnwlf,TotalSweeps,IterSweeps,flatWL()); //1.0*acceptdiff/attemptdiff, 1.0*acceptrot/attemptrot  );
	for(i=0;i<D1BINS;++i)
	{
		//Printing Out 1D Data
		//fprintf(ofp,"%g\t%g\t%g\n",i/invdWLD1+WLD1min,(wllng[i]-maxg),wlH[i]);
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

	//this value times Ebond-WLEbondmin, gives the bin number;
	//WLinvdEbond=(1.0*Lbond)/(WLEbondmax-WLEbondmin);
	//WLinvdEnonbond=(1.0*Lnonbond)/(WLEnonbondmax-WLEnonbondmin);

	//Dynamic Memory Allocation - this has to be done because an input file is used.
	//WL 2D Arrays - Histogram, Density of States, Mask, and Rawmask
	//allocate storage for an array of pointers
  
	wlH = (int*)malloc( D1BINS * sizeof(int) );
	wlHd = (unsigned short*)malloc( D1BINS * sizeof(unsigned short) );
	wlHi = (unsigned short*)malloc( D1BINS * sizeof(unsigned short) );
	wllng = (double*)malloc( D1BINS * sizeof(double ) );
	wllngd = (double*)malloc( D1BINS * sizeof(double ) );
	wllngi = (double*)malloc( D1BINS * sizeof(double ) );
	mask = (int*)malloc( D1BINS * sizeof(int ) );

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
		wllngd[i]=0.0;
		wllngi[i]=0.0;
		wlH[i]=0;
		wlHd[i]=0;
		wlHi[i]=0;
		mask[i]=1;
		attemptrot[i]=0;
		acceptrot[i]=0;

	};
	for(i = 0 ; i<1000;i++)
		list[i] = false;

	//Run the standard MC routine until the configuration has energy within the WL simulation energy range
	pT=100;
	while( (currEtot*invN> WLD1max)  )   // || (currEtot*invN<WLD1min) )
    	{
		mchybrid(1);
		fprintf(stderr,"rank%d Relaxing:  %g of %g \n", myrank, currEtot*invN,WLD1max);
    	};

	if(currEtot*invN<WLD1min){ // extend ground energy;
		fprintf(stderr,"currEtot/N: %g  \n",currEtot*invN);
		//exit(-1);
	};
	LOWESTE = 0; //(int) ((currEtot*invN-WLD1min)*invdWLD1);
}

void freeWL(){
	int i;
	free(wlH);
	free(wllng);
	free(wlHi);
	free(wllngi);
	free(wlHd);
	free(wllngd);
	free(mask);
	free(attemptrot);
	free(acceptrot);



}
//attempts to run one WL hybrid move per bin
void sweepWL(int sweeps, int mode)
{
	for(int i=0;i<sweeps;++i)
   	{
		//wlhybrid();
		Rot(mode);
	};
}


//excutes all of the combined MC moves, performs one step of each type of move, and one sweep of diffusion moves
void wlhybrid(void)
{
	int i,j,k;
	for(i=0;i<N*N*N;++i){
		Rot(1);
	} 
	
}

void global_update(){
	int i,j;
	double dg;
	FILE * fp;
	double w;
    	double sum =0.0;

	if(wllng[0] != 0.0) {

		w = (PERW*wllng[0]);     // assuming it's monotonically increase;
		fp = fopen("after_global_updata.dat", "w");
		for(i=0;i<D1BINS;++i){
			if(wllng[i] > w){
				dg = KAPA*lnwlf * exp( -LAMDA/(wllng[i]-w) );
				wllng[i] += dg;
			}
			fprintf(fp,"%g\t%18.10e\n",i/invdWLD1+WLD1min,wllng[i]);
		}
		fclose(fp);
	}else{ // first iteration;
		j = (int)((1-PERW)*D1BINS);
		for(i=j+1;i<D1BINS;++i){
			dg = KAPA*lnwlf * exp( -LAMDA/(i-j) );
			wllng[i] += dg;
		}
	
	}
//	KAPA = sqrt(KAPA);
}


//routine to encapsulate the WangLandau algorithm, returns 1 if accepted and 0 if rejected
int WangLandau(double Ei, double Ef) //, double Enbi, double Enbf)
{
	int iti;//indices of initial config
	int fti;//indices of final config
	double lngi,lngf;  //values for the density of states, intial and final
	double R;

	//Primary direction index
	iti=(int) ((Ei*invN-WLD1min)*invdWLD1);
	fti=(int) ((Ef*invN-WLD1min)*invdWLD1);

	//This statement simply prints out the lowest 20 energy configurations;
	if(fti < LOWESTE+300 && fti >= LOWESTE)
	{	
		if(! list[fti - LOWESTE]){
			//write_xyz(fti - LOWESTE);
			list[fti - LOWESTE] = true;
		}
		//fprintf(stderr,"New Lowest E config %g\n\n",LOWESTE);
	};

	attemptrot[iti] += 1;
    	assert ( iti>=0 && iti < D1BINS);   

	if((fti<0)||(fti>=D1BINS)||(mask[fti]==0))
    	{
		/*if(fti < 0 ){
			printf("Emin not low enough! Ei = %g, Ef=%g\n", Ei, Ef);
			exit(-1);
		}
		if(fti >=D1BINS){*/
			wlHi[iti]+=1;
			//wllngi[iti]+=lnwlf;	
			//wlHi[iti]+= (iti > E_0)? 1.0/(1.0+k_f): 1.0;
			wllngi[iti]+= (iti < E_0)? lnwlf+lngk_f: lnwlf;	
			return 0;
		//}

    	}
	else
    	{
		//inside of bounds

		lngi=wllngi[iti] + wllng[iti];
		lngf=wllngi[fti] + wllng[fti];
      		R=exp(lngi-lngf);
		if(randd1()<R)
		{
			//accept
			acceptrot[iti] += 1;
	    		wlHi[fti]+=1;
			//wllngi[fti]+=lnwlf;		
			//wlHi[fti]+= (fti > E_0)? 1.0/(1.0+k_f): 1.0;
			wllngi[fti]+= (fti < E_0)? lnwlf+lngk_f: lnwlf;	
			return 1;
		}
		else
		{
			//reject
			wlHi[iti]+=1;
			//wllngi[iti]+=lnwlf;	
			//wlHi[iti]+= (iti > E_0)? 1.0/(1.0+k_f): 1.0;
			wllngi[iti]+= (iti < E_0)? lnwlf+lngk_f: lnwlf;	
			return 0;
		};
    	};

}


