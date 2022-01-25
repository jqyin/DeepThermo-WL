#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#include "InputOutput.h"
#include "rand.h"
#include "water.h"
#include "pt.h"
#include "wanglandau.h"
#include "MD_trial.h" 

//#define NDEBUG

int nprocs;
int myrank;
void ini_walker();

int main(int argc, char *argv[])
{
 
	int i,j,m;
	int STARTED;  //Shows whether or not the restart has been applied
	FILE *ofp_run;
	char s[512];
	double tmp_flat;
	int nT;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	
	//ini random number;
	srand(atoi(argv[4])*myrank);
	shelltimeseed(rand());

	//Error message if the number of arguments is incorrect
	if (argc !=5 && argc!=6) ErrorMsg(0, "");

	//Reads Input Parameters
	ReadInput(argv[1]);

    //Initializes the System
	if(argc==5 )
		initialize( NULL);
	else
		initialize(argv[5]);


	//Metropolis Sampling 
	if(MetropolisSampling == 1)
	{
        nT = atoi(argv[2]);
		ini_T(MTi,MTf,nT);
		if(argc==5 )
			ini_sys(nT,NULL);
		else
			ini_sys(nT,argv[5]);

		parallel_tempering(nT,MDROP,MSAMPS,MSEP);
		exit(0);	
	}


	//if JA greater than zero, code is currently not ready
/*	if(JA > 0.0)
	{
		fprintf(stderr,"You have to adjust the energy calculations if JA > 0.0\n");
		exit(1);
	};
*/
	//Initialize the Wang-Landau sampling parameters
	initWL();
	// synchronize all the walker;
	 ini_walker();

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
		if(myrank == 0)
			fprintf(stderr,"ProductionRun is not equal to 1 or 0! Is set to: %d \n\n",ProductionRun);
		exit(1);
	};
	
	if(myrank == 0)	{
		ofp_run=fopen("run.dat","a");
		fprintf(ofp_run,"#seeds: %d,%d,%d\n",314159265,362436069,atoi(argv[4]));
	}
	//Initialization
	numf=1;  //Keeps track of the number of interations  
	STARTED = 0;
    
	//Main Wang-Landau sampling loop
	for(lnwlf=ModFactorInit;lnwlf>ModFactorFinal;lnwlf=lnwlf/IterationFactor)
    {

		IterSweeps=0;	//Keeps track of the number of sweeps per iteration
		resetWL();		//Resets the WL histgram      
		reset_XYZ();	//Resets monomer positions by subtracting away the center of mass
		tmp_flat=0.0;	//Resets the flatness for the loop below
		MPI_Barrier(MPI_COMM_WORLD);
#ifdef  GLOBAL
		if(myrank == 0)
			global_update();
		MPI_Bcast(wlH, D1BINS, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
		MPI_Bcast(wllng, D1BINS, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
#endif		
	  	  
		//Decides if the code should be restarted or start from scratch
		if( (atoi(argv[2]) == 1) && (STARTED == 0) )
		{
			read_restart();  //Loads in the previously simulated data from the last checkpoint	  
			STARTED = 1;
			tmp_flat=flatWL();  //Checks to see if the previous restart file already had a flat histogram
			
			if( (lnwlf <= ModFactorFinal) && (tmp_flat >= Flatness) )
			{
				if(myrank == 0)
					fprintf(stderr,"Run has already reached the final modification factor!\n\n");
				exit(1);
			};
		};
  
		//Writes out the restart file at the beginning of each iteration
		if(myrank==0)
			write_restart();

		//This while loop is executed until the histogram is sampled
		while(  tmp_flat <= Flatness  )
		{
			//Performs a number () of WL sweeps
			sweepWL( D1BINS );
			IterSweeps+=1;		
			TotalSweeps+=1;
			tmp_flat=flatWL();
	 		
			if( (IterSweeps % 100 ) == 0)
			{// after a long run, still not flat;
				MPI_Reduce(attemptdiff, attdiff, D1BINS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Reduce(attemptrot, attrot, D1BINS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Reduce(acceptdiff, actdiff, D1BINS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Reduce(acceptrot, actrot, D1BINS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);

				MPI_Reduce(acceptmd, actmd, D1BINS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Reduce(attemptmd, attmd, D1BINS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);

				MPI_Reduce(acceptV, actV, D1BINS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Reduce(attemptV, attV, D1BINS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);

				MPI_Reduce(Lc, Lca, D1BINS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
				if(myrank==0){
					write_DOS_H();
					write_restart();
				}
			};

		};
				MPI_Reduce(attemptdiff, attdiff, D1BINS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Reduce(attemptrot, attrot, D1BINS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Reduce(acceptdiff, actdiff, D1BINS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Reduce(acceptrot, actrot, D1BINS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);

				MPI_Reduce(acceptmd, actmd, D1BINS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Reduce(attemptmd, attmd, D1BINS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);

				MPI_Reduce(acceptV, actV, D1BINS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Reduce(attemptV, attV, D1BINS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Reduce(Lc, Lca, D1BINS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if(myrank==0){
			write_DOS_H();
			write_restart();
		}
		//thermoqs();
		if(myrank == 0){	
			fprintf(ofp_run,"%g\t%d\t%d\t%g\t%g\n",lnwlf,TotalSweeps*D1BINS,IterSweeps*D1BINS,tmp_flat,numbelow_flat);
			fflush(ofp_run);
		}
		numf=numf+1;  //Used in write_DOS_H() to lable each DOS
	};
	if(myrank == 0)
		fclose(ofp_run); 			

	//Prints out the final normalized density of states, thermodynamics, and order parameters
	//write_DOS_H();
	//write_restart();
	if(myrank==0){
		write_normDOS();
		thermoqs();
	} 
	MPI_Finalize();

	return 0;

}

void ini_walker(){
	int i,j;
	double *bufferx, *buffery,*bufferz;

	bufferx = malloc(sizeof(double)*3*N);
	buffery = malloc(sizeof(double)*3*N);
	bufferz = malloc(sizeof(double)*3*N);
	//ini to the same water conf.
	if(myrank == 0){

	//	MPI_Pack_size(1, MPI_DOUBLE,MPI_COMM_WORLD,&wsize);
	//	wsize *=6;
	//	buffer = malloc(wsize);

		for(i=0;i<N;i++){
			bufferx[3*i] = Wmol[i].Ox;
			bufferx[3*i+1] = Wmol[i].H1x;
			bufferx[3*i+2] = Wmol[i].H2x;

			buffery[3*i] = Wmol[i].Oy;
			buffery[3*i+1] = Wmol[i].H1y;
			buffery[3*i+2] = Wmol[i].H2y;

			bufferz[3*i] = Wmol[i].Oz;
			bufferz[3*i+1] = Wmol[i].H1z;
			bufferz[3*i+2] = Wmol[i].H2z;
		}
	}

	MPI_Bcast(bufferx, 3*N, MPI_DOUBLE,  0, MPI_COMM_WORLD);
	MPI_Bcast(buffery, 3*N, MPI_DOUBLE,  0, MPI_COMM_WORLD);
	MPI_Bcast(bufferz, 3*N, MPI_DOUBLE,  0, MPI_COMM_WORLD);
	

	if(myrank != 0){
		for(i=0;i<N;i++){
			 Wmol[i].Ox=bufferx[3*i] ;
			  Wmol[i].H1x=bufferx[3*i+1];
			 Wmol[i].H2x=bufferx[3*i+2] ;

			 Wmol[i].Oy=buffery[3*i] ;
			  Wmol[i].H1y=buffery[3*i+1];
			 Wmol[i].H2y= buffery[3*i+2];

			 Wmol[i].Oz=bufferz[3*i] ;
			Wmol[i].H1z=bufferz[3*i+1]  ;
			 Wmol[i].H2z=bufferz[3*i+2] ;

		}

	currEtot=Etot(Wmol);

	for(i=0;i<N;i++){
			Wmol[i].cx = (Wmol[i].Ox *Mo +(Wmol[i].H1x+Wmol[i].H2x)*Mh ) /(2*Mh+Mo) ;
			Wmol[i].cy = ( Wmol[i].Oy *Mo +(Wmol[i].H1y+Wmol[i].H2y)*Mh) /(2*Mh+Mo) ;
			Wmol[i].cz = (Wmol[i].Oz *Mo +(Wmol[i].H1z+Wmol[i].H2z)*Mh) /(2*Mh+Mo) ;
		}

	for(i=0;i<N;i++)
		for(j=0;j<N;j++)
			Epold[i][j] = Epold[j][i]  = Ep[i][j];

	}


	//write_mol2(myrank, Wmol);

	free(bufferx);
	free(buffery);
	free(bufferz);
}

