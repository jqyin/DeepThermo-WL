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

	//Error message if the number of arguments is incorrect
	if (argc != 3 && myrank == 0) ErrorMsg(0, "");

	//Reads Input Parameters
	ReadInput(argv[1]);
  
  	//Initializes the System
	ini_sys();

	//Initialize the Wang-Landau sampling parameters
	initWL();

	while(  TotalSweeps < (int)MSAMPS  )
	{		
			for(i=0;i<N*N*N;++i){
				Rot();
			} 
			if( (TotalSweeps % (int)MSEP ) == 0)
			{
				MPI_Allreduce(wllngi, wllngd, D1BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				MPI_Allreduce(wlHi, wlHd, D1BINS, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			        for(i=0;i<D1BINS;i++){
					wlH[i] += wlHd[i];
					wllng[i] += wllngd[i];
				}	
				memset(wllngi, 0, D1BINS*sizeof(double));
				memset(wlHi, 0, D1BINS*sizeof(int));
			};
			TotalSweeps++;

			if(TotalSweeps % (int)MDROP == 0){
				i=0;
				while(wlH[i] == 0) i++;
				if(LOWESTE == i){
					if(lnwlf > 1e-9)
						lnwlf=lnwlf/IterationFactor;
					resetWL();
					for(i=0; i<10;i++)
						list[i] = false;
				}else if (LOWESTE > i)
					LOWESTE = i;

				if(myrank == 0)
					printf("TotalSweeps = %d, lnwlf = %g, Low = %d\n", TotalSweeps, lnwlf, LOWESTE);
			}

	};
	write_DOS_H();
	freeWL();
	freePT();
	MPI_Finalize();
}



