#ifndef MAIN_H
#define MAIN_H

#include "model.hpp"

#ifdef CRUSHER
#define GPUperNode 8
#else
#define GPUperNode 6
#endif

//Alloy 
#define NE 4   // number of elements 
#define SH 6   // number of shells
#define T_scale 0.08618 // kT -> mev
#define E_scale 13605.69301 // Ry -> mev
#define MAX_NEIGHBORS 300

enum SamplingMode {metropolis, WLprior, WLdos, WLproduction};

struct MPIState{
	int nprocs;	//world size
	int myrank;	//world rank
};

struct PTState{
        int MetropolisSampling;  //Turns Metropolis Sampling On (1) and Off (0) in the input file
        double MTi;		//Initial Temperature used in simple Seq() loop
        double MTf;		//Final Tempertaure used in simple Seq() loop
        double MdT;		//Tempertaure Increment used in simple Seq() loop
        double MSAMPS;	//Number of Samples taking in the Seq() function
        double MSEP;		//Number of Point seperating actual taken data
        double MDROP;	//Number of MC Steps dropped before each temperature run
        int Restart;	//0 for no 1 for yes
        double* T;	//temperature list 
        double pT;	//temperature for the myrank
        int att;	//Number of attempted swap 
        int acc;	//Number of accepted swap 
};

struct WLState{
        double dWLD1;	//Bin Widths for Primary sampling direction
        double WLD1max,WLD1min;	//sampling boundaries for PRIMARY (Energy) direction
        int numf;	//Used to label DOS for each modification factor
        int TotalSweeps,IterSweeps;	//Number of sweeps (total and per iteration)

        double Flatness;	//flatness of histogram	
        double ModFactorInit;	//inital lnf
        double IterationFactor;	//reduce factor 
        double ModFactorFinal;	//final lnf
        int ProductionBinSamps;	//number of sampls for production 
        double numbelow_flat;  //gives the number of bins below the flatness criteria
        int D1BINS;	//Number of Bins in the DOS and Histogram
        double invdWLD1;	//Stores inverse bin widths

        double *wlH;	//the Wang-Landau accumulated histogram
        unsigned short *wlHd, *wlHi;	//the Wang-Landau accumulated histogram
        double *wllng, *wllng_prior, *wllngd, *wllngi;  //the Wang-Landau natural log of the density of states
        double lnwlf;  //natural log of the Wang-Landau update factor, f
        int *mask;	//mask for unaccessible E 
        int *acceptrot, *attemptrot; //number of attempted and accepted moves 

        int LOWESTE;	//cutoff E to print out configuration
        bool list[1000]; //list of print out configurations 

};

struct AlloyState{
	const char* element[5]={"HEA", "Mo","Nb","Ta","W"};
	// MoNbTaW
	const double reglin_intercept = -1.2702430255548436; //Ry
	// MoNbTaVW
	//const double reglin_intercept = -1.7414431589249322; //Ry
	// MoNbTaTiW
	//const double reglin_intercept = -1.5945920830411398; //Ry
	const double Z_R = 0.1;	//radius of random walk in z space
        int E_0;	//cutoff E for applying dl proposal

        int VAE_D, SHIFT, PAD;	//vae input dimension, shift, and padding
        float* inputConfig;	//input configuration 

        heaModel model[2];	//vae model for encoder and decoder
        double* op, *op2;	//vae order parameter 

        int N, N_2, N_3;	//dimenstion for simulation cell
        double invN;		//inverse N
        double currEtot;	//current energy 

        int attd;		//number of attempted move  
        int accd;		//number of accepted move
        short* Atom; 		//Atom configuration 
        short* Atomo; 		
        int* elist[NE]; 	//element list
        int** inputPos; 	//input positions 
        double J[NE][NE][SH];	//interaction coupling 
        int W[NE][NE][SH];	//probablity for atom pairs 
        int NS[SH];		//number of atoms in each shell
	double Dist[SH];	//shell distance
        int NT[NE];		//number of atoms for each element
        int Nneighbors; 	//number of neighbors

        double TTi,TTf,dTT;	//temperature range for thermodynamics 
        int NBINTERACTION;	 
        struct {
	double x;
	double y;
	double z;
	}nlist[MAX_NEIGHBORS];	//neighbor list 
};



extern MPIState mpiState; 
extern PTState ptState; 
extern WLState wlState; 
extern AlloyState alloyState; 

#endif
