#ifndef WANGLANDAU_H
#define WANGLANDAU_H

#include "main.hpp"
#ifdef DEFINE_GLOBALS
#define GLOBAL
#else
#define GLOBAL extern
#endif 

#define PERW 0.8
#define LAMDA 1.0
#define KAPA  1000


GLOBAL double dWLD1;  //Bin Widths for Primary sampling direction
GLOBAL double WLD1max,WLD1min;  //sampling boundaries for PRIMARY (Energy) direction
GLOBAL int numf;  //Used to label DOS for each modification factor
GLOBAL int TotalSweeps,IterSweeps;  //Number of sweeps (total and per iteration)
GLOBAL int nummoves;

GLOBAL double Flatness;
GLOBAL int NumSweepsFlat;
GLOBAL double ModFactorInit;
GLOBAL double IterationFactor;
GLOBAL double ModFactorFinal;
GLOBAL int ProductionRun;
GLOBAL int ProductionBinSamps;

GLOBAL double numbelow_flat;  //gives the number of bins below the flatness criteria

//Parameters which are found from input parameters
GLOBAL int D1BINS;  //Number of Bins in the DOS and Histogram
GLOBAL double invdWLD1;  //Stores inverse bin widths

//For 2D WL simulations
GLOBAL double *wlH;  //the Wang-Landau accumulated histogram
//GLOBAL int *wlHd, *wlHi;  //the Wang-Landau accumulated histogram
GLOBAL unsigned short *wlHd, *wlHi;  //the Wang-Landau accumulated histogram
GLOBAL double *wllng, *wllng_prior, *wllngd, *wllngi;  //the Wang-Landau natural log of the density of states
GLOBAL double lnwlf;  //natural log of the Wang-Landau update factor, f
GLOBAL int *mask;
GLOBAL int *acceptrot, *attemptrot;

//Associated with the production run
GLOBAL int GRBINS;
GLOBAL double dGR;
GLOBAL double GRmax,GRmin;

GLOBAL int LOWESTE;
GLOBAL bool list[1000];

void initWL(void); //initializes for WL
void freeWL();
void sweepWL(int sweeeps, SamplingMode mode);//runs sweeps WL attempts
void wlhybrid();
void resetWL(void);//resets H histogram array
double flatWL(SamplingMode); //returns H_min/H_avg

void write_DOS_H(void);
void read_DOS_H(void);
void readg(void);
void readmask(void);
void writemask(void);
void write_restart(void);
void read_restart(void);
void global_update(int, double);
//Wang-Landau routines
int WangLandau(double Ei, double Ef, SamplingMode mode);//, double Enbi, double Enbf);

#endif

