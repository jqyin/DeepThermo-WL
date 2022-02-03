#ifndef WANGLANDAU_H
#define WANGLANDAU_H

//#define PROF

#ifdef DEFINE_GLOBALS
#define GLOBAL
#else
#define GLOBAL extern
#endif 
//#define GROUND  // extend ground energy;

//#define GLOBAL_UPDATE // global update;
#define PERW 0.8
#define LAMDA 1.0
#define  KAPA  100


//for delay WL;
//#define DELAY
#define dwl 1


//Input Parameters
GLOBAL int label; 

GLOBAL long long  tMC;
GLOBAL int Ed[dwl];
GLOBAL double lnwd[dwl];

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
GLOBAL double *wlH, *wlHt, *wlHi, *prof;  //the Wang-Landau accumulated histogram
GLOBAL double *wlHu, *wlHd, *wlHui, *wlHdi;
GLOBAL double *wllng, *wllngt, *wllngi, *gexact;  //the Wang-Landau natural log of the density of states
GLOBAL double *HRg;  //Array used to keep track of radius gyration (eventually holds the average)
GLOBAL double *HEEdist;  //Array used to keep track of end to end distance
GLOBAL double *Hcore; 
GLOBAL double **gr;  
GLOBAL double lnwlf;  //natural log of the Wang-Landau update factor, f
GLOBAL int *mask;

//Associated with the production run
GLOBAL int GRBINS;
GLOBAL double dGR;
GLOBAL double GRmax,GRmin;

GLOBAL double LOWESTE;

void production_run(void);
void histfill_prun(int pbin);
void init_production_run(void);
void initWL(void); //initializes for WL
void freeWL();
void sweepWL(int sweeeps);//runs sweeps WL attempts
void resetWL(void);//resets H histogram array
double flatWL(void); //returns H_min/H_avg

void write_DOS_H(void);
void read_DOS_H(void);
void write_normDOS(void);
void write_mov_g(int frame);
void readg(void);
void writeEERgyr(void);
void readmask(void);
void writemask(void);
void write_restart(void);
void read_restart(void);
void read_D(void); 
void write_D(void);

//Wang-Landau routines
int WangLandau(double Ei, double Ef);//, double Enbi, double Enbf);

void wlhybrid(void);
void seqMC(void);
int Metropolis(double Ei, double Ef);
void global_update(void);
void extendg(int );
void polint(const  double* xa, const double* ya, int n, double x, double* y,double*dy);

void wlrandpivot(void);
void wlcrankshaft(void);
void wldiff(void);
void wlreptation(void);
void wlcutjoin(void);

#endif

