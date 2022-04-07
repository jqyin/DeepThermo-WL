#ifndef WANGLANDAU_H
#define WANGLANDAU_H

#include "main.hpp"
#ifdef DEFINE_GLOBALS
#define GLOBAL
#else
#define GLOBAL extern
#endif 

extern MPIState mpiState;
extern PTState ptState; 
extern WLState wlState; 
extern AlloyState alloyState; 

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

