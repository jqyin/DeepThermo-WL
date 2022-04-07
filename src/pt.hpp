#ifndef PT_H
#define PT_H
#include "main.hpp"
#ifdef DEFINE_GLOBALS
#define GLOBAL
#else
#define GLOBAL extern
#endif 

#define CHPT_STEPS 100
#define TIMER 2


extern MPIState mpiState; 
extern PTState ptState; 


void ini_T(double Ti, double Tf, int nT);
void ini_sys();
void parallel_tempering(int nT,double DROPI,double SAMPS, double SEP, int i, SamplingMode m);
void swap(bool);
void freePT();
int Metropolis(double Ei, double Ef);
void read_state();
void write_state();
void mchybrid(SamplingMode);

#endif
