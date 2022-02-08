
#ifndef PT_H
#define PT_H

#ifdef DEFINE_GLOBALS
#define GLOBAL
#else
#define GLOBAL extern
#endif 

#define CHPT_STEPS 100
#define TIMER 2
//Input File Metropolis Variables
GLOBAL int MetropolisSampling;  //Turns Metropolis Sampling On (1) and Off (0) in the input file
GLOBAL int IRUN;  //Number of independent runs;
GLOBAL double MTi;		//Initial Temperature used in simple Seq() loop
GLOBAL double MTf;		//Final Tempertaure used in simple Seq() loop
GLOBAL double MdT;		//Tempertaure Increment used in simple Seq() loop
GLOBAL double MSAMPS;	//Number of Samples taking in the Seq() function
GLOBAL double MSEP;		//Number of Point seperating actual taken data
GLOBAL double MDROP;	//Number of MC Steps dropped before each temperature run
GLOBAL int iState;	//Initial state 0 for disorder 1 for order
GLOBAL int Restart;	// 0 for no 1 for yes


GLOBAL double* T;
GLOBAL double pT;
GLOBAL int att;
GLOBAL int acc;

void ini_T(double Ti, double Tf, int nT);
void ini_sys();
void parallel_tempering(int nT,double DROPI,double SAMPS, double SEP, int i, int m);
void swap(bool);
void freePT();
int Metropolis(double Ei, double Ef);
void read_state();
void write_state();
void mchybrid(int);

#endif
