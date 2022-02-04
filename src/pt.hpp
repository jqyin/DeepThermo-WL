
#ifndef PT_H
#define PT_H

#define ORDER 1
#define DISORDER 2

#ifdef DEFINE_GLOBALS
#define GLOBAL
#else
#define GLOBAL extern
#endif 

//Input File Metropolis Variables
GLOBAL int MetropolisSampling;  //Turns Metropolis Sampling On (1) and Off (0) in the input file
GLOBAL double MTi;		//Initial Temperature used in simple Seq() loop
GLOBAL double MTf;		//Final Tempertaure used in simple Seq() loop
GLOBAL double MdT;		//Tempertaure Increment used in simple Seq() loop
GLOBAL double MSAMPS;	//Number of Samples taking in the Seq() function
GLOBAL double MSEP;		//Number of Point seperating actual taken data
GLOBAL double MDROP;	//Number of MC Steps dropped before each temperature run


GLOBAL double* T;
GLOBAL double pE, pT;
GLOBAL int att;
GLOBAL int acc;

void ini_T(double Ti, double Tf, int nT);
void ini_sys();
void mchybrid();
void parallel_tempering(int nT,double DROPI,double SAMPS, double SEP, int i);
void swap(bool);
void freePT();
int Metropolis(double Ei, double Ef);
void read_state();
void write_state();

#endif
