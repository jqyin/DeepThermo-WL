#include"water.h"
#ifndef PT_H
#define PT_H
//Input File Metropolis Variables
int MetropolisSampling;  //Turns Metropolis Sampling On (1) and Off (0) in the input file
double MTi;		//Initial Temperature used in simple Seq() loop
double MTf;		//Final Tempertaure used in simple Seq() loop
double MdT;		//Tempertaure Increment used in simple Seq() loop
double MSAMPS;	//Number of Samples taking in the Seq() function
double MSEP;		//Number of Point seperating actual taken data
double MDROP;	//Number of MC Steps dropped before each temperature run


double* pT;
double* pE;
double* attemps;
double* accepts;
struct Water** pW;

void ini_T(double Ti, double Tf, int nT);
void ini_sys(int nT,const char*);
void parallel_tempering(int nT,double DROPI,double SAMPS, double SEP);
void swap(int nT);
void mcdiff(int i);
int Metropolis(double Ei, double Ef,int n );
#endif