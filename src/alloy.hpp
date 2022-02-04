#ifndef  HEISENBERG_H
#define HEISENBERG_H


#ifdef DEFINE_GLOBALS
#define GLOBAL
#else
#define GLOBAL extern
#endif 

#define SH 8   // number of shells
//#define Latt_Const 6.73 // A
#define T_scale 0.08618 // kT -> mev


GLOBAL bool* cluster;
GLOBAL double* Mavg, *Mavg2, *Mavg4;

 GLOBAL int N, N_2, N_3; //linear dimenstion;
 GLOBAL double invN;
 GLOBAL double currEtot;

GLOBAL int attd;
GLOBAL int accd;
 GLOBAL double D, DD; // step length;
GLOBAL double* S;
GLOBAL int* Type; 
GLOBAL double J[5][5][15];
GLOBAL int NS[10];

 GLOBAL double TTi,TTf,dTT;
 GLOBAL int NBINTERACTION; //dummy
GLOBAL char* atom;

// GLOBAL int* acceptrot,*attemptrot;

void initialize(int state);
void ini_conf(int state);
void ini_alloy();
void	ini_coupling();
double Etot();
inline void noffset(int i, int j, int k, int offi, int offj, int offk, int*nn, int cnt);
inline void  neighbor(int i, int j, int k, int* nn);
double Esite(int i, int j, int k);

void wolff(int i, int j, int k, double rx, double ry, double rz);
void Rot();
void Vol();
void Mag(double*);
void write_pos();
void write_mol2(int frame);



// p. b. c.
//#define BC(I) ( (I>=N)?(I-N):( (I<0)?(I+N):I ) )  
//#define BCb(I) ( (I>=(N-1) )?(I-(N-1)):( (I<0)?(I+(N-1)):I ) )  
#define Pi 3.14159265

#endif
