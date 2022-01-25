#ifndef  MD_TRIAL
#define MD_TRIAL

//scale

#define Ul 3.154   // A , length scale= (A/C)^1/6
//#define Um 2.987*10^23 // gram, mass scale= mass of a water molecule
#define Ue 0.155 // k cal /mole, energy scale
#define b 183.5  //  e^2 / Ue/ Ul
#define Ut 78.2 // K, temperature scale

#define  DT 0.00005 // time step, in unit of 1.66*10^(-12) s,  for regular MD dt=0.0005;
#define  Nt 5 
#define NtMax 100 //maximum time steps

#define Ixx 0.0098
#define Iyy 0.0034
#define Izz 0.0064

struct Quat{
	double q0;
	double q1;
	double q2;
	double q3;
};

struct Vec{
	double x;
	double y;
	double z;
};

struct Quat* quat;
struct Vec *F, *Tr, *V, *J;
struct Water* WmolOld; 

int* attemptmd;
int* acceptmd;
int* actmd, *attmd;


void ini_MD();
void ini_quat();
void FT();
void MD();
int constrain_MD();
double ini_p();
double gaussrand();
#endif
