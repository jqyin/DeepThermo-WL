#ifndef  HEISENBERG_H
#define HEISENBERG_H
#include <cstdint>
#ifdef DEFINE_GLOBALS
#define GLOBAL
#else
#define GLOBAL extern
#endif 
#include "model.hpp"
#define NE 4
#define SH 6   // number of shells
#define O_SH 1   // number of shells
#define Latt_Const 6.73 // A
#define T_scale 0.08618 // kT -> mev
#define E_scale 13605.662285 // Ry -> mev
#define MAX_NEIGHBORS 300

//static const char *  element[]={"HEA", "Fe"};

//static const char* element[]={"HEA", "Mo","Nb","Ta","Ti","W"};
//static const int base_energy[NE]={-8090,-7632,-31232,-1704,-32312};
//static const int encode[NE][NE-1]={{0,0,0,0},{1,0,0,0},
//				   {0,1,0,0},{0,0,1,0},
//				   {0,0,0,1}};

static const char* element[]={"HEA", "Mo","Nb","Ta","W"};
static const int base_energy[NE]={-8090,-7632,-31233,-32312};
static const int encode[NE][NE-1]={{0,0,0},{1,0,0},
				   {0,1,0},{0,0,1}};


GLOBAL heaModel models[NE];
GLOBAL bool* cluster;
GLOBAL double* Mavg;

GLOBAL int N, N_2, N_3; //linear dimenstion;
GLOBAL double invN;
GLOBAL double currEtot;

GLOBAL int attd;
GLOBAL int accd;
GLOBAL double D, DD; // step length;
GLOBAL short* Atom; 
//GLOBAL uint8_t* Atom; 
GLOBAL int** inputPos; 
GLOBAL double J[NE][NE][SH];
GLOBAL int NS[SH];
GLOBAL double Dist[SH];
GLOBAL int NT[NE];
GLOBAL int Nneighbors; // number of neighbors

GLOBAL double TTi,TTf,dTT;
GLOBAL int NBINTERACTION; //dummy
GLOBAL struct {
	double x;
	double y;
	double z;
}nlist[MAX_NEIGHBORS];

// GLOBAL int* acceptrot,*attemptrot;

void initialize();
void ini_conf();
void ini_apos();
void ini_alloy(int state);
void	ini_coupling();
double Eatom(int);
double Etot();
inline void noffset(int i, int j, int k, int offi, int offj, int offk, int*nn, int cnt);
 void  neighbor(int i, int j, int k, int* nn);
double Esite(int i, int j, int k);

void wolff(int i, int j, int k, double rx, double ry, double rz);
void Rot();
void Vol();
void O(double* op);
void write_pos();
void write_mol2(int frame);



// p. b. c.
//#define BC(I) ( (I>=N)?(I-N):( (I<0)?(I+N):I ) )  
//#define BCb(I) ( (I>=(N-1) )?(I-(N-1)):( (I<0)?(I+(N-1)):I ) )  
#define Pi 3.14159265

#endif
