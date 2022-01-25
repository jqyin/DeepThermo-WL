#ifndef  WATER_H
#define WATER_H

#define TIP4P

#ifdef SPCE
#define Qo  -0.8476   
#define Qh  0.4238
#define Angle 109.47
#define Alj  2633409.6 // 629.4*4.184*1000  kJoules/mol
#define Clj   2617.092 //625.5* 4.184   kJoules/mol
#define Roh  1.0
#endif

#ifdef TIP3P
#define Qo  -0.834  
#define Qh  0.417
#define Angle 104.52
#define Alj  2435088 // 582*4.184*1000  kJoules/mol
#define Clj  2489.48 //595* 4.184   kJoules/mol
#define Roh  0.9572
#endif

#ifdef TIP4P
//TIP4P
#define Qm  -1.04   
#define Qh  0.52
#define Angle 104.52
#define Alj  2510400 // 600*4.184*1000  kJoules/mol
#define Clj  2552.24 //610* 4.184   kJoules/mol
#define Roh  0.9572
#define Rom 0.15
#define Roc  0.06555840
#endif

#define Rconstrain 3.0
#define Pres 0.1   //  =  T * 1.0  atm 
#define Lmax 0.4


#define Mo 15.999
#define Mh 1.0079
#define ke2 1389.35  // e*e/(4*pi*epsilong)*Na  kJoules/mol 
#define Rg 0.008314472 // Na * Kb/1000  (gas constant)/1000   kJoules/K/mol

#define ProbD 0.5
#define Pi 3.14159265

int N;  //number of water molecure;
int NBINTERACTION;  //Selects the non-bonded interaction (1 = LJ potential, 2 = Quasi LJ potential)
double D,DD,D0,DD0;//D is for pivot angle (1.0 is 2pi), DC is for crank angle, and DD is for diffussion between +/-DD in each dir
double* Dr, *Df;
int* acceptdiff,*attemptdiff;
int* acceptrot,*attemptrot;
int* actdiff, *attdiff;
int* actrot, *attrot;
double** Ep, ** Epold; 

double invN,rootinvN; //inverse parameters, 1/N, sqrt(1/N)
double currEtot; //current total energy, total bonbonded energy, and temperature
double TTi,TTf,dTT;


double T, Lc0, Lca0;
double* Lc, *Lca;
int* acceptV, *attemptV;
int* actV, * attV;

struct Water{
	double Ox;
	double Oy;
	double Oz;

	double H1x;
	double H1y;
	double H1z;
	
	double H2x;
	double H2y;
	double H2z;

	double cx;
	double cy;
	double cz;

	
};
struct Water* Wmol; 
struct Water* Wold;
void ini_water(int, const char*);
double Etot(const struct Water*);
double Eupdate(int i,const struct Water* iWmol);
int constrain_ini(const struct Water* Wclus, int n);
int constrain(int n, int flag);
void DiffandRot();
void Diff();
void Rot();
void Vol();

void reset_XYZ(); 
void thermoqs(void);
void initialize(const char*);
void write_mol2(int frame,const struct Water*);
void writeinc(int frame, double energy);

void sphere(double ox, double oy, double oz, double * npos);
void cone(double ox, double oy, double oz, double h1x, double h1y, double h1z, double* npos);


void count_Hbond(void);


void orderqs(void);
void grcalculate(int pbin);
double coredensity(void);
double Rgyr2(void);
double EEdist(void);


#endif
