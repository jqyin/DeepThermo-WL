#ifndef  HEISENBERG_H
#define HEISENBERG_H


#ifdef DEFINE_GLOBALS
#define GLOBAL
#else
#define GLOBAL extern
#endif 

//#define VIR

#define Latt_Const 2.86 // A
#define Viron0 -1.5371890 // Viron(equibrium)
#define T_scale 1  //0.00008618 // kT -> ev

//#define Const_K 1.0
#define R_constrain 0.288
#define V_cutoff 1
#define N_cutoff 2 //cutoff box, N_cutoff <= N/2 for P. B. C.;
#define r_cutoff 4.743 // cutoff radius, r_cutoff <= N_cutoff;
#define N_J 2000 //number of entries for J;
/*#define sizeHistMag 2000
int* histMag[sizeHistMag];
int* histMag2[sizeHistMag];
int* histMag4[sizeHistMag];
int* cntMag;
*/

GLOBAL double* NN;
GLOBAL int* cntNN;

 GLOBAL double* Mavg, *Mavg2, *Mavg4;

 GLOBAL int N; //linear dimenstion;
 GLOBAL double Nl;
 GLOBAL double Lbc;  //max distance due to p.b.c.
 GLOBAL double Rc_min;
 GLOBAL int Nbc; 
 GLOBAL double invN;
 GLOBAL double currEtot;

 GLOBAL double D, DD; // step length;
 GLOBAL int** S;  // spin array;
 GLOBAL double**** Sb; // spin in body center;
 GLOBAL double**** pos; // spin postion;
 GLOBAL double**** posb; // central spin position;

 GLOBAL double* Jr;// J(r);
 GLOBAL double Rmin, Rmax, invRwidth;
 GLOBAL double**** vol;//magnitude of spin
 GLOBAL double**** volb;//magnitude of central spin

 GLOBAL double T, TTi,TTf,dTT;
 GLOBAL int NBINTERACTION; //dummy
//Input File Metropolis Variables
 GLOBAL int MetropolisSampling;  //Turns Metropolis Sampling On (1) and Off (0) in the input file
 GLOBAL double MTi;		//Initial Temperature used in simple Seq() loop
 GLOBAL double MTf;		//Final Tempertaure used in simple Seq() loop
 GLOBAL double MdT;		//Tempertaure Increment used in simple Seq() loop
 GLOBAL double MSAMPS;	//Number of Samples taking in the Seq() function
 GLOBAL double MSEP;		//Number of Point seperating actual taken data
 GLOBAL double MDROP;	//Number of MC Steps dropped before each temperature run

 GLOBAL int* acceptrot,*attemptrot;

void initialize(const char* filename);
void ini_conf(int Read_conf, const char* filename);
//void	ini_coupling();
void count_nn(int t, FILE*);
double Etot();
double Esite(int i, int );
double Esiteb(int i, int j, int k);
double Eplaqueta(int i, int j, int k);
double Eplaquetab(int i, int j, int k);
void Rot();
void Vol();
void Mag(int);
void Tavg();
void write_conf();
void write_mol2(int frame);
//void write_Mdist( double Mag[sizeHistMag],  double T);
//void write_raw_data();
void thermoqs();
void freeS();
double ****create_a_foo ( int max_x, int max_y, int max_r, int max_c );

//////////improve on model 
double Viron(double r);
double mlength(double v);
double Volume(int i, int j, int k);
double Volumeb(int i, int j, int k);

// p. b. c.
#define BC(I) ( (I>=N)?(I-N):( (I<0)?(I+N):I ) )  
//#define BCb(I) ( (I>=(N-1) )?(I-(N-1)):( (I<0)?(I+(N-1)):I ) )  
#define Pi 3.14159265

#endif
