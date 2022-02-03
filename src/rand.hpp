#ifndef RAND_H
#define RAND_H

#ifdef DEFINE_GLOBALS
#define GLOBAL
#else
#define GLOBAL extern
#endif 
/*structure for storage of seeds:*/
struct seed_type 
{ 
unsigned int i;
unsigned int j;
unsigned int k;
};

GLOBAL int RSEED;

unsigned int kiss(void);
void gettimeseed(void);
void shelltimeseed(unsigned int tseed);
double randd1(void);
void gaussian(double m, double sigma, double *x1, double *x2);

#endif



