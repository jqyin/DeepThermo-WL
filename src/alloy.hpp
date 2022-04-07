#ifndef  HEISENBERG_H
#define HEISENBERG_H
#include <cstdint>
#ifdef DEFINE_GLOBALS
#define GLOBAL
#else
#define GLOBAL extern
#endif 
#include "main.hpp"

extern MPIState mpiState; 
extern AlloyState alloyState; 

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
void BondSwap(SamplingMode);
void vae_update(SamplingMode);
void decode(float* z);
void walk(float* npos);
void Vol();
double L1();
void write_pos();
void write_xyz(int frame);
void thermoqs();
void ini_W();
void OrderParameter(int idx);

#endif
