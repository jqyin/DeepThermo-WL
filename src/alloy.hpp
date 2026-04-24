#pragma once

#include "sim_context.hpp"

void initialize();
void ini_conf();
void ini_apos();
void ini_alloy(int state);
void ini_coupling();
double Eatom(int);
double Etot();
void neighbor(int i, int j, int k, int* nn);
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
