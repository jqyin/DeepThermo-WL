#pragma once

#include "sim_context.hpp"

void initWL(void);
void sweepWL(int sweeps, SamplingMode mode);
void wlhybrid();
void resetWL(void);
double flatWL(SamplingMode);

void write_DOS_H(void);
void read_DOS_H(void);
void readg(void);
void readmask(void);
void writemask(void);
void write_restart(void);
void read_restart(void);
void global_update(int, double);
int WangLandau(double Ei, double Ef, SamplingMode mode);
