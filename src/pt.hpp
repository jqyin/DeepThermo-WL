#pragma once

#include "sim_context.hpp"

constexpr int CHPT_STEPS = 100;
constexpr int TIMER = 2;

void ini_T(double Ti, double Tf, int nT);
void ini_sys();
void parallel_tempering(int nT, double DROPI, double SAMPS, double SEP, int i, SamplingMode m);
void swap(bool);
void freePT();  // retained for now; its body empties in favour of RAII in Task 7.
int Metropolis(double Ei, double Ef);
void read_state();
void write_state();
void mchybrid(SamplingMode);
