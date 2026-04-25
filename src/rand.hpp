#pragma once

// KISS RNG + Gaussian helper. Used by the Wang–Landau acceptance tests
// (randd1) and the lattice shuffle. Independent from the std::mt19937
// in pt.cc that drives replica-exchange acceptance.

struct seed_type {
    unsigned int i;
    unsigned int j;
    unsigned int k;
};

extern int RSEED;

unsigned int kiss(void);
void gettimeseed(void);
void shelltimeseed(unsigned int tseed);
double randd1(void);
void gaussian(double m, double sigma, double* x1, double* x2);
