#pragma once

// Legacy KISS RNG + Gaussian helper. Used alongside the Mersenne Twister in
// random.h. The two share no state and are independent generators; KISS is
// used for the Wang–Landau acceptance tests (randd1) and the shuffle, while
// the MT is used for replica-exchange / Metropolis acceptance inside pt.cc.

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
