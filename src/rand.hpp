/**
 * @file rand.hpp
 * @brief KISS RNG and Gaussian helper used by the WL acceptance path.
 *
 * KISS supplies the pseudo-random stream consumed by ::randd1 in the
 * Wang–Landau acceptance test and by the lattice shuffle. It is
 * independent of the `std::mt19937` engine in pt.cc which drives
 * replica-exchange acceptance.
 */

#pragma once

/**
 * @brief KISS RNG seed triple.
 *
 * Internal state of the KISS generator. Not part of the public API;
 * exposed only because the ::kiss / ::shelltimeseed wrappers rely on
 * the structure layout.
 */
struct seed_type {
    unsigned int i;
    unsigned int j;
    unsigned int k;
};

/// @brief Last seed value installed (used for run.dat metadata).
extern int RSEED;

/**
 * @brief Advance the KISS generator and return its next 32-bit value.
 */
unsigned int kiss();

/**
 * @brief Seed KISS from the wall clock and `RSEED`.
 *
 * Used when no user-supplied seed is given. Production runs prefer
 * ::shelltimeseed for reproducibility.
 */
void gettimeseed();

/**
 * @brief Seed KISS deterministically from a single user-supplied integer.
 *
 * @param tseed Seed value (typically `argv[2]` of `hea-wl`).
 */
void shelltimeseed(unsigned int tseed);

/**
 * @brief Uniform random double in (0, 1).
 *
 * Excludes the endpoints (0 and 1 are rejected and re-drawn).
 */
double randd1();

/**
 * @brief Box–Muller Gaussian sampler.
 *
 * @param m Mean of the output Gaussian.
 * @param sigma Standard deviation.
 * @param[out] x1 First sample.
 * @param[out] x2 Second sample.
 */
void gaussian(double m, double sigma, double* x1, double* x2);
