/**
 * @file wanglandau.hpp
 * @brief Wang–Landau histogram management and acceptance machinery.
 *
 * Implements the inner half of paper Algo. 1: ::WangLandau handles the
 * per-step transition and DOS update, ::sweepWL drives a fixed number
 * of attempts, ::flatWL evaluates the histogram-flatness convergence
 * criterion, and ::initWL allocates the per-bin arrays in ::wlState.
 */

#pragma once

#include "sim_context.hpp"

/**
 * @brief Allocate WL histograms and seed the energy window.
 *
 * Computes `D1BINS`, calls `WLState::resize`, sizes the order-parameter
 * accumulators on `alloyState`, and runs Metropolis sweeps until
 * `currEtot/N <= WLD1max` so the WL loop starts inside its energy
 * window.
 */
void initWL();

/**
 * @brief Run a fixed number of WL bond-swap attempts.
 * @param sweeps Number of `BondSwap` calls.
 * @param mode Acceptance regime (one of `WLprior`, `WLdos`, `WLproduction`).
 */
void sweepWL(int sweeps, SamplingMode mode);

void wlhybrid();

/**
 * @brief Zero out the working histogram between WL outer iterations.
 */
void resetWL();

/**
 * @brief Histogram flatness ratio (min / avg) for the convergence test.
 *
 * @param mode Sampling regime; affects which DOS array participates in the
 *             acceptance probability but not the flatness computation
 *             itself.
 * @return min(H[i] over masked bins) / mean(H[i] over masked bins with H>0).
 */
double flatWL(SamplingMode mode);

/**
 * @brief Write the per-iteration DOS + histogram to `DOS_H_iter<numf>.dat`.
 */
void write_DOS_H();

void read_DOS_H();

/**
 * @brief Read a previously-saved DOS from `g.dat` for restart.
 *
 * Silently no-ops if the file is missing.
 */
void readg();

/**
 * @brief Read a previously-saved bin mask from `mask.dat`.
 */
void readmask();

/**
 * @brief Write the bin mask to `mask.dat`.
 */
void writemask();

void write_restart();
void read_restart();

/**
 * @brief Run the prior-DOS estimation phase (paper §III-B "global update").
 *
 * @param sweeps Number of bond-swap attempts per pass.
 * @param flat_target Convergence threshold on the prior histogram.
 *
 * Iterates `sweepWL` + `vae_update` until the prior histogram passes
 * `flat_target`, then folds the prior into `wllng`.
 */
void global_update(int sweeps, double flat_target);

/**
 * @brief Append cumulative VAE global-move acceptance to vae.dat.
 *
 * Collective -- every rank must call it. Columns: total sweeps, lnwlf,
 * attempts, accepts, ratio. This is the number that says whether the VAE
 * surrogate is earning its inference cost against plain local moves.
 */
void report_vae(int sweeps, double lnwlf);

/**
 * @brief Wang–Landau acceptance test for one trial energy change.
 *
 * @param Ei Initial total energy.
 * @param Ef Final total energy after the proposed move.
 * @param mode Acceptance regime.
 * @return 1 if accepted, 0 if rejected.
 *
 * Updates the appropriate histogram (`wlHi`) and ln(g) (`wllngi`) for
 * the resulting bin. Out-of-range or masked-out targets are rejected
 * but still increment the source bin (paper §III-A boundary handling).
 */
int WangLandau(double Ei, double Ef, SamplingMode mode);
