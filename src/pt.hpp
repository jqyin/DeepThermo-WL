/**
 * @file pt.hpp
 * @brief Parallel tempering driver and Metropolis warm-up loop.
 */

#pragma once

#include "sim_context.hpp"

constexpr int CHPT_STEPS = 100; ///< Sweeps between checkpoint attempts.
constexpr int TIMER = 2;        ///< Wall-clock budget (hours) before checkpoint-then-exit.

/**
 * @brief Build the geometric temperature ladder for the PT ranks.
 *
 * @param Ti Lowest temperature.
 * @param Tf Highest temperature.
 * @param nT Number of replicas (= MPI world size).
 *
 * Populates `ptState.T[0..nT-1]` in geometric progression and assigns
 * `ptState.pT = T[myrank]`.
 *
 * @throws via ::deepthermo::die if `Ti * Tf` underflows.
 */
void ini_T(double Ti, double Tf, int nT);

/**
 * @brief Allocate `Atom` / `Atomo`, seed RNGs, and call ::initialize.
 */
void ini_sys();

/**
 * @brief Drive the parallel-tempering warm-up loop (paper §III-C stage 1).
 *
 * @param nT     Number of replicas.
 * @param DROPI  Equilibration sweeps to discard.
 * @param SAMPS  Number of post-equilibration samples.
 * @param SEP    Sweep count between recorded samples.
 * @param irun   Run-suffix for output filenames (`stat<irun>.dat`, ...).
 * @param mode   Acceptance regime (typically `metropolis` for the warm-up).
 *
 * Skips replica-exchange swaps automatically when `nprocs == 1`.
 */
void parallel_tempering(int nT, double DROPI, double SAMPS, double SEP,
                        int irun, SamplingMode mode);

/**
 * @brief One round of replica-exchange swap proposals.
 *
 * @param even If true, swap pairs (0↔1, 2↔3, ...); otherwise (1↔2, 3↔4, ...).
 *
 * Each rank either proposes a swap upward or accepts an incoming proposal
 * from below; ::ini_W is run after either path to refresh the pair tally.
 */
void swap(bool even);

/**
 * @brief Standard Metropolis acceptance test.
 *
 * @param Ei Initial total energy.
 * @param Ef Trial total energy.
 * @return 1 if accepted, 0 otherwise.
 */
int Metropolis(double Ei, double Ef);

/**
 * @brief Read a per-rank checkpoint state from `state<rank>.input`.
 * @throws via ::deepthermo::die if the file is missing or its energy
 *         disagrees with the recomputed energy by more than 1e-5.
 */
void read_state();

/**
 * @brief Write the per-rank state to `state<rank>.input` for restart.
 */
void write_state();

/**
 * @brief Run one full lattice's worth of bond-swap attempts.
 * @param mode Acceptance regime.
 */
void mchybrid(SamplingMode mode);
