/**
 * @file alloy.hpp
 * @brief Lattice initialisation, energy evaluation, and VAE-driven moves.
 *
 * All free functions read and mutate the global ::alloyState and the
 * other state objects through their reference aliases. `Etot` is the
 * EPI-model energy from paper Eq. 1; `vae_update` is the global MC move
 * defined in paper §III-C / Algo. 1 lines 6–9.
 */

#pragma once

#include "sim_context.hpp"

/**
 * @brief Initialise lattice scratch arrays and VAE input geometry.
 *
 * Reads `coupling_file`, sizes pair tables, computes `VAE_D` / `SHIFT`
 * / `PAD`, and allocates `inputConfig`. Call after ::ReadInput.
 */
void initialize();

void ini_conf();

/**
 * @brief Build the per-site neighbor index list `inputPos`.
 *
 * Must be called after `Atom` is populated; also computes the initial
 * pair-count tally `W` via ::ini_W.
 */
void ini_apos();

/**
 * @brief Populate the lattice with atoms drawn from the configured composition.
 *
 * @param state Pass 0 for a fresh shuffled init; nonzero leaves the
 *              ordered placement intact (used for tests / debugging).
 */
void ini_alloy(int state);

/**
 * @brief Read the per-shell pair couplings from `coupling_file`.
 *
 * Allocates `nlist`, `J`, `W`, `NS`, `Dist`, `NT`, `elist`,
 * `inputPos` based on `NE` / `SH` / N. Sized for at most `max_neighbors`.
 */
void ini_coupling();

double Eatom(int);

/**
 * @brief Total EPI energy of the current configuration (paper Eq. 1).
 * @return Energy in Ry, including the `reglin_intercept` constant term.
 *
 * @note Reads `alloyState.W` and `alloyState.J`; assumes they have been
 *       maintained incrementally by ::BondSwap / ::vae_update or
 *       initialised by ::ini_W.
 */
double Etot();

/**
 * @brief Generate the per-shell offsets for one site, with periodic wrap.
 * @param i Site lattice index along x.
 * @param j Site lattice index along y.
 * @param k Site lattice index along z.
 * @param[out] nn Flat int array sized at least `Nneighbors * 3`.
 */
void neighbor(int i, int j, int k, int* nn);

double Esite(int i, int j, int k);

void wolff(int i, int j, int k, double rx, double ry, double rz);

/**
 * @brief Try one nearest-neighbor pair swap.
 *
 * Picks two neighboring sites, computes the energy delta, and accepts
 * or rejects according to `mode`'s acceptance rule (Metropolis or one
 * of the WL variants).
 *
 * @param mode Acceptance regime (paper §III-C).
 */
void BondSwap(SamplingMode mode);

/**
 * @brief Try the VAE-generated global update (paper Algo. 1 lines 6–9).
 *
 * Encodes the current configuration to latent space, walks by `Z_R`,
 * decodes back to a candidate configuration, applies a conservation
 * check, evaluates the new energy, and accepts or rejects.
 *
 * @param mode Acceptance regime.
 */
void vae_update(SamplingMode mode);

/**
 * @brief Run the encoder + decoder and write the result to `Atom`.
 *
 * @param[in,out] z Latent-space coordinate. Decoded from this value.
 *
 * @note Used as a building block by ::vae_update; rarely called directly.
 */
void decode(float* z);

/**
 * @brief Take one Gaussian-distributed step in the latent space.
 * @param[in,out] npos 3-vector of latent coordinates; updated in place.
 *
 * Step magnitude is `alloyState.Z_R`.
 */
void walk(float* npos);

void Vol();

/**
 * @brief L1 norm of the latent encoding of the current configuration.
 *
 * Used as a scalar order-parameter proxy for `OrderParameter` and
 * `parallel_tempering` statistics.
 *
 * @return Sum of |z_0| + |z_1| + |z_2|.
 */
double L1();

/**
 * @brief Dump the current configuration to compos.dat (rank 0 only).
 */
void write_pos();

/**
 * @brief Append the current configuration to a per-rank xyz snapshot.
 *
 * Emits one standard two-header-line xyz frame: the atom count, then a
 * comment line carrying `Eng`, `MC_step` and the short-range-order table,
 * then one `element x y z` line per site. Frames accumulate in the file, so
 * repeated calls with the same `group` build a trajectory that
 * `vae-modeling/preprocessing/create_vae_input.py` can consume directly.
 *
 * Driven by AlloyState::snapshot_stride / AlloyState::snapshot_lowe; both
 * are off by default.
 *
 * @param step  MC step recorded on the comment line.
 * @param group Trajectory id used in the filename (`snap_<group>_<rank>.xyz`);
 *              0 for PT samples, 1 for the one-shot low-energy dumps.
 */
void write_xyz(int step, int group);

/**
 * @brief Compute thermodynamic quantities from the converged DOS.
 *
 * Writes therm.dat with columns `T U C F S M chi` over the
 * configured temperature grid. Paper Eqs. 8–11.
 */
void thermoqs();

/**
 * @brief Recompute the pair-count tally W from the current `Atom`.
 */
void ini_W();

/**
 * @brief Update the order-parameter accumulator at the given energy bin.
 *
 * @param idx Bin index into `alloyState.op` / `alloyState.op2`.
 *
 * Called from the WL acceptance path under `WLproduction`.
 */
void OrderParameter(int idx);
