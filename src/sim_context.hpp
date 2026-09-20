/**
 * @file sim_context.hpp
 * @brief Owned, RAII state container for a DeepThermo run.
 *
 * Replaces the four extern singletons (`mpiState` / `ptState` / `wlState` /
 * `alloyState`) that previously lived at file scope. Each sub-state is an
 * ordinary class with `std::vector`-backed storage; ctors allocate, dtors
 * free, no manual `initWL`/`freeWL`/`freePT` dance.
 *
 * All previously-`#define`'d alloy constants (NE, SH, MAX_NEIGHBORS, Z_R,
 * `reglin_intercept`, `E_scale`, `T_scale`, `gpus_per_node`) are now plain
 * fields, populated by parameter.cc from the TOML config file. Switching
 * alloy family no longer requires recompiling.
 *
 * The translation units still reference the four sub-states by their bare
 * names (`alloyState`, `wlState`, `ptState`, `mpiState`); these are
 * reference aliases to the members of a single global ::deepthermo_sim
 * instance, defined once in sim_context.cc.
 */

#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

/**
 * @namespace deepthermo
 * @brief Top-level namespace for the simulation engine.
 */
namespace deepthermo {

/**
 * @brief Sampling regime for a Wang–Landau / Metropolis step.
 *
 * Selects how `wanglandau::WangLandau` and `alloy::BondSwap` weight
 * acceptance and which density-of-states accumulator is updated. See
 * paper §III-C (DeepThermo sampling) for the regime semantics.
 */
enum SamplingMode {
    metropolis,    ///< Plain Metropolis acceptance (PT warm-up).
    WLprior,       ///< Prior-DOS estimation (paper Algo. 1 line 7).
    WLdos,         ///< Main WL iteration on the working DOS.
    WLproduction   ///< Production sweep — accumulates order parameters.
};

/**
 * @brief Lattice offset (i,j,k) for a neighbor site.
 *
 * `coupling.input` produces one of these per shell occupant; values are
 * applied with periodic wrap inside ::neighbor.
 */
struct Neighbor {
    int x = 0; ///< Offset along x.
    int y = 0; ///< Offset along y.
    int z = 0; ///< Offset along z.
};

// ---------------------------------------------------------------- MPI state

/**
 * @brief MPI rank metadata + per-rank GPU binding.
 *
 * Populated in `main` from `MPI_Comm_size`/`MPI_Comm_rank`;
 * `gpus_per_node` defaults to the platform-specific value
 * (4 for Perlmutter, 8 for Frontier) baked in at compile time and
 * is overridable from the config.
 */
class MPIState {
public:
    int nprocs = 0;        ///< World size from `MPI_Comm_size`.
    int myrank = 0;        ///< World rank from `MPI_Comm_rank`.
    int gpus_per_node = 1; ///< Per-rank GPU index = `myrank % gpus_per_node`.
};

// ---------------------------------------------------------- Parallel tempering

/**
 * @brief Parallel-tempering / Metropolis warm-up parameters and counters.
 *
 * Field names with the `M*` prefix come from the original tf.keras code
 * and refer to the standalone Metropolis "seq" loop used as a warm-up
 * before Wang–Landau (paper §III-C stage 1).
 */
class PTState {
public:
    int MetropolisSampling = 0; ///< 1 if running stand-alone Metropolis only.

    double MTi = 0.0;    ///< Metropolis seq() initial temperature.
    double MTf = 0.0;    ///< Metropolis seq() final temperature.
    double MdT = 0.0;    ///< Metropolis seq() temperature increment.
    double MSAMPS = 0.0; ///< Number of MC samples per `parallel_tempering` call.
    double MSEP = 0.0;   ///< Sweeps separating recorded samples.
    double MDROP = 0.0;  ///< Equilibration sweeps dropped at the start.

    int Restart = 0; ///< Whether this run resumes from `mc<rank>.input`.

    std::vector<double> T; ///< Geometric temperature ladder (size = nprocs).
    double pT = 0.0;       ///< Temperature assigned to this rank.

    int att = 0; ///< Attempted replica-exchange swaps on this rank.
    int acc = 0; ///< Accepted swaps on this rank.
};

// ---------------------------------------------------------------- Wang–Landau

/**
 * @brief Wang–Landau histograms, DOS, and iteration controls.
 *
 * `dWLD1` is the per-atom energy bin width; the engine internally bins
 * total energy with `invdWLD1 = 1 / (dWLD1 * invN)` (so `D1BINS = (e_max
 * - e_min) * invdWLD1`). The unsigned-short increment array (`wlHi`,
 * `wlHd`) is what gets allreduced — keeping the per-step message under
 * 2 KB so SHARP applies (paper §IV-B).
 */
class WLState {
public:
    double dWLD1 = 0.0;    ///< Bin width in per-atom energy units.
    double WLD1max = 0.0;  ///< Upper energy boundary.
    double WLD1min = 0.0;  ///< Lower energy boundary.
    int D1BINS = 0;        ///< Number of bins along the primary axis.
    double invdWLD1 = 0.0; ///< Inverse bin width (total energy).

    int numf = 1;             ///< Outer-iteration index (1 -> ModFactorInit).
    int TotalSweeps = 1;      ///< Cumulative sweep counter across all iters.
    int IterSweeps = 0;       ///< Sweep counter within the current iter.
    double Flatness = 0.0;    ///< Min/avg histogram threshold for "flat".
    double ModFactorInit = 0.0;  ///< Initial ln(f) (paper Algo. 1 line 1).
    double ModFactorFinal = 0.0; ///< Loop terminates once ln(f) drops below.
    double IterationFactor = 0.0; ///< Divides ln(f) per outer iter (typically 2).
    int ProductionBinSamps = 0;   ///< Sweeps in the post-WL production loop.
    double numbelow_flat = 0.0;   ///< Diagnostic — bins under the flat threshold.
    double lnwlf = 0.0;           ///< Current ln(f).

    std::vector<double> wlH;            ///< Accumulated histogram (across iters).
    std::vector<unsigned short> wlHd;   ///< Per-step rank-summed histogram delta.
    std::vector<unsigned short> wlHi;   ///< Per-step local histogram delta.
    std::vector<double> wllng;          ///< ln(g(E)).
    std::vector<double> wllng_prior;    ///< Prior DOS (paper §III-B).
    std::vector<double> wllngd;         ///< Reduced (allreduced) DOS scratch.
    std::vector<double> wllngi;         ///< Per-step local DOS delta.
    std::vector<int> mask;              ///< 0/1 mask: bins ineligible for the flatness check.
    std::vector<int> attemptrot;        ///< Per-bin attempt count (acceptance diagnostic).
    std::vector<int> acceptrot;         ///< Per-bin accept count.

    int LOWESTE = 0;                    ///< Cutoff bin index used by xyz dumping.
    std::vector<bool> print_list;       ///< One-shot guards for `write_xyz`.

    /**
     * @brief Allocate / reset all histogram and DOS arrays.
     * @param bins Number of bins along the primary energy direction.
     *
     * Sets `D1BINS = bins` and zero-initialises every array in this
     * object. `mask` is initialised to 1 across the board.
     */
    void resize(int bins);
};

// ------------------------------------------------------------------ Alloy

/**
 * @brief Lattice configuration, EPI couplings, and VAE-input scratch.
 *
 * The pair tables `J[i][j][s]` and `W[i][j][s]` are stored flat as
 * length `NE*NE*SH` `std::vector`s; access them through ::J_at and
 * ::W_at to avoid hand-rolled index arithmetic.
 */
class AlloyState {
public:
    /**
     * @brief Element table.
     *
     * Index 0 is the composite placeholder (`"HEA"`); indices 1..NE map
     * to actual species (e.g. `{"HEA", "Mo", "Nb", "Ta", "W"}`). The
     * `Atom` array stores 0..NE-1 so callers index `element[Atom[s]+1]`.
     */
    std::vector<std::string> element{"HEA", "Mo", "Nb", "Ta", "W"};

    int NE = 4;             ///< Number of chemical components.
    int SH = 6;             ///< Number of neighbor shells.
    int max_neighbors = 300;///< Hard cap for `nlist` allocation.

    /**
     * @brief Energy offset for the EPI model (paper Eq. 1, "+ const").
     *
     * Units: Ry. Different alloy families have different intercepts
     * (MoNbTaW: -1.27024..., MoNbTaVW: -1.74144..., MoNbTaTiW: -1.59459...).
     */
    double reglin_intercept = -1.2702430255548436;

    /// Latent-space random-walk radius (paper §III-B "δZ").
    double Z_R = 0.1;

    /// Energy-bin index below which `vae_update` fires (shrinks each WL iter).
    int E_0 = 0;

    int VAE_D = 0; ///< VAE input grid dimension (BCC-encoded, paper §III-B).
    int SHIFT = 0; ///< BCC coordinate shift (= N - 1).
    int PAD = 0;   ///< Padding to round VAE_D up to a multiple of 16.

    /// Per-call one-hot input buffer for the encoder. Size = VAE_D^3 * NE.
    std::vector<float> inputConfig;

    std::vector<double> op;  ///< Order parameter accumulator (size = D1BINS).
    std::vector<double> op2; ///< Order parameter squared accumulator.

    int N = 0;       ///< Lattice side length.
    int N_2 = 0;     ///< N^2.
    int N_3 = 0;     ///< N^3 (number of lattice sites).
    double invN = 0.0; ///< 1 / N^3.

    /// Current total energy (Ry units). Updated by ::Etot and incremental moves.
    double currEtot = 0.0;

    int attd = 0; ///< Attempted atom-swap moves on this rank.
    int accd = 0; ///< Accepted atom-swap moves on this rank.

    std::vector<short> Atom;   ///< Per-site atom type (size N^3).
    std::vector<short> Atomo;  ///< Pre-VAE checkpoint of `Atom` for revert.

    /// `elist[t]`: indices of sites currently occupied by element `t`.
    std::vector<std::vector<int>> elist;

    /// `inputPos[site]`: flattened neighbor site indices ordered by shell.
    std::vector<std::vector<int>> inputPos;

    /// Pair couplings J[i][j][s] flattened as NE x NE x SH (use ::J_at).
    std::vector<double> J;
    /// Pair-count tally W[i][j][s] flattened as NE x NE x SH (use ::W_at).
    std::vector<int> W;

    std::vector<int> NS;      ///< NS[s] = neighbors in shell s.
    std::vector<double> Dist; ///< Dist[s] = shell distance.
    std::vector<int> NT;      ///< NT[t] = count of element t.
    int Nneighbors = 0;       ///< Total neighbors across all shells.

    std::vector<Neighbor> nlist; ///< Per-neighbor offsets, size Nneighbors.

    double TTi = 0.0; ///< Thermodynamics output: initial temperature.
    double TTf = 0.0; ///< Thermodynamics output: final temperature.
    double dTT = 0.0; ///< Thermodynamics output: temperature step.

    int NBINTERACTION = 0; ///< Interaction selector (legacy).

    /**
     * @brief Per-element composition fractions used by ::ini_alloy.
     *
     * Sums to ~1.0. Replaces the legacy composition.input file.
     */
    std::vector<double> composition;

    std::string coupling_file = "coupling.input"; ///< Pair-coupling table path.
    std::string model_dir = "./models";           ///< VAE artefact directory.

    /**
     * @brief Capture training snapshots for the VAE (see ::write_xyz).
     *
     * `snapshot_stride > 0` appends every Nth parallel-tempering sample to
     * `snap_0_<rank>.xyz`; 0 (the default) disables PT capture entirely.
     * `snapshot_lowe` re-enables the one-shot dump taken the first time the
     * WL walker reaches each of the ten lowest energy slots, gated by
     * WLState::print_list, and lands in `snap_1_<rank>.xyz`.
     *
     * Both default to off so production runs pay no I/O cost; they exist to
     * generate the training set consumed by `vae-modeling`.
     */
    int snapshot_stride = 0;
    bool snapshot_lowe = false;

    /**
     * @brief Allocate the NE-dependent arrays (J, W, NS, Dist, NT).
     *
     * Call once after ::ReadInput finalises `NE` and `SH`. Safe to call
     * repeatedly; clears existing contents.
     */
    void resize_pair_arrays();

    /// @brief Mutable accessor into the flattened J table.
    inline double& J_at(int i, int j, int s)       { return J[(i * NE + j) * SH + s]; }
    /// @brief Read accessor into the flattened J table.
    inline double  J_at(int i, int j, int s) const { return J[(i * NE + j) * SH + s]; }
    /// @brief Mutable accessor into the flattened W table.
    inline int&    W_at(int i, int j, int s)       { return W[(i * NE + j) * SH + s]; }
    /// @brief Read accessor into the flattened W table.
    inline int     W_at(int i, int j, int s) const { return W[(i * NE + j) * SH + s]; }

    /// @brief Total entry count of J / W (= NE*NE*SH).
    inline std::size_t pair_count() const {
        return static_cast<std::size_t>(NE) * NE * SH;
    }
};

/**
 * @brief Physical unit conversions used by the energy model.
 *
 * The MC engine carries energies in Ry internally and converts to meV
 * at integration / printout points.
 */
struct Constants {
    double T_scale = 0.08618;       ///< kT in meV per Kelvin (= k_B * 1).
    double E_scale = 13605.69301;   ///< Rydberg -> meV conversion factor.
};

}  // namespace deepthermo

/**
 * @brief Forward-declared inference backend.
 *
 * Defined in `src/backend/inference_backend.hpp`. Forward-declared here
 * so SimContext can hold a `unique_ptr<InferenceBackend>` without every
 * translation unit pulling in LibTorch / TF / SmartRedis headers.
 */
class InferenceBackend;

namespace deepthermo {

/**
 * @brief Aggregate owner of all simulation state.
 *
 * Single instance ::deepthermo_sim is defined in sim_context.cc; the
 * four state references (::mpiState, ::ptState, ::wlState, ::alloyState)
 * alias its members for the existing `alloyState.foo` access pattern.
 *
 * Non-copyable, non-movable: the inference backend's GPU-attached
 * resources are not safe to duplicate.
 */
class SimContext {
public:
    SimContext();
    ~SimContext();
    SimContext(const SimContext&) = delete;
    SimContext& operator=(const SimContext&) = delete;

    Constants constants;           ///< Unit conversions.
    MPIState mpi;                  ///< MPI rank metadata.
    PTState pt;                    ///< Parallel-tempering state.
    WLState wl;                    ///< Wang–Landau histograms / DOS.
    AlloyState alloy;              ///< Lattice + couplings + VAE scratch.
    std::unique_ptr<InferenceBackend> backend; ///< Selected encoder/decoder.
};

}  // namespace deepthermo

/// @brief Single global SimContext owning the run's state.
extern deepthermo::SimContext deepthermo_sim;
/// @brief Reference alias for `deepthermo_sim.mpi`.
extern deepthermo::MPIState&   mpiState;
/// @brief Reference alias for `deepthermo_sim.pt`.
extern deepthermo::PTState&    ptState;
/// @brief Reference alias for `deepthermo_sim.wl`.
extern deepthermo::WLState&    wlState;
/// @brief Reference alias for `deepthermo_sim.alloy`.
extern deepthermo::AlloyState& alloyState;

// Bring the WL sampling-mode enum into the global namespace so the
// existing call sites can keep using bare names. New code may also write
// `deepthermo::SamplingMode` etc. directly.
using deepthermo::SamplingMode;
using deepthermo::metropolis;
using deepthermo::WLprior;
using deepthermo::WLdos;
using deepthermo::WLproduction;
