// sim_context.hpp — owned, RAII state for a DeepThermo run.
//
// Replaces the four extern singletons (mpiState/ptState/wlState/alloyState)
// that previously lived at file scope. Each sub-state is an ordinary class
// with vector-backed storage; ctors allocate, dtors free, no manual
// initWL/freeWL/freePT dance.
//
// All previously-#define'd alloy constants (NE, SH, max_neighbors, Z_R,
// reglin_intercept, E_scale, T_scale, gpus_per_node) are now plain fields,
// populated by parameter.cc from the config file. Changing alloy family no
// longer requires recompiling.
//
// Allocations are vector-backed; ctors allocate, dtors free. There are no
// initWL/freeWL/freePT calls — just construct a SimContext and let it
// clean up at scope exit.

#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace deepthermo {

// WL modes, used across alloy / wanglandau / pt.
enum SamplingMode { metropolis, WLprior, WLdos, WLproduction };

struct Neighbor {
    int x = 0;
    int y = 0;
    int z = 0;
};

// ---------------------------------------------------------------- MPI state

class MPIState {
public:
    int nprocs = 0;
    int myrank = 0;

    // GPUs visible to each rank; platform default set by caller.
    int gpus_per_node = 1;
};

// ---------------------------------------------------------- Parallel tempering

class PTState {
public:
    int MetropolisSampling = 0;

    // Metropolis seq() loop parameters
    double MTi = 0.0, MTf = 0.0, MdT = 0.0;
    double MSAMPS = 0.0, MSEP = 0.0, MDROP = 0.0;

    int Restart = 0;

    std::vector<double> T;  // per-replica temperature ladder
    double pT = 0.0;        // this rank's temperature

    int att = 0;  // attempted swaps
    int acc = 0;  // accepted swaps
};

// ---------------------------------------------------------------- Wang–Landau

class WLState {
public:
    // Binning
    double dWLD1 = 0.0;
    double WLD1max = 0.0;
    double WLD1min = 0.0;
    int D1BINS = 0;
    double invdWLD1 = 0.0;

    // Iteration control
    int numf = 1;
    int TotalSweeps = 1;
    int IterSweeps = 0;
    double Flatness = 0.0;
    double ModFactorInit = 0.0;
    double ModFactorFinal = 0.0;
    double IterationFactor = 0.0;
    int ProductionBinSamps = 0;
    double numbelow_flat = 0.0;
    double lnwlf = 0.0;

    // Histogram / DOS arrays (sized in resize())
    std::vector<double> wlH;
    std::vector<unsigned short> wlHd;
    std::vector<unsigned short> wlHi;
    std::vector<double> wllng;
    std::vector<double> wllng_prior;
    std::vector<double> wllngd;
    std::vector<double> wllngi;
    std::vector<int> mask;
    std::vector<int> attemptrot;
    std::vector<int> acceptrot;

    int LOWESTE = 0;
    std::vector<bool> print_list;  // replaces fixed-size list[1000]

    void resize(int bins);
};

// ------------------------------------------------------------------ Alloy

class AlloyState {
public:
    // Element table. Index 0 is a composite placeholder ("HEA"); indices
    // 1..NE map to the actual species. Configured from the input file.
    std::vector<std::string> element{"HEA", "Mo", "Nb", "Ta", "W"};

    // Formerly compile-time constants.
    int NE = 4;             // number of chemical components
    int SH = 6;             // number of neighbor shells
    int max_neighbors = 300;

    // Energy offset for the EPI model (paper §II, "+ const"). Units: Ry.
    double reglin_intercept = -1.2702430255548436;

    // Latent-space random-walk radius (paper §III-B "δZ").
    double Z_R = 0.1;

    // Energy-range cutoff index below which vae_update fires.
    int E_0 = 0;

    // VAE tensor geometry; set by initialize().
    int VAE_D = 0;
    int SHIFT = 0;
    int PAD = 0;
    std::vector<float> inputConfig;

    // VAE order-parameter accumulators (sized to D1BINS).
    std::vector<double> op;
    std::vector<double> op2;

    // Lattice
    int N = 0;
    int N_2 = 0;
    int N_3 = 0;
    double invN = 0.0;

    // Current energy
    double currEtot = 0.0;

    // Atom-swap acceptance counters
    int attd = 0;
    int accd = 0;

    // Configurations and scratch
    std::vector<short> Atom;   // size N^3
    std::vector<short> Atomo;  // checkpoint copy

    // elist[t]: indices of sites currently occupied by element t.
    std::vector<std::vector<int>> elist;

    // inputPos[site]: flattened neighbor site indices ordered by shell.
    std::vector<std::vector<int>> inputPos;

    // Pair couplings J[i][j][s] and pair counts W[i][j][s] — flattened as
    // NE x NE x SH. Use J_at / W_at accessors.
    std::vector<double> J;
    std::vector<int> W;

    std::vector<int> NS;      // NS[s] = neighbors in shell s
    std::vector<double> Dist; // Dist[s] = shell distance
    std::vector<int> NT;      // NT[t] = count of element t
    int Nneighbors = 0;

    std::vector<Neighbor> nlist;  // size Nneighbors (<= max_neighbors)

    // Thermodynamics output grid
    double TTi = 0.0, TTf = 0.0, dTT = 0.0;

    int NBINTERACTION = 0;

    // Per-element initial composition fractions (sums ~1.0). Used by
    // ini_alloy() instead of the legacy composition.input file.
    std::vector<double> composition;

    // External files referenced by the run. Both default to the legacy
    // names; they are overridable via [model] / [lattice] in config.toml.
    std::string coupling_file = "coupling.input";
    std::string model_dir = "./models";

    // Allocates the NE-dependent arrays (J, W, NS, Dist, NT). Call once
    // NE and SH are finalised by ReadInput(). Safe to call repeatedly.
    void resize_pair_arrays();

    // Flat-index accessors for the NE x NE x SH tables.
    inline double& J_at(int i, int j, int s)       { return J[(i * NE + j) * SH + s]; }
    inline double  J_at(int i, int j, int s) const { return J[(i * NE + j) * SH + s]; }
    inline int&    W_at(int i, int j, int s)       { return W[(i * NE + j) * SH + s]; }
    inline int     W_at(int i, int j, int s) const { return W[(i * NE + j) * SH + s]; }
    inline std::size_t pair_count() const {
        return static_cast<std::size_t>(NE) * NE * SH;
    }
};

// Unit conversions.
struct Constants {
    double T_scale = 0.08618;       // kT  -> meV
    double E_scale = 13605.69301;   // Ry  -> meV
};

}  // namespace deepthermo

// Forward-declared so SimContext can hold one without pulling the
// heavyweight inference headers into every translation unit.
class InferenceBackend;

namespace deepthermo {

class SimContext {
public:
    SimContext();
    ~SimContext();
    SimContext(const SimContext&) = delete;
    SimContext& operator=(const SimContext&) = delete;

    Constants constants;
    MPIState mpi;
    PTState pt;
    WLState wl;
    AlloyState alloy;
    std::unique_ptr<InferenceBackend> backend;
};

}  // namespace deepthermo

// Legacy aliases. The existing translation units refer to the four states by
// bare names (alloyState, wlState, ptState, mpiState); these references keep
// those call sites working while the module-by-module rewrite lands. See
// src/sim_context.cc for the single definition of `deepthermo_sim` that owns
// the state, plus the aliases below.
extern deepthermo::SimContext deepthermo_sim;
extern deepthermo::MPIState&   mpiState;
extern deepthermo::PTState&    ptState;
extern deepthermo::WLState&    wlState;
extern deepthermo::AlloyState& alloyState;

// For the still-existing printf+exit error path in parameter.cc and similar.
// Will be replaced by a proper die() / exceptions in Task 7.
using deepthermo::SamplingMode;
using deepthermo::metropolis;
using deepthermo::WLprior;
using deepthermo::WLdos;
using deepthermo::WLproduction;
