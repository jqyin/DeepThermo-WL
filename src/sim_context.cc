#include "sim_context.hpp"

namespace deepthermo {

void WLState::resize(int bins) {
    D1BINS = bins;
    wlH.assign(bins, 0.0);
    wlHd.assign(bins, 0);
    wlHi.assign(bins, 0);
    wllng.assign(bins, 0.0);
    wllng_prior.assign(bins, 0.0);
    wllngd.assign(bins, 0.0);
    wllngi.assign(bins, 0.0);
    mask.assign(bins, 1);
    attemptrot.assign(bins, 0);
    acceptrot.assign(bins, 0);
    print_list.assign(bins, false);
}

void AlloyState::resize_pair_arrays() {
    J.assign(pair_count(), 0.0);
    W.assign(pair_count(), 0);
    NS.assign(SH, 0);
    Dist.assign(SH, 0.0);
    NT.assign(NE, 0);
    nlist.clear();
    nlist.reserve(max_neighbors);
    elist.assign(NE, std::vector<int>{});
}

}  // namespace deepthermo

// The single owner. main() never touches state directly; everything flows
// through this object and the reference aliases below.
deepthermo::SimContext deepthermo_sim;

deepthermo::MPIState&   mpiState   = deepthermo_sim.mpi;
deepthermo::PTState&    ptState    = deepthermo_sim.pt;
deepthermo::WLState&    wlState    = deepthermo_sim.wl;
deepthermo::AlloyState& alloyState = deepthermo_sim.alloy;
