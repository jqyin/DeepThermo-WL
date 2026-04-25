#include "doctest.h"

#include <cstring>
#include <vector>

#include "alloy.hpp"
#include "sim_context.hpp"

// Etot reads alloyState directly. The fixture below seeds a tiny 2x2x2
// system with one nontrivial pair coupling and a known atom layout, then
// checks the formula produces the expected scalar.

namespace {

void seed_minimal_system() {
    auto& a = alloyState;
    a.NE = 2;
    a.SH = 1;
    a.N = 2;
    a.N_2 = 4;
    a.N_3 = 8;
    a.invN = 1.0 / 8.0;
    a.reglin_intercept = 0.0;
    a.resize_pair_arrays();

    // 4 unit-step nearest neighbors per site (we declare exactly 4 in this
    // shell). Equal coupling between species 0-1.
    a.NS[0] = 4;
    a.J_at(0, 1, 0) = -1.0;
    a.J_at(1, 0, 0) = -1.0;

    // All sites in species 0 — W_at(0,0,0) gets every neighbor pair, the
    // off-diagonal stays at 0, so Etot should be 0 (the J_at(0,0,0) term
    // is zero).
    a.Atom.assign(a.N_3, 0);
    // Manually populate W to avoid pulling in ini_W's neighbor-list
    // dependency (which needs inputPos to be populated by ini_apos).
    a.W_at(0, 0, 0) = a.NS[0] * a.N_3;
    a.W_at(0, 1, 0) = 0;
    a.W_at(1, 0, 0) = 0;
    a.W_at(1, 1, 0) = 0;
}

}  // namespace

TEST_CASE("Etot is zero for a single-species configuration") {
    seed_minimal_system();
    CHECK(Etot() == doctest::Approx(0.0));
}

TEST_CASE("Etot picks up reglin_intercept * N^3") {
    seed_minimal_system();
    alloyState.reglin_intercept = 0.125;  // Ry, contrived
    CHECK(Etot() == doctest::Approx(0.125 * 8));
}

TEST_CASE("Etot mixes pair-coupling correctly") {
    seed_minimal_system();
    // Manually flip half the off-diagonal weight on, scale by J = -1.
    alloyState.W_at(0, 0, 0) = 16;  // half the pairs are 0-0
    alloyState.W_at(0, 1, 0) = 16;  // half are 0-1
    // E += (W/NS/N_3) * J  for the i<j pair (i=0, j=1)
    //    = (16 / 4 / 8) * (-1.0)
    //    = -0.5
    // then * N_3 = -4
    CHECK(Etot() == doctest::Approx(-4.0));
}
