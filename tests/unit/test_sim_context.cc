#include "doctest.h"

#include "sim_context.hpp"

TEST_CASE("AlloyState resize_pair_arrays sizes flat tables") {
    deepthermo::AlloyState a;
    a.NE = 5;
    a.SH = 6;
    a.max_neighbors = 64;
    a.resize_pair_arrays();

    CHECK(a.J.size() == static_cast<std::size_t>(a.NE * a.NE * a.SH));
    CHECK(a.W.size() == static_cast<std::size_t>(a.NE * a.NE * a.SH));
    CHECK(a.NS.size() == static_cast<std::size_t>(a.SH));
    CHECK(a.Dist.size() == static_cast<std::size_t>(a.SH));
    CHECK(a.NT.size() == static_cast<std::size_t>(a.NE));
    CHECK(a.elist.size() == static_cast<std::size_t>(a.NE));
}

TEST_CASE("AlloyState J_at / W_at indexers round-trip") {
    deepthermo::AlloyState a;
    a.NE = 3;
    a.SH = 2;
    a.resize_pair_arrays();

    a.J_at(0, 1, 1) = 7.5;
    a.J_at(2, 0, 0) = -3.0;
    a.W_at(1, 1, 1) = 42;

    CHECK(a.J_at(0, 1, 1) == doctest::Approx(7.5));
    CHECK(a.J_at(2, 0, 0) == doctest::Approx(-3.0));
    CHECK(a.W_at(1, 1, 1) == 42);
    CHECK(a.J_at(0, 0, 0) == doctest::Approx(0.0));  // untouched
}

TEST_CASE("WLState resize sizes histogram + DOS arrays") {
    deepthermo::WLState w;
    w.resize(123);
    CHECK(w.D1BINS == 123);
    CHECK(w.wlH.size() == 123);
    CHECK(w.wlHd.size() == 123);
    CHECK(w.wllng.size() == 123);
    CHECK(w.mask.size() == 123);
    for (auto m : w.mask) CHECK(m == 1);
    for (auto h : w.wlH) CHECK(h == 0.0);
    CHECK(w.print_list.size() == 123);
}

TEST_CASE("SimContext is non-copyable but constructible") {
    deepthermo::SimContext ctx;
    ctx.alloy.NE = 4;
    ctx.alloy.SH = 6;
    ctx.alloy.resize_pair_arrays();
    CHECK(ctx.alloy.J.size() == 96u);
    // explicit non-copy / non-move guarantee
    static_assert(!std::is_copy_constructible_v<deepthermo::SimContext>);
    static_assert(!std::is_copy_assignable_v<deepthermo::SimContext>);
}
