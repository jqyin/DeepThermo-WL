#include "doctest.h"

#include <set>

#include "rand.hpp"

TEST_CASE("KISS RNG: shelltimeseed makes randd1 reproducible") {
    shelltimeseed(42);
    const double a1 = randd1();
    const double a2 = randd1();
    const double a3 = randd1();

    shelltimeseed(42);
    CHECK(randd1() == doctest::Approx(a1));
    CHECK(randd1() == doctest::Approx(a2));
    CHECK(randd1() == doctest::Approx(a3));
}

TEST_CASE("KISS RNG: outputs in (0, 1)") {
    shelltimeseed(7);
    for (int i = 0; i < 10000; ++i) {
        const double r = randd1();
        CHECK(r > 0.0);
        CHECK(r < 1.0);
    }
}

TEST_CASE("KISS RNG: two different seeds give different first draws") {
    shelltimeseed(1);
    const double a = randd1();
    shelltimeseed(2);
    const double b = randd1();
    CHECK(a != b);
}
