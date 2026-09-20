#include "doctest.h"

#include <stdexcept>
#include <string>

#include "parameter.hpp"
#include "sim_context.hpp"

namespace {
std::string fixture(const char* name) {
    return std::string{DEEPTHERMO_TESTS_FIXTURE_DIR} + "/" + name;
}
}

TEST_CASE("ReadInput populates SimContext from a valid TOML") {
    ReadInput(fixture("config_valid.toml").c_str());

    CHECK(alloyState.N == 4);
    CHECK(alloyState.NE == 4);
    CHECK(alloyState.SH == 6);
    CHECK(alloyState.element.size() == 5);  // [HEA, Mo, Nb, Ta, W]
    CHECK(alloyState.element[1] == "Mo");
    CHECK(alloyState.composition.size() == 4);
    CHECK(alloyState.composition[0] == doctest::Approx(0.25));

    CHECK(wlState.dWLD1 == doctest::Approx(0.005));
    CHECK(wlState.WLD1min == doctest::Approx(-1.30));
    CHECK(wlState.WLD1max == doctest::Approx(-1.20));
    CHECK(wlState.Flatness == doctest::Approx(0.6));
    CHECK(wlState.ProductionBinSamps == 5);

    CHECK(ptState.MetropolisSampling == 0);
    CHECK(ptState.MTi == doctest::Approx(10.0));
    CHECK(ptState.MTf == doctest::Approx(2000.0));

    CHECK(alloyState.TTi == doctest::Approx(100));
    CHECK(alloyState.TTf == doctest::Approx(3000));

    CHECK(alloyState.model_dir == "./models");
    CHECK(alloyState.coupling_file == "coupling.input");
}

TEST_CASE("ReadInput defaults composition to uniform when omitted") {
    ReadInput(fixture("config_no_composition.toml").c_str());
    CHECK(alloyState.composition.size() == static_cast<std::size_t>(alloyState.NE));
    for (auto v : alloyState.composition) {
        CHECK(v == doctest::Approx(1.0 / alloyState.NE));
    }
}

TEST_CASE("ReadInput throws on missing required field") {
    CHECK_THROWS_AS(ReadInput(fixture("config_missing_required.toml").c_str()),
                    std::runtime_error);
}

TEST_CASE("ReadInput throws on out-of-range value") {
    CHECK_THROWS_AS(ReadInput(fixture("config_bad_range.toml").c_str()),
                    std::runtime_error);
}

TEST_CASE("[output] snapshot capture defaults to off") {
    ReadInput(fixture("config_valid.toml").c_str());
    CHECK(alloyState.snapshot_stride == 0);
    CHECK(alloyState.snapshot_lowe == false);
}

TEST_CASE("[output] snapshot capture parses when present") {
    ReadInput(fixture("config_snapshots.toml").c_str());
    CHECK(alloyState.snapshot_stride == 10);
    CHECK(alloyState.snapshot_lowe == true);
}

TEST_CASE("ReadInput throws on negative snapshot_stride") {
    CHECK_THROWS_AS(ReadInput(fixture("config_bad_snapshot.toml").c_str()),
                    std::runtime_error);
}

TEST_CASE("ReadInput throws on malformed TOML") {
    CHECK_THROWS_AS(ReadInput(fixture("config_garbled.toml").c_str()),
                    std::runtime_error);
}
