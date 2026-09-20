// parameter.cc — TOML-driven configuration loader.
//
// Replaces the legacy line-by-line key/value parser. Sections:
//
//   [lattice]              — system size / chemistry / pair-coupling table
//   [wang_landau]          — WL binning, flatness, modification factor schedule
//   [parallel_tempering]   — PT warm-up / Metropolis-only sweep parameters
//   [thermodynamics]       — temperature grid for thermoqs() output
//   [model]                — VAE model directory
//   [output]               — optional training-snapshot capture
//
// All values land directly in the SimContext fields. Range checks are
// inline; failures throw std::runtime_error which main.cc reports + aborts.

#include <stdexcept>
#include <string>
#include <vector>

#define TOML_HEADER_ONLY 1
#include "toml.hpp"

#include "parameter.hpp"
#include "sim_context.hpp"

namespace {

[[noreturn]] void fail(const std::string& msg) {
    throw std::runtime_error("config: " + msg);
}

template <typename T>
T require(const toml::table& tbl, std::string_view section, std::string_view key) {
    auto node = tbl[section][key];
    auto v = node.value<T>();
    if (!v) {
        fail(std::string{section} + "." + std::string{key} + " missing or wrong type");
    }
    return *v;
}

template <typename T>
T optional(const toml::table& tbl, std::string_view section, std::string_view key, T fallback) {
    return tbl[section][key].value_or(fallback);
}

void load_lattice(const toml::table& tbl) {
    auto& a = alloyState;
    a.N             = static_cast<int>(require<int64_t>(tbl, "lattice", "N"));
    a.NE            = static_cast<int>(optional<int64_t>(tbl, "lattice", "NE", 4));
    a.SH            = static_cast<int>(optional<int64_t>(tbl, "lattice", "SH", 6));
    a.max_neighbors = static_cast<int>(optional<int64_t>(tbl, "lattice", "max_neighbors", 300));
    a.NBINTERACTION = static_cast<int>(optional<int64_t>(tbl, "lattice", "nb_interaction", 1));
    a.reglin_intercept =
        optional<double>(tbl, "lattice", "reglin_intercept", -1.2702430255548436);
    a.Z_R = optional<double>(tbl, "lattice", "Z_R", 0.1);

    if (a.N < 1) fail("lattice.N must be >= 1");
    if (a.NE < 1) fail("lattice.NE must be >= 1");
    if (a.SH < 1) fail("lattice.SH must be >= 1");
    if (a.NBINTERACTION < 1 || a.NBINTERACTION > 4) {
        fail("lattice.nb_interaction must be in [1,4]");
    }

    if (auto* arr = tbl["lattice"]["elements"].as_array()) {
        a.element.clear();
        a.element.push_back("HEA");
        for (auto& el : *arr) {
            auto s = el.value<std::string>();
            if (!s) fail("lattice.elements entries must be strings");
            a.element.push_back(*s);
        }
        if (static_cast<int>(a.element.size()) != a.NE + 1) {
            fail("lattice.elements length must equal NE");
        }
    }

    if (auto* arr = tbl["lattice"]["composition"].as_array()) {
        a.composition.clear();
        for (auto& v : *arr) {
            auto d = v.value<double>();
            if (!d) fail("lattice.composition entries must be numeric");
            a.composition.push_back(*d);
        }
    }
    if (a.composition.empty()) {
        a.composition.assign(a.NE, 1.0 / a.NE);
    }
    if (static_cast<int>(a.composition.size()) != a.NE) {
        fail("lattice.composition length must equal NE");
    }

    a.coupling_file =
        optional<std::string>(tbl, "lattice", "coupling_file", "coupling.input");
}

void load_wang_landau(const toml::table& tbl) {
    auto& w = wlState;
    w.dWLD1            = require<double>(tbl, "wang_landau", "bin_width");
    w.WLD1min          = require<double>(tbl, "wang_landau", "e_min");
    w.WLD1max          = require<double>(tbl, "wang_landau", "e_max");
    w.Flatness         = require<double>(tbl, "wang_landau", "flatness");
    w.ModFactorInit    = require<double>(tbl, "wang_landau", "mod_factor_init");
    w.IterationFactor  = require<double>(tbl, "wang_landau", "iteration_factor");
    w.ModFactorFinal   = require<double>(tbl, "wang_landau", "mod_factor_final");
    w.ProductionBinSamps =
        static_cast<int>(require<int64_t>(tbl, "wang_landau", "production_bin_samps"));

    if (w.dWLD1 <= 1e-6 || w.dWLD1 > 10.0) fail("wang_landau.bin_width out of range");
    if (w.Flatness <= 0.0 || w.Flatness >= 1.0) fail("wang_landau.flatness must be in (0,1)");
    if (w.ModFactorInit <= 0.0) fail("wang_landau.mod_factor_init must be > 0");
    if (w.IterationFactor <= 0.0) fail("wang_landau.iteration_factor must be > 0");
    if (w.ModFactorFinal <= 0.0) fail("wang_landau.mod_factor_final must be > 0");
    if (w.ProductionBinSamps < 1) fail("wang_landau.production_bin_samps must be >= 1");
}

void load_thermo(const toml::table& tbl) {
    auto& a = alloyState;
    a.TTi = require<double>(tbl, "thermodynamics", "T_init");
    a.TTf = require<double>(tbl, "thermodynamics", "T_final");
    a.dTT = require<double>(tbl, "thermodynamics", "dT");
    if (a.TTi < 0 || a.TTf < 0 || a.dTT <= 0) {
        fail("thermodynamics temperatures must be >= 0 and dT > 0");
    }
}

void load_pt(const toml::table& tbl) {
    auto& p = ptState;
    p.MetropolisSampling =
        optional<bool>(tbl, "parallel_tempering", "metropolis_sampling", false) ? 1 : 0;
    p.MTi    = require<double>(tbl, "parallel_tempering", "T_init");
    p.MTf    = require<double>(tbl, "parallel_tempering", "T_final");
    p.MdT    = optional<double>(tbl, "parallel_tempering", "dT", 0.1);
    p.MSAMPS = require<double>(tbl, "parallel_tempering", "samples");
    p.MSEP   = require<double>(tbl, "parallel_tempering", "sep");
    p.MDROP  = require<double>(tbl, "parallel_tempering", "drop");
    p.Restart = static_cast<int>(optional<int64_t>(tbl, "parallel_tempering", "restart", 0));

    if (p.MTi < 0 || p.MTf < 0 || p.MdT < 0) fail("parallel_tempering temperatures must be >= 0");
    if (p.MSAMPS < 0 || p.MSEP < 0 || p.MDROP < 0) {
        fail("parallel_tempering counts must be >= 0");
    }
}

void load_model(const toml::table& tbl) {
    alloyState.model_dir = optional<std::string>(tbl, "model", "dir", "./models");
}

void load_output(const toml::table& tbl) {
    auto& a = alloyState;
    a.snapshot_stride =
        static_cast<int>(optional<int64_t>(tbl, "output", "snapshot_stride", 0));
    a.snapshot_lowe = optional<bool>(tbl, "output", "snapshot_lowe", false);

    if (a.snapshot_stride < 0) fail("output.snapshot_stride must be >= 0");
}

}  // namespace

void ReadInput(const char* filename) {
    toml::table tbl;
    try {
        tbl = toml::parse_file(filename);
    } catch (const toml::parse_error& e) {
        fail(std::string{"could not parse "} + filename + ": " +
             std::string{e.description()});
    }

    load_lattice(tbl);
    load_wang_landau(tbl);
    load_pt(tbl);
    load_thermo(tbl);
    load_model(tbl);
    load_output(tbl);
}
