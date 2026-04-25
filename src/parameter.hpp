/**
 * @file parameter.hpp
 * @brief TOML configuration loader.
 */

#pragma once

/**
 * @brief Parse a config.toml and populate the global ::deepthermo_sim fields.
 *
 * Accepted top-level tables:
 *   - `[lattice]` — N, NE, SH, elements, composition, coupling_file, ...
 *   - `[wang_landau]` — bin_width, e_min, e_max, flatness, mod_factor_*,
 *     iteration_factor, production_bin_samps
 *   - `[parallel_tempering]` — metropolis_sampling, T_init, T_final, dT,
 *     samples, sep, drop, restart
 *   - `[thermodynamics]` — T_init, T_final, dT
 *   - `[model]` — dir
 *
 * @param filename Path to the TOML file.
 * @throws std::runtime_error On parse failure, missing required field, or
 *         out-of-range value.
 */
void ReadInput(const char* filename);
