/**
 * @file diagnostic.hpp
 * @brief Uniform error-reporting helpers.
 *
 * `die()` prints a formatted message to stderr and aborts the MPI job
 * with a non-zero exit code. It does not throw; callers that want
 * recoverable failure should `throw std::runtime_error` directly.
 */

#pragma once

#include <cstdarg>
#include <cstdio>
#include <cstdlib>

#include "mpi.h"

namespace deepthermo {

/**
 * @brief Internal implementation shared by ::die / ::die_with_code.
 *
 * @param rank_only Print on this rank only, or all ranks if negative.
 * @param exit_code Code passed to `MPI_Abort` / `std::exit`.
 * @param fmt printf-style format string.
 * @param ap Variadic argument list, already started by the caller.
 *
 * @note `[[noreturn]]` — never returns to the caller.
 */
[[noreturn]] inline void die_impl(int rank_only, int exit_code,
                                  const char* fmt, va_list ap) {
    int rank = 0;
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (initialized) MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank_only < 0 || rank == rank_only) {
        std::fprintf(stderr, "[deepthermo] fatal: ");
        std::vfprintf(stderr, fmt, ap);
        std::fprintf(stderr, "\n");
        std::fflush(stderr);
    }
    if (initialized) {
        MPI_Abort(MPI_COMM_WORLD, exit_code);
    } else {
        std::exit(exit_code);
    }
    std::abort();
}

/**
 * @brief Print a fatal error on rank 0 and abort the MPI job.
 *
 * @param fmt printf-style format string.
 *
 * @par Example
 * @code
 *   deepthermo::die("could not open coupling file: %s", path.c_str());
 * @endcode
 */
[[noreturn]] inline void die(const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    die_impl(0, 1, fmt, ap);
    va_end(ap);
}

/**
 * @brief Like ::die but with a caller-chosen exit code.
 *
 * @param exit_code Process exit / MPI_Abort code.
 * @param fmt printf-style format string.
 */
[[noreturn]] inline void die_with_code(int exit_code, const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    die_impl(0, exit_code, fmt, ap);
    va_end(ap);
}

}  // namespace deepthermo
