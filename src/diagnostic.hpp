// diagnostic.hpp — uniform error-reporting helpers.
//
// die(rank, fmt, ...) prints a formatted message to stderr (only on rank 0
// of MPI_COMM_WORLD if rank == 0, otherwise on every rank) and aborts the
// MPI job with a non-zero exit code. It does not throw; callers that want
// recoverable failure should throw std::runtime_error directly.

#pragma once

#include <cstdarg>
#include <cstdio>
#include <cstdlib>

#include "mpi.h"

namespace deepthermo {

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
    // Unreachable; satisfies [[noreturn]].
    std::abort();
}

[[noreturn]] inline void die(const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    die_impl(0, 1, fmt, ap);
    va_end(ap);
}

[[noreturn]] inline void die_with_code(int exit_code, const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    die_impl(0, exit_code, fmt, ap);
    va_end(ap);
}

}  // namespace deepthermo
