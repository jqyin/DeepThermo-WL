/**
 * @file allreduce.h
 * @brief Ring-allreduce alternative to MPI_Allreduce.
 *
 * Bandwidth-optimal scatter-reduce + allgather over a ring topology;
 * Baidu's reference implementation. Enabled at build time with
 * `-DDEEPTHERMO_ENABLE_RING_ALLREDUCE` (CMake option).
 */

#ifndef BAIDU_ALLREDUCE_COLLECTIVES_H_
#define BAIDU_ALLREDUCE_COLLECTIVES_H_

#include <cstddef>
#include <mpi.h>

/**
 * @brief Ring allreduce on a buffer of `length` doubles.
 *
 * @param[in]  data   Input buffer, identical length on every rank.
 * @param[in]  length Number of doubles in `data`.
 * @param[out] output Pointer that will be set to a freshly allocated
 *                    buffer holding the reduced result.
 * @param[in]  rank   This rank's index.
 * @param[in]  size   Total number of ranks.
 */
void RingAllreduce(double* data, std::size_t length, double** output,
                   int rank, int size);

#endif  // BAIDU_ALLREDUCE_COLLECTIVES_H_
