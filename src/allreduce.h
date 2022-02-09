#ifndef BAIDU_ALLREDUCE_COLLECTIVES_H_
#define BAIDU_ALLREDUCE_COLLECTIVES_H_ value

#include <cstddef>
#include <mpi.h>


// The ring allreduce. The lengths of the data chunks passed to this function
// must be the same across all MPI processes. The output memory will be
// allocated and written into `output`.
void RingAllreduce(double* data, size_t length, double** output);


#endif /* ifndef BAIDU_ALLREDUCE_COLLECTIVES_H_ */
