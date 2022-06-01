SOURCE_DIR?=.
BIN_DIR?=bin
MKDIR=mkdir -p

#gpu target
#https://readthedocs.org/projects/jan20-test2/downloads/pdf/latest/
#If running on gfx906/vega20-type GPU(MI50, MI60)HCC_AMDGPU_TARGET=gfx906
AMDGPU=gfx90a


#hipcc
HIPCC=${ROCM_PATH}/bin/hipcc

#BLAS
ROCBLAS_LIB=rocblas
ROCBLAS_INC=${ROCM_PATH}/rocblas/include

#CONV
MIOPEN_LIB?=MIOpen
MIOPEN_INC?=${ROCM_PATH}/miopen/include

#DeepBench
DEEPBENCH_INC=${SOURCE_DIR}/../kernels

RCCL_INC_DIR=${ROCM_PATH}/rccl/include
RCCL_LIB_DIR=${ROCM_PATH}/rccl/lib
RCCL_LIB=rccl


all: conv 

#OPT=-g -O0 -fsanitize=undefined -fno-omit-frame-pointer
OPT=-O3

conv:
	$(MKDIR) $(BIN_DIR)
	$(HIPCC) ${SOURCE_DIR}/conv_bench_rocm.cpp -o $(BIN_DIR)/conv_bench -I$(DEEPBENCH_INC) -l$(MIOPEN_LIB) $(OPT) -std=c++11 --amdgpu-target=${AMDGPU}

rnn:
	$(MKDIR) $(BIN_DIR)
	$(HIPCC) ${SOURCE_DIR}/rnn_bench_rocm.cpp -o $(BIN_DIR)/rnn_bench -I$(DEEPBENCH_INC) -l$(MIOPEN_LIB) $(OPT) -std=c++11 --amdgpu-target=${AMDGPU}

gemm:
	$(MKDIR) $(BIN_DIR)
	$(HIPCC) ${SOURCE_DIR}/gemm_bench.cpp -o $(BIN_DIR)/gemm_bench -I$(DEEPBENCH_INC) -l$(ROCBLAS_LIB) $(OPT) -std=c++11 --amdgpu-target=${AMDGPU}

gemm_fp16:
	$(MKDIR) $(BIN_DIR)
	$(HIPCC) ${SOURCE_DIR}/gemm_bench_half.cpp -o $(BIN_DIR)/gemm_fp16_bench -I$(DEEPBENCH_INC) -l$(ROCBLAS_LIB) $(OPT) -std=c++11 --amdgpu-target=${AMDGPU}

rccl_single:
	$(MKDIR) $(BIN_DIR)
	$(HIPCC) ${SOURCE_DIR}/rccl_single_all_reduce.cpp -o $(BIN_DIR)/rccl_single_all_reduce -I$(DEEPBENCH_INC) -I$(RCCL_INC_DIR) -L$(RCCL_LIB_DIR) -l$(RCCL_LIB) $(OPT) -std=c++11 --amdgpu-target=${AMDGPU}

clean:
	rm -rf $(BIN_DIR)
                            
