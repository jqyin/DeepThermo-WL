#!/bin/bash
module load rocm

rm -rf ${SMARTBENCH_LOGDIR}/build
mkdir -p ${SMARTBENCH_LOGDIR}/build && pushd ${SMARTBENCH_LOGDIR}/build

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda3-latest-Linux-x86_64.sh
BUILD=$(pwd)/miniconda
bash ./Miniconda3-latest-Linux-x86_64.sh -b -p $BUILD
export PATH=${BUILD}/bin:$PATH
conda create -y --prefix=$BUILD/tf_env python=3.8
source ${BUILD}/etc/profile.d/conda.sh
conda activate ${BUILD}/tf_env
conda update -n base -c defaults conda 
conda install -y pyyaml
pip install typing_extensions 
pip install numpy==1.19.5
pip install keras-preprocessing
export TEST_TMPDIR=/tmp
rm -rf /tmp/_bazel*
conda install -y bazel

git clone -b r2.5-rocm-enhanced https://github.com/ROCmSoftwarePlatform/tensorflow-upstream.git tensorflow
pushd tensorflow 
sed -i 's/\/usr\/bin\/python3/$(which python)/g' build_rocm_python3
sed -i 's/bazel build/bazel build --local_cpu_resources=32 --local_ram_resources=HOST_RAM*.1 /g' build_rocm_python3
export HCC_AMDGPU_TARGET=gfx90a 
export TF_ROCM_AMDGPU_TARGETS=gfx90a
module swap PrgEnv-cray PrgEnv-gnu
export CC=cc
export CXX=CC
module swap gcc gcc/10.3.0
module swap PrgEnv-gnu PrgEnv-cray

./build_rocm_python3 $ROCM_PATH &> log.tf
popd

git clone -b v0.21.3 --recursive https://github.com/horovod/horovod 
pushd horovod 
export HOROVOD_GPU=ROCM
export HOROVOD_WITH_TENSORFLOW=0 
export HOROVOD_WITH_PYTORCH=0 
export HOROVOD_WITHOUT_MXNET=1
export HOROVOD_ROCM_HOME=${ROCM_PATH}
export HOROVOD_GPU_ALLREDUCE=NCCL
export LIBPATHS=${ROCM_PATH}/lib:${ROCM_PATH}/hcc/lib:${ROCM_PATH}/lib64:${ROCM_PATH}/hsa/lib:${ROCM_PATH}/rocblas/lib
python setup.py install &> log.hvd
popd


tar -xf /opt/cray/pe/craype-dl-plugin-py3/21.04.1/wheel/dl_comm-21.4.1.tar.gz
pushd dl_comm-21.4.1
export CRAYPE_ML_PLUGIN_BASEDIR=/opt/cray/pe/craype-dl-plugin-py3/21.04.1
python setup.py install &> log.dl-plugin
popd

pip install transformers
python -c "import tensorflow as tf; print(tf.__version__)"
popd


