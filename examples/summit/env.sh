module load gcc
export TF_DIR=/sw/summit/open-ce/anaconda-base/envs/open-ce-1.4.0-py39-0/lib/python3.9/site-packages/tensorflow
export REDIS_DIR=/gpfs/alpine/scratch/junqi/stf218/mcml-wl/verification/deepthermo-wl/envs/redis
export LD_LIBRARY_PATH=${TF_DIR}:${REDIS_DIR}/lib/backends/redisai_tensorflow:$LD_LIBRARY_PATH
export PATH=${REDIS_DIR}/bin:$PATH
