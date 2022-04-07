module load open-ce/1.4.0-py38-0
module load gcc/9.1.0
ROOT=/gpfs/alpine/scratch/junqi/stf011/smartsim
conda activate $ROOT/smartsim-env
export LD_LIBRARY_PATH=$ROOT/smartsim-env/lib/python3.8/site-packages/torch/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$ROOT/smartsim/smartsim/lib/backends/redisai_torch:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$ROOT/smartsim/smartsim/lib/backends/redisai_tensorflow:$LD_LIBRARY_PATH
export PATH=/gpfs/alpine/scratch/junqi/stf011/smartsim/smartsim/smartsim/bin:$PATH
