# DeepThermo: Deep Learning Accelerated Parallel Monte Carlo Sampling for Thermodynamics Evaluation of High Entropy Alloys 
This repository provides code for distributed Monte Carlo Sampling with Deep Learning generated proposals. The sampling methods include both Parallel Tempering and Wang-Landau sampling.  

## Software Requirements

- TensorFlow >= 2.4
- SmartRedis >= 0.2 
- Horovod >= 0.22
- gcc >= 9.0 
- MPI 

## Build 
Sampling code is under `src` with Makefiles provided for `Summit` and `Crusher`. 

3D convolution kernel benchmark code is under `kernels` for Nvidia `cuDNN` and AMD `MIOpen`.

Step-by-step quickstart can be followed for [MC sampling](https://code.ornl.gov/jqyin/deepthermo-wl/-/blob/vae-proposal/examples/README.md) and [kernel benchmarking](https://code.ornl.gov/jqyin/deepthermo-wl/-/blob/vae-proposal/kernels/README.md)

## VAE model training
VAE modeling is based on `tf.keras` using `Horovod` for distributed data parallel, with a reproducible Code Ocean [capsule](https://doi.org/10.24433/CO.8787331.v1).   

Three pretrained models for MoNbTaW of size 10<sup>3</sup>, 16<sup>3</sup>, 20<sup>3</sup> are provided under `models`. 

## Monte Carlo Sampling of High Entropy Alloys
Two example input and job scripts are provided for MoNbTaW on `Summit` and `Crusher` under `examples`, with [plot utility](https://code.ornl.gov/jqyin/deepthermo-wl/-/blob/vae-proposal/utils/plots.ipynb) avaliable.  




