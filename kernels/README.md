## 3D convolution kernel benchmark

Modified from [DeepBench](https://github.com/baidu-research/DeepBench)

### Software Requirements:

- cuda >= 10.1
- rocm >= 4.5 

### Quickstart

- Clone the repo
```bash
git clone https://code.ornl.gov/jqyin/deepthermo-wl
cd kernels
```

- Submit the job script that will complie the source code and run the benchmark

Summit 
```bash
cd nvidia
bsub conv3D.lsf
```
Crusher
```bash 
cd amd
sbatch conv3D.sb
```

- The generated log file contains the run time (usec) for various input parameters

```bash
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   w      h      c      n      k      f_w      f_h    pad_w  pad_h    stride_w  stride_h    fwd_time (usec)  bwd_inputs_time (usec)  bwd_params_time (usec)  total_time (usec)   fwd_algo
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    2      2     64      1     64      3      3      1       1         2         2            76                      76                      76                228   ConvolutionFwdAlgoGEMM
```
The above shows the batch (n) 1 input (w, h, c) of (2x2x64) with filter size (3), padding (1), and stride (2) takes 76 usec for forward pass with ConvolutionFwdAlgoGEMM algorithm     

