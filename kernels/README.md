## 3D convolution kernel benchmark

Modified from [DeepBench](https://github.com/baidu-research/DeepBench).
This is an independent microbenchmark — it does not link against the
DeepThermo simulation engine and uses its own Makefile under each
sub-directory.

### Requirements

- CUDA / cuDNN (Perlmutter, NVIDIA A100), or
- ROCm / MIOpen (Frontier, AMD MI250X)

### Quickstart

Perlmutter:
```bash
cd kernels/nvidia
sbatch conv3D.sb
```

Frontier:
```bash
cd kernels/amd
sbatch conv3D.sb
```

### Output

The launchers write `log.fp32` containing the per-shape timing table:

```
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   w      h      c      n      k      f_w      f_h    pad_w  pad_h    stride_w  stride_h    fwd_time (usec)  bwd_inputs_time (usec)  bwd_params_time (usec)  total_time (usec)   fwd_algo
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    2      2     64      1     64      3      3      1       1         2         2            76                      76                      76                228   ConvolutionFwdAlgoGEMM
```

`fwd_time / bwd_inputs_time / bwd_params_time` are in microseconds; the
shape columns describe one VAE conv layer (paper Table I) at the indicated
input geometry.
