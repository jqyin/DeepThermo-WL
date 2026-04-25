## DeepThermo example runs

Two ready-to-go configurations: `perlmutter/` (NVIDIA A100, NERSC) and
`frontier/` (AMD MI250X, OLCF). Both target MoNbTaW at N=10. Switch alloy
family by editing `config.toml` (`elements`, `composition`, `NE`,
`reglin_intercept`) — no recompile needed.

## Frontier

```bash
# 1. Build the binary (LibTorch ROCm)
export TORCH_INSTALL_PREFIX=/ccs/proj/.../libtorch-rocm
../../scripts/build-frontier.sh

# 2. Make sure the TorchScript model files are in place
mkdir -p models
# encoder_MoNbTaW.pt and decoder_MoNbTaW.pt are produced by
# vae-modeling/vae/src/export_torchscript.py once you have a trained
# vae.pt. See ../vae-modeling/README.md for the training pipeline.
cp /path/to/encoder_MoNbTaW.pt models/
cp /path/to/decoder_MoNbTaW.pt models/

# 3. Submit
cd examples/frontier
sbatch hea-wl.sb
```

## Perlmutter

```bash
export TORCH_INSTALL_PREFIX=/global/.../libtorch-cuda
../../scripts/build-perlmutter.sh

mkdir -p models
cp /path/to/encoder_MoNbTaW.pt models/
cp /path/to/decoder_MoNbTaW.pt models/

cd examples/perlmutter
sbatch hea-wl.sb
```

## Other backends

The default is LibTorch. If you have an existing TensorFlow SavedModel or
SmartRedis-served `.pb`, configure with the matching backend:

```bash
# TF C++ — uses the SavedModel directories under <model_dir>/encoder and
# <model_dir>/decoder.
cmake -B build -DDEEPTHERMO_BACKEND=tf -DTF_DIR=/path/to/tensorflow ...

# SmartRedis — co-locate a redis-server per GPU on each node before
# launching hea-wl; encoder_<alloy>.pb / decoder_<alloy>.pb live under
# <model_dir>/.
cmake -B build -DDEEPTHERMO_BACKEND=redis -DREDIS_DIR=/path/to/smartredis ...
```

The pre-shipped TF SavedModels for MoNbTaW at N=10/16/20 are under the
top-level `models/` directory.

## Outputs

After a successful run the working directory contains:

| file               | meaning                                                             |
|--------------------|---------------------------------------------------------------------|
| `run.dat`          | per-WL-iteration log: lnf, total sweeps, iter sweeps, flatness, etc |
| `comm.dat`         | allreduce timing for the WL inner loop                              |
| `DOS_H_iter*.dat`  | density of states + histogram per WL iteration                      |
| `therm.dat`        | thermodynamic quantities U, C, F, S, M, χ over the temperature grid |
| `compos.dat`       | atom-position dump after PT warm-up                                 |
| `infer.dat`        | VAE inference timing (only on the second WL iter, rank 0)           |

`utils/plots.ipynb` consumes these to reproduce the paper's Figs. 12 / 15.
