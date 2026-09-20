# DeepThermo: Deep Learning Accelerated Parallel Monte Carlo Sampling for High Entropy Alloys

DeepThermo is an MPI C++ engine for distributed Monte Carlo sampling of
high-entropy alloys, with large MC moves proposed by a pre-trained
variational autoencoder. It implements parallel tempering and Wang–Landau
sampling and scales to thousands of GPUs (paper: Yin, Wang, Shankar,
*DeepThermo*, IPDPS 2023, [DOI 10.1109/IPDPS54959.2023.00041](https://doi.org/10.1109/IPDPS54959.2023.00041)).

This is the **v2** layout (April 2026 refactor):

- Inference is pluggable. Three backends — **LibTorch (default)**,
  TensorFlow C++, SmartRedis/RedisAI — sit behind one
  `InferenceBackend` interface.
- Build system is CMake. Targets OLCF Frontier (AMD MI250X / ROCm) and
  NERSC Perlmutter (NVIDIA A100 / CUDA); both produce the same source.
- Configuration is a single `config.toml` (see `examples/`).
- VAE training and TorchScript export live in the
  [`vae-modeling`](https://code.ornl.gov/jqyin/deepthermo) submodule on
  branch `torch` (PyTorch + DDP, replacing the old tf.keras/Horovod stack).

## Build

Frontier:
```bash
git clone --recurse-submodules https://code.ornl.gov/jqyin/deepthermo-wl
cd deepthermo-wl
export TORCH_INSTALL_PREFIX=/path/to/libtorch-rocm
./scripts/build-frontier.sh
```

Perlmutter:
```bash
export TORCH_INSTALL_PREFIX=/path/to/libtorch-cuda
./scripts/build-perlmutter.sh
```

Generic (laptop / CI):
```bash
cmake -B build -DCMAKE_PREFIX_PATH=/path/to/libtorch-cpu -DDEEPTHERMO_BACKEND=torch
cmake --build build --parallel
```

Backend selection: `-DDEEPTHERMO_BACKEND={torch,tf,redis}`. The default is
`torch`. The TF and SmartRedis paths still work for users with existing
SavedModel / `.pb` artifacts; see `examples/README.md`.

## Run

```bash
cd examples/perlmutter        # or examples/frontier
sbatch hea-wl.sb              # produces run.dat, therm.dat, DOS_H_iter*.dat
```

A run needs three files in the working directory:
- `config.toml` — the simulation configuration (TOML).
- `coupling.input` — per-shell pair couplings (tabular).
- `models/encoder_<alloy>.pt`, `models/decoder_<alloy>.pt` —
  TorchScript artefacts produced by the submodule's `export_torchscript.py`.

## VAE training

**[`docs/vae-workflow.md`](docs/vae-workflow.md) walks the whole pipeline
end to end** — collecting configurations, preprocessing, training, export,
and measuring whether the trained global move actually helps — with the
pitfalls that bite in practice (the padded grid, element ordering,
unreachable WL windows, PT that does not mix).

In brief: PyTorch + `torch.distributed` DDP. See the submodule README and
`vae-modeling/vae/run-{frontier,perlmutter}.sb` for ready-to-run launchers.
Stage 1 (PT warm-up to collect training configs) is the same as the
production run; stage 2 (training) is in the submodule; stage 3 (WL
sampling with the trained VAE) runs the C++ binary built here.

## Tests

```bash
ctest --test-dir build --output-on-failure
```
Two suites: `deepthermo_unit` (doctest, ~0.5 s) and `deepthermo_smoke`
(end-to-end pipeline on a tiny config, ~3 s; soft-skips if `python3 +
torch + mpirun` aren't on PATH).

## Layout

```
src/                 simulation engine + InferenceBackend impls
src/backend/         {torch,tf,smartredis}_backend + factory
third_party/         vendored tomlplusplus, doctest (single-header)
cmake/               backend selection helpers
docs/                vae-workflow.md tutorial (+ generated Doxygen output)
scripts/             build wrappers per platform
examples/{frontier,perlmutter}/   ready-to-run configs + Slurm launchers
kernels/{nvidia,amd}              standalone 3D-conv kernel microbenchmarks
tests/{unit,smoke}                ctest suites
utils/               plots.ipynb post-processing notebook + freeze helper
vae-modeling/        submodule: PyTorch VAE training + TorchScript export
```

## Citing

```bibtex
@inproceedings{deepthermo2023,
  author    = {Junqi Yin and Feiyi Wang and Mallikarjun Shankar},
  title     = {DeepThermo: Deep Learning Accelerated Parallel Monte Carlo
               Sampling for Thermodynamics Evaluation of High Entropy Alloys},
  booktitle = {IPDPS},
  year      = {2023},
  doi       = {10.1109/IPDPS54959.2023.00041}
}
```
