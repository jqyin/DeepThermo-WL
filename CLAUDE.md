# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

DeepThermo-WL is an MPI C++ engine for distributed Monte Carlo sampling (parallel tempering + Wang–Landau) of high-entropy alloys, with large MC moves proposed by a pre-trained VAE. The repo went through a v2 refactor (April 2026); the layout, build system, input format, and default backend all changed at once. If you are reading older commit messages they reference dropped concepts (the makefile system, the `wlrun1D.input` parser, the four extern-globals pattern, Summit/Crusher targets, the inline-defined `model.hpp`/`client.hpp`).

## Workflow tutorial

`docs/vae-workflow.md` is the end-to-end VAE pipeline as runnable commands
(collect → preprocess → split → train → export → WL run → verify), with the
measured numbers from a real pass and the known pitfalls. Start there for
anything touching the training loop or the global move.

## Paper reference

Yin, Wang, Shankar, *DeepThermo*, IPDPS 2023, DOI 10.1109/IPDPS54959.2023.00041. PDF on this machine at `~/Downloads/DeepThermo_Deep_Learning_Accelerated_Parallel_Monte_Carlo_Sampling_for_Thermodynamics_Evaluation_of_High_Entropy_Alloys.pdf`. Run `pdftotext -layout <pdf> /tmp/deepthermo_paper.txt` (poppler is at `/opt/homebrew/bin`) when you need the body. Algorithm 1 is the spec for `main.cc`'s outer WL loop; Table I is the encoder/decoder graph spec — preserve graph node names (`input_1`, `input_2`, `Identity`) when retraining.

## Build

CMake-only since v2:

```
# Frontier (ROCm + LibTorch-rocm)
TORCH_INSTALL_PREFIX=<libtorch> ./scripts/build-frontier.sh

# Perlmutter (CUDA + LibTorch-cuda)
TORCH_INSTALL_PREFIX=<libtorch> ./scripts/build-perlmutter.sh

# Generic / laptop / CI
cmake -B build -DCMAKE_PREFIX_PATH=<libtorch> -DDEEPTHERMO_BACKEND=torch
cmake --build build --parallel
```

Top-level `CMakeLists.txt` exposes:
- `DEEPTHERMO_BACKEND` ∈ `{torch (default), tf, redis}`
- `DEEPTHERMO_PLATFORM` ∈ `{generic (default), frontier, perlmutter}`
- `DEEPTHERMO_BUILD_TESTS` (default ON)
- `DEEPTHERMO_ENABLE_GLOBAL_UPDATE` (default ON; emits `-DGLOBAL_UPDATE` so the WL outer loop calls `global_update()`)
- `DEEPTHERMO_ENABLE_RING_ALLREDUCE` (default OFF; the ring path in `allreduce.cc` instead of `MPI_Allreduce`)

The CMake produces `deepthermo_core` (static lib, no `main.cc`) and a thin `hea-wl` executable. The static lib is what `tests/` links against.

## Run

```
mpirun -n <ranks> ./hea-wl config.toml <seed>
```

Working directory must contain `config.toml`, `coupling.input`, and a `models/` dir holding `encoder_<alloy>.pt` + `decoder_<alloy>.pt` (TorchScript). For the TF / SmartRedis backends, model files are TF SavedModel dirs / frozen `.pb` instead.

`examples/{frontier,perlmutter}/` are the reference run dirs. The legacy `examples/{summit,crusher}/` directories are gone — Summit is decommissioned and Crusher succeeded by Frontier.

## Architecture

### State

All simulation state lives in one `deepthermo::SimContext` instance, declared in `src/sim_context.hpp` and defined in `src/sim_context.cc`. SimContext aggregates four sub-classes — `MPIState`, `PTState`, `WLState`, `AlloyState` — plus a `unique_ptr<InferenceBackend>`. Allocations are `std::vector`-backed; ctor/dtor handle setup / teardown so there are no `freeWL`/`freePT` calls.

For backwards compatibility with existing access patterns the four sub-state objects are also exposed as bare references (`alloyState`, `wlState`, `ptState`, `mpiState`), defined once in `sim_context.cc`. New code should still reference them by these names — they are *the* canonical handles.

Constants previously hard-coded as `#define`s (`NE`, `SH`, `MAX_NEIGHBORS`, `Z_R`, `reglin_intercept`, `gpus_per_node`, the element table) are now runtime fields populated by `parameter.cc` from `config.toml`. Pair tables (`J[NE][NE][SH]`, `W[NE][NE][SH]`) are flat `std::vector<…>` with `J_at(i,j,s)` / `W_at(i,j,s)` accessors. Switching alloy family no longer requires recompiling.

### Backend abstraction

```
src/backend/
  inference_backend.hpp     pure-virtual interface (load / encode / decode / name)
  torch_backend.{hpp,cc}    LibTorch (default)
  tf_backend.{hpp,cc}       TensorFlow C++ SavedModel
  smartredis_backend.{hpp,cc}  SmartRedis / RedisAI
  backend_factory.cc        compile-time-selected make_backend()
```

`alloy.cc::encode/decode` build the one-hot input buffer in NHWDC layout and hand it to `deepthermo_sim.backend->encode/decode`. Each backend internally adapts to its native layout — for LibTorch the exported wrapper modules permute NHWDC↔NCDHW, so the C++ side never deals with channels-first.

### Pipeline

`main.cc` follows paper §III-C three-stage pipeline:
1. **Stage 1 — PT warm-up.** `parallel_tempering(..., metropolis)` runs Metropolis sweeps and replica-exchange swaps to bracket the energy window. Stage 2 (VAE training) is in the `vae-modeling` submodule, branch `torch`; the trained TorchScript artefacts feed back into stage 3.
2. **Stage 3 — WL outer loop.** Mod factor `lnwlf` halves each iter (paper Algo. 1 line 1). Inside each iter `sweepWL(WLdos)` runs until histogram flatness ≥ `Flatness`. The θ gate at `main.cc` (`TotalSweeps % pow(2, numf) == 0 && currEtot < (E_0/invdWLD1 + WLD1min)*N_3`) decides when `vae_update(WLdos)` fires. Each sweep all-reduces `wlHi → wlHd` and folds into the global `wlH`/`wllng`. `E_0` shrinks every iter (paper "lower 20% of energy range" cutoff narrows).
3. **Production run.** Fixed `ProductionBinSamps` sweeps with `WLproduction` mode; accumulates `alloyState.op` / `op2` for `thermoqs()` to write `therm.dat`.

### Module map

- `alloy.cc` — lattice, neighbor shells from `coupling.input`, `Etot`/`Esite`, local moves (`BondSwap`, `wolff`), VAE-driven global move (`vae_update` → `encode` → walk in latent space → `decode` → conservation check). `NE`/`SH` taken from `alloyState` per call.
- `wanglandau.cc` — WL bins, DOS/histogram update, flatness test, mask I/O (`g.dat`/`mask.dat`).
- `pt.cc` — temperature ladder, replica-exchange step (`swap`), Metropolis warm-up. Uses `std::mt19937` (replaced the old header-defined `random.h` Mersenne Twister in v2).
- `parameter.cc` — TOML loader; populates SimContext fields. Throws `std::runtime_error` on parse / validation failure; `main.cc` catches and `MPI_Abort`s.
- `diagnostic.hpp` — `deepthermo::die(fmt, …)` helper used at fatal error sites across the engine.
- `rand.{cc,hpp}` — KISS RNG used by the WL acceptance tests and lattice shuffle.
- `allreduce.cc` — Baidu-style ring allreduce; enabled with `-DDEEPTHERMO_ENABLE_RING_ALLREDUCE`.

### Variant builds

The pre-v2 `-DTime_Series` / `-DHIST` build variants were removed in the cleanup pass — they had bit-rotted (referenced bare `currEtot`, indexed scalar `M` as `M[NE]`) and there is no record of any user actively building them. If per-step trajectories are needed, add them as a fresh feature on top of the new SimContext, not as a resurrection.

## Input format (TOML)

```
[lattice]              N, NE, SH, elements, composition, coupling_file,
                       max_neighbors, nb_interaction, reglin_intercept, Z_R
[wang_landau]          bin_width, e_min, e_max, flatness, mod_factor_init,
                       iteration_factor, mod_factor_final, production_bin_samps
[parallel_tempering]   metropolis_sampling (bool), T_init, T_final, dT,
                       samples, sep, drop, restart
[thermodynamics]       T_init, T_final, dT
[model]                dir
[output]               snapshot_stride, snapshot_lowe
```

Only `[lattice].N`, the WL set, the PT temperatures + counts, and the thermo grid are required. Everything else has a default (uniform composition, MoNbTaW intercept, etc.).

`[output]` drives VAE training-set capture via `alloy.cc::write_xyz`; both knobs default to off so production runs pay no I/O cost. `snapshot_stride = n` appends every nth PT sample to `snap_0_<rank>.xyz`; `snapshot_lowe = true` captures the first visit to each of the ten lowest WL energy slots (gated by `wlState.print_list`) into `snap_1_<rank>.xyz`. Frames are standard two-header-line xyz with `Eng`/`MC_step`/`SRO` on the comment line, consumed by `vae-modeling/preprocessing/create_vae_input.py`.

## Tests

```
ctest --test-dir build --output-on-failure
```

- `deepthermo_unit` — doctest, ~0.5 s. Covers TOML parser, SimContext sizing, `Etot` on a hand-seeded 2³ system, KISS RNG determinism. No MPI runtime needed.
- `deepthermo_smoke` — `tests/smoke/run.sh`. Exports a fresh random-weight TorchScript VAE via the submodule's `export_torchscript.py`, drops a tiny TOML config + `coupling.input` into a tmpdir, runs `mpirun -n 1 ./hea-wl`, and verifies `run.dat`/`therm.dat`/`comm.dat`/`DOS_H_iter*.dat` materialise. Soft-skips if `python3 + torch + mpirun` aren't on PATH. The smoke config short-circuits the WL outer loop (`mod_factor_init < mod_factor_final`) — it's a pipeline regression check, not WL-physics validation.

## Paper ↔ code map (v2 names)

- **Algo. 1 line 5, θ gate** → `main.cc` (`TotalSweeps % pow(2, wlState.numf) == 0 && alloyState.currEtot < (alloyState.E_0/wlState.invdWLD1 + wlState.WLD1min) * alloyState.N_3`).
- **Algo. 1 line 7, δZ random walk** → `alloy.cc::walk` with step radius `alloyState.Z_R` (default 0.1, configurable via `[lattice].Z_R`).
- **Algo. 1 line 9, conservation check** → tail of `alloy.cc::decode` after argmax.
- **Eq. 1 EPI energy** → `alloyState.J_at(i,j,s)` table loaded from `coupling.input` + `alloyState.reglin_intercept`; computed in `alloy.cc::Etot`.
- **Eq. 5 WL transition** → `wanglandau.cc::WangLandau`.
- **Unsigned-short H(E) allreduce** → `wlState.wlHi`/`wlHd`; `wllng` rebuilt locally after reduce. Paper §IV-B explains the message-size budget that puts H(E) under 2 KB so SHARP applies.
- **GPUs per node** → `mpiState.gpus_per_node`, defaulted from the `DEEPTHERMO_PLATFORM_*` macro at MPI init (4 for Perlmutter, 8 for Frontier).

## Submodule

`vae-modeling` is on branch `torch` (PyTorch + DDP). The `tfv2` branch retains the legacy tf.keras + Horovod stack for users still on the TF backend. The C++ TorchBackend loads files named `encoder_<alloy_tag>.pt` / `decoder_<alloy_tag>.pt` (e.g. `encoder_MoNbTaW.pt`); these come out of `vae-modeling/vae/src/export_torchscript.py`. Wrapper modules in that file internally permute NHWDC↔NCDHW so the C++ side never knows about layout.
