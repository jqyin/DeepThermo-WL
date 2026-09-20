# End-to-end VAE workflow

How to go from a freshly built engine to a Wang–Landau run driven by a trained
VAE global move, and how to tell whether the surrogate is actually helping.

This is the paper's three-stage pipeline (§III-C) written out as commands:
PT warm-up to collect configurations, VAE training in the `vae-modeling`
submodule, then WL sampling with the trained model.

Every number quoted below came from one real pass on a laptop (8 CPU cores,
MoNbTaW at N=10). They are there to tell you what "working" looks like, not as
targets to reproduce exactly.

**Contents**

1. [Prerequisites](#0-prerequisites)
2. [Collect training configurations](#1-collect-training-configurations)
3. [Preprocess to one-hot](#2-preprocess-to-one-hot)
4. [Split train/val](#3-split-trainval)
5. [Train](#4-train)
6. [Export to TorchScript](#5-export-to-torchscript)
7. [Run WL with the global move](#6-run-wl-with-the-global-move)
8. [Did it work?](#7-did-it-work)
9. [Pitfalls](#pitfalls)

---

## 0. Prerequisites

Build the engine and confirm the test suite is green:

```bash
cmake -B build -DCMAKE_PREFIX_PATH="$(python3 -c 'import torch; print(torch.utils.cmake_prefix_path)')" -DDEEPTHERMO_BACKEND=torch
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

`DEEPTHERMO_ENABLE_GLOBAL_UPDATE` must be ON (it is by default) or the WL loop
never calls `vae_update()` and step 6 measures nothing.

Make sure the submodule is checked out — the training code lives there:

```bash
git submodule update --init --recursive
```

### The bootstrapping problem

You need a loadable model *before* you have a trained one. `parallel_tempering`
computes its order parameter with `L1()` → `encode()`, so the engine loads the
encoder during PT even when you are only collecting data.

Export a random-weight model to get started. The configurations you collect
come from Metropolis/PT physics and are unaffected by the random weights —
only the order-parameter column is meaningless, and you do not use it here.

```bash
DT=/path/to/DeepThermo-WL          # repo root

mkdir -p run/models && cd run
PYTHONPATH="$DT/vae-modeling/vae" python3 - models <<'PY'
import sys, torch
from pathlib import Path
from src.export_torchscript import export
from src.model import VAE, VAEConfig
out = Path(sys.argv[1])
GRID, NE = 32, 4                      # GRID must be VAE_D for your lattice; see Pitfalls
cfg = VAEConfig(grid_size=GRID, num_elements=NE, latent_dim=3)
w = out / "bootstrap.pt"; torch.save(VAE(cfg).state_dict(), w)
export(w, "MoNbTaW", grid_size=GRID, num_elements=NE, latent_dim=3, out_dir=out)
PY
```

---

## 1. Collect training configurations

Snapshot capture is off by default so production runs pay no I/O cost. Turn it
on with the `[output]` section:

```toml
[output]
snapshot_stride = 1      # append every nth PT sample to snap_0_<rank>.xyz
snapshot_lowe   = false  # first visit to each of the 10 lowest WL slots -> snap_1_<rank>.xyz
```

For a collection-only run, short-circuit the WL outer loop
(`mod_factor_init < mod_factor_final`) and set `production_bin_samps = 1`, so
all the work is the PT warm-up. Keep the WL energy window wide and coarse —
its only job is to let `initWL()`'s relax loop exit immediately.

```toml
[parallel_tempering]
metropolis_sampling = false
T_init  = 10.0
T_final = 2000.0
samples = 400            # frames per rank
sep     = 10             # sweeps between samples (decorrelation)
drop    = 100            # equilibration sweeps

[wang_landau]
bin_width            = 0.5
e_min                = -1.35
e_max                = -1.00
mod_factor_init      = 0.001   # < mod_factor_final: skips the WL loop entirely
mod_factor_final     = 0.5
production_bin_samps = 1
```

Then run it. Each rank sits at its own ladder temperature, so the union across
ranks spans the energy range the surrogate has to model — use as many ranks as
you can.

```bash
mpirun -n 8 ../build/hea-wl config.toml 17
```

A full reference config is in
[`examples/frontier/config.toml`](../examples/frontier/config.toml).

**Check before moving on.** Energies should spread monotonically across the
ladder, and the span should bracket the WL window you intend to sample later:

```bash
for f in snap_0_*.xyz; do echo -n "$f: "; sed -n '2p' $f | grep -o "Eng = [-0-9.]*"; done
```

```
snap_0_0.xyz: Eng = -1281.0321     # T=10 K, ordered
...
snap_0_7.xyz: Eng = -1277.6024     # T=2000 K, disordered
```

---

## 2. Preprocess to one-hot

`--elements` must list species in the **same order** as `[lattice].elements` in
`config.toml`: the one-hot channel is the index into that table, so a different
order silently trains the model on permuted species.

```bash
cd vae-modeling
python3 preprocessing/create_vae_input.py \
    --xyz_file /path/to/run/snap_0_*.xyz \
    --elements Mo Nb Ta W \
    --skip_frames 50 \
    --save_file train_n10
```

```
lattice N=10 natoms=1000 grid=32^3 offset=11 elements=['Mo', 'Nb', 'Ta', 'W']
frames: 2800 (from 8 file(s), skipping 50 each)
wrote train_n10.npy shape=(2800, 32, 32, 32, 4) dtype=int8
train/export with --grid_size 32 --num_elements 4
```

Note the last line — carry that `--grid_size` through to training and export.
The grid comes from the engine's padding rule (`utils/geometry.py`), not from
the coordinate range in the file. See [Pitfalls](#the-grid-is-padded).

---

## 3. Split train/val

`train_vae.py` takes `--train_file` and `--val_file` separately. Do not pass
the same file twice: `VAEDataset` converts the whole array to float32 up front,
so you would pay the memory twice and validate on training data.

**Deduplicate before splitting.** Cold replicas repeat configurations (frozen
dynamics), and a repeat straddling the split would let validation score frames
the model memorised. Stratify by rank so validation spans the whole ladder
rather than landing entirely in the hottest replica.

```bash
python3 - <<'PY'
import numpy as np, hashlib
rng = np.random.default_rng(6)
a = np.load("train_n10.npy", mmap_mode="r")
n = a.shape[0]; per = n // 8              # 8 = number of PT ranks

seen = {}
for i in range(n):
    seen.setdefault(hashlib.md5(np.asarray(a[i]).tobytes()).digest(), i)
keep = np.array(sorted(seen.values()))
print(f"unique: {len(keep)}/{n}")

tr, va = [], []
for r in range(8):
    blk = rng.permutation(keep[(keep >= r*per) & (keep < (r+1)*per)])
    cut = max(1, int(round(0.1 * len(blk))))
    va.append(blk[:cut]); tr.append(blk[cut:])
tr = np.sort(np.concatenate(tr)); va = np.sort(np.concatenate(va))
np.save("train_split.npy", np.asarray(a[tr]))
np.save("val_split.npy",   np.asarray(a[va]))
print(f"train {len(tr)}  val {len(va)}")
PY
```

The unique count is itself a diagnostic. Ours was 2302/2800, and the loss was
concentrated at the cold end — see
[PT may not mix](#parallel-tempering-may-not-mix).

---

## 4. Train

```bash
python3 vae/train_vae.py \
    --train_file train_split.npy --val_file val_split.npy \
    --epochs 100 --batch_size 32 --latent_dim 3 --lr 5e-4 --seed 6 \
    --out_dir checkpoints
```

Multi-GPU is the same entry point under `torchrun`; see the submodule README.

```
epoch   1/100  train=3338835.8074  val=5029.6387
epoch  10/100  train=   1892.3     val=1849.1
epoch  50/100  train=   1130.0     val=1134.9
epoch 100/100  train=   1076.2     val=1088.9
```

That epoch-1 train figure is an initialisation transient, not a bug: the loss
is sum-reduced BCE over `32³ × 4 = 131072` voxels plus an early KL spike, and
the epoch average includes it. It recovers within one epoch. Watch the
train/val *gap* instead — ours closed to 1.2%, i.e. no overfitting.

100 epochs took ~30 min on 8 CPU cores (~19 s/epoch) for 2072 samples.

---

## 5. Export to TorchScript

```bash
python3 -m vae.src.export_torchscript \
    --weights checkpoints/vae.pt \
    --alloy_tag MoNbTaW \
    --grid_size 32 --num_elements 4 --latent_dim 3 \
    --out_dir exported
```

`--grid_size` is the value preprocessing printed in step 2. The wrappers
permute NHWDC↔NCDHW internally, so the C++ side never deals with layout.
Verify the round trip before handing it to the engine:

```bash
python3 -c "
import torch
e = torch.jit.load('exported/encoder_MoNbTaW.pt')
d = torch.jit.load('exported/decoder_MoNbTaW.pt')
x = torch.zeros(1,32,32,32,4)          # NHWDC, as the C++ hands it over
z = e(x); y = d(z)
print(tuple(z.shape), tuple(y.shape))  # (1, 3) (1, 32, 32, 32, 4)
"
```

Copy `encoder_<alloy>.pt` / `decoder_<alloy>.pt` into the run directory's
`models/`, replacing the bootstrap pair from step 0.

---

## 6. Run WL with the global move

Now run the real thing: WL outer loop enabled (`mod_factor_init >
mod_factor_final`) over the production energy window.

```toml
[wang_landau]
bin_width            = 0.0055
e_min                = -1.2808
e_max                = -1.2770
flatness             = 0.6
mod_factor_init      = 1.0
iteration_factor     = 2.0
mod_factor_final     = 1.0e-6
production_bin_samps = 10
```

```bash
mpirun -n 8 ../build/hea-wl config.toml 42
```

---

## 7. Did it work?

`vae.dat` records the VAE global-move acceptance: sweeps, `lnwlf`, cumulative
attempts, cumulative accepts, ratio. Cumulative, so difference successive rows
for the rate within an interval.

```bash
awk '{print $1, $2, $5}' vae.dat | tail -5
```

This is the number that matters. It is tracked separately from `attd`/`accd`
(which count `BondSwap`) because BondSwap fires N³ times per sweep and
`vae_update` at most once — sharing counters buries the signal under a 1000:1
dilution.

Our result, identical config and seed, varying only the model:

| | attempts | accepts | acceptance |
|---|---|---|---|
| trained VAE | 29146 | 17275 | **59.3%** |
| random weights | 35992 | 326 | 0.9% |

Acceptance decays from 63% to 41% as `lnwlf` goes 1 → 0.016. That is expected:
the DOS flattens and the WL criterion tightens.

**High acceptance alone is not success.** If the decoder reconstructs its input
almost exactly, the "global" move barely changes anything and is accepted
trivially. Check the move size too:

| | sites changed per proposed move | acceptance |
|---|---|---|
| trained | 2.1% (21 of 1000) | 59.3% |
| random weights | 74.8% (748 of 1000) | 0.9% |
| local `BondSwap` | 0.2% (2 sites) | — |

A trained move is ~10× larger than a local swap *and* usually accepted; the
random one proposes near-total re-randomisation and is therefore always
rejected. `Z_R` (`[lattice].Z_R`, the latent step radius) trades move size
against acceptance — that trade is now measurable.

A cheap offline sanity check, before spending the WL run: reconstruction
accuracy the way `alloy.cc::decode` uses the model (argmax per site, then the
composition repair loop). We saw 93.6% for the trained model vs 25.0% —
exactly chance for four species — for random weights.

---

## Pitfalls

### The grid is padded

`ini_sys()` inflates an N³ lattice onto a padded cubic grid:

```
SHIFT = N - 1
VAE_D = 2*(N-1) + SHIFT + 1          # == 3N - 2
PAD   = (ceil(VAE_D/16)*16 - VAE_D) / 2    # truncating integer divide
VAE_D = VAE_D + 2*PAD
```

Each site goes at `coord + SHIFT + PAD`. So N=4 spans 10 cells but the engine
feeds a **16³** grid; N=10 spans 28 and feeds **32³**.

Sizing the grid from the coordinate range in an xyz file gives the *unpadded*
span and produces tensors the engine cannot use. `utils/geometry.py` mirrors
the rule; use it rather than re-deriving:

```bash
cd vae-modeling
python3 -c "from utils.geometry import engine_grid; print(engine_grid(10))"
# (32, 11)  -> grid 32, offset 11
```

The divide truncates, so `VAE_D` is **not** always a multiple of 16 — N=5 gives
15. Copy the formula; do not "round up to 16".

| N | 4 | 5 | 6 | 10 | 16 | 20 |
|---|---|---|---|----|----|----|
| grid | 16 | 15 | 16 | 32 | 48 | 64 |

### Element order must match

The one-hot channel is the index into `[lattice].elements` (the engine's
`element[]` table has a `"HEA"` sentinel at index 0, so `element[t+1]` is the
symbol for `Atom` value `t`). Pass `--elements` in that same order.

### A WL window must be reachable

`flatWL()` returns `min/avg` over masked bins, so the flatness loop cannot exit
until **every** masked bin has been visited at least once. A window whose lower
edge sits below the reachable ground state never converges, and nothing warns
you — the run simply spins.

Pick the window from a PT run's observed energies (`ptEmin`/`ptEmax`, printed
at startup) rather than guessing.

### Parallel tempering may not mix

Swap acceptance goes as `exp(-Δβ·ΔE)` and `ΔE` is extensive, so large lattices
on few replicas do not exchange at all. Our N=10 (1000 sites) run on 8 ranks
spanning two decades gave:

```
swap prob:  T=10: 2.0%   21: 1.3%   46: 1.5%   97: 0.5%
            T=207: 0     441: 0     938: 0
```

Zero above 207 K. With local moves also frozen at the cold end (`Rot prob`
0.09% at T=10 K vs 75% at T=2000 K), the two coldest replicas produced 78% and
58% duplicate frames — so the **ordered, low-energy configurations the θ gate
actually targets were the thinnest part of the training set**.

Check `misc0.dat` after collecting. If swap acceptance is near zero, the fixes
are more replicas (a denser ladder), a much larger `sep`, or
`[output] snapshot_lowe = true` during a WL run, which targets low-energy
configurations by construction.

### Snapshots are standard xyz

`write_xyz` emits exactly two header lines per frame — the atom count, then a
comment line carrying `Eng`, `MC_step` and the short-range-order table. Frames
accumulate in one file per rank per trajectory group, which is what
`create_vae_input.py` expects. VMD and ASE read them directly.

---

## See also

- [`CLAUDE.md`](../CLAUDE.md) — architecture, the paper ↔ code map, TOML reference
- [`vae-modeling/README.md`](../vae-modeling/README.md) — training stack details
- Yin, Wang, Shankar, *DeepThermo*, IPDPS 2023,
  [10.1109/IPDPS54959.2023.00041](https://doi.org/10.1109/IPDPS54959.2023.00041)
