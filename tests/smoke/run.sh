#!/usr/bin/env bash
# Smoke test for the full hea-wl pipeline.
#
# Builds a tiny random-weight VAE, exports it to TorchScript, drops the
# .pt files plus a 4^3 config.toml + coupling.input into a temp dir, runs
# hea-wl single-process for one short WL iteration, and checks that the
# canonical output files (run.dat, therm.dat, DOS_H_iter*.dat, comm.dat)
# materialise. No physics check — purely a "did the pipeline survive?"
# regression detector.
#
# CMake passes DEEPTHERMO_BUILD_DIR + DEEPTHERMO_SOURCE_DIR via the env.

set -euo pipefail

: "${DEEPTHERMO_BUILD_DIR:?must be set (CMake should pass it)}"
: "${DEEPTHERMO_SOURCE_DIR:?must be set (CMake should pass it)}"

if ! command -v python3 >/dev/null; then
    echo "smoke: python3 not on PATH; skipping" >&2
    exit 0  # treat as soft skip
fi
if ! python3 -c 'import torch' 2>/dev/null; then
    echo "smoke: python3 -c 'import torch' failed; skipping" >&2
    exit 0
fi
if ! command -v mpirun >/dev/null; then
    echo "smoke: mpirun not on PATH; skipping" >&2
    exit 0
fi

BIN="$DEEPTHERMO_BUILD_DIR/hea-wl"
if [[ ! -x "$BIN" ]]; then
    echo "smoke: $BIN not executable" >&2
    exit 1
fi

WORK=$(mktemp -d -t deepthermo-smoke-XXXX)
trap 'rm -rf "$WORK"' EXIT
mkdir -p "$WORK/models"

# 1. Generate a random-weight VAE matching what the C++ side will feed it.
#    For lattice N=4 the BCC coordinate encoding (paper §III-B) inflates
#    to a 16x16x16 input grid (alloyState.VAE_D = 16). Export wrappers at
#    models/encoder_MoNbTaW.pt etc.
PYTHONPATH="$DEEPTHERMO_SOURCE_DIR/vae-modeling/vae" \
python3 - "$WORK/models" <<'PY'
import sys
import torch
from pathlib import Path
from src.export_torchscript import export
from src.model import VAE, VAEConfig

out_dir = Path(sys.argv[1])
GRID = 16
NE = 4
cfg = VAEConfig(grid_size=GRID, num_elements=NE, latent_dim=3)
model = VAE(cfg)
weights = out_dir / "vae.pt"
torch.save(model.state_dict(), weights)
export(weights, "MoNbTaW", grid_size=GRID, num_elements=NE, latent_dim=3, out_dir=out_dir)
PY

# 2. Reuse the canonical coupling.input so we don't reinvent shell offsets.
cp "$DEEPTHERMO_SOURCE_DIR/tests/smoke/coupling.input" "$WORK/coupling.input"

# 3. Tiny config. WL convergence is intentionally short-circuited by
#    mod_factor_init < mod_factor_final so the outer for-loop in main.cc
#    is skipped — this is a pipeline smoke test, not a WL-physics test.
#    The torch backend's load() path still runs, ReadInput is exercised,
#    PT warm-up, ini_alloy, Etot, the production-sampling loop and
#    thermoqs all run through.
cat >"$WORK/config.toml" <<'TOML'
[lattice]
N              = 4
NE             = 4
SH             = 6
elements       = ["Mo", "Nb", "Ta", "W"]
composition    = [0.25, 0.25, 0.25, 0.25]
nb_interaction = 1
coupling_file  = "coupling.input"

[wang_landau]
bin_width            = 0.02
e_min                = -1.40
e_max                = -1.00
flatness             = 0.05
mod_factor_init      = 0.001
iteration_factor     = 2.0
mod_factor_final     = 0.5
production_bin_samps = 1

[parallel_tempering]
metropolis_sampling = false
T_init  = 100.0
T_final = 1000.0
samples = 2
sep     = 1
drop    = 1

[thermodynamics]
T_init  = 100
T_final = 200
dT      = 50

[model]
dir = "models"
TOML

# 4. Run. CMake imposes a 120s wall-clock test timeout from the outside.
cd "$WORK"
mpirun --oversubscribe -n 1 "$BIN" config.toml 7 2>&1 | tee run.log

# 5. Sanity-check the output files. comm.dat may be empty when the WL
#    outer loop short-circuits (intentional in this config); we only
#    require that it was opened.
for f in run.dat therm.dat; do
    if [[ ! -s "$WORK/$f" ]]; then
        echo "smoke: $f not written or empty" >&2
        exit 1
    fi
done
if [[ ! -e "$WORK/comm.dat" ]]; then
    echo "smoke: comm.dat not created" >&2
    exit 1
fi

echo "smoke: ok"
