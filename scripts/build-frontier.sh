#!/bin/bash
# Build DeepThermo on OLCF Frontier (AMD MI250X, ROCm).
#
# Usage: ./scripts/build-frontier.sh [cmake-args...]
#
# By default builds with the LibTorch (ROCm) backend. Override with
#   BACKEND=tf ./scripts/build-frontier.sh -DTF_DIR=...
#   BACKEND=redis ./scripts/build-frontier.sh -DREDIS_DIR=...
set -euo pipefail

BACKEND="${BACKEND:-torch}"
BUILD_DIR="${BUILD_DIR:-build-frontier}"

module load PrgEnv-gnu
module load rocm
module load cmake

SOURCE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$SOURCE_DIR"

if [[ "$BACKEND" == "torch" ]]; then
    : "${TORCH_INSTALL_PREFIX:?Set TORCH_INSTALL_PREFIX to the LibTorch-ROCm install root}"
    cmake -B "$BUILD_DIR" \
        -DCMAKE_PREFIX_PATH="$TORCH_INSTALL_PREFIX" \
        -DCMAKE_CXX_COMPILER=CC \
        -DDEEPTHERMO_BACKEND=torch \
        -DDEEPTHERMO_PLATFORM=frontier \
        "$@"
else
    cmake -B "$BUILD_DIR" \
        -DCMAKE_CXX_COMPILER=CC \
        -DDEEPTHERMO_BACKEND="$BACKEND" \
        -DDEEPTHERMO_PLATFORM=frontier \
        "$@"
fi

cmake --build "$BUILD_DIR" --parallel
echo "built: $BUILD_DIR/hea-wl"
