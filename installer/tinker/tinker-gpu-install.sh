#!/usr/bin/env bash
# Tinker9-GPU installer — standalone, for advanced users.
# Target: Ubuntu 22.04 WSL2 + NVIDIA GPU (default: RTX 4070, CC 8.9)
# Repo:   https://github.com/TinkerTools/tinker-gpu
#
# Usage:  bash installer/tinker/tinker-gpu-install.sh
# Env overrides:
#   COMPUTE_CAPABILITY=89   GPU compute capability × 10 (89 = RTX 4070)
#   PREFIX=~/tinker9        Installation prefix (no sudo needed)
#   SRC_ROOT=/tmp/tinker9   Build directory root
#   JOBS=$(nproc)           Parallel build jobs
#   BUILD_TYPE=Release      CMake build type
#   PREC=mixed              Floating-point precision: double/single/mixed
set -euo pipefail

REPO_URL="https://github.com/TinkerTools/tinker-gpu"
SRC_ROOT="${SRC_ROOT:-/tmp/tinker9-build}"
REPO_DIR="${SRC_ROOT}/tinker-gpu"
BUILD_DIR="${SRC_ROOT}/build"
PREFIX="${PREFIX:-$HOME/tinker9}"
HOME_BIN="$HOME/bin"
PROFILE_FILE="$HOME/.tinker9_env"
BUILD_TYPE="${BUILD_TYPE:-Release}"
JOBS="${JOBS:-$(nproc)}"
COMPUTE_CAPABILITY="${COMPUTE_CAPABILITY:-89}"   # RTX 4070 = 8.9 → 89
PREC="${PREC:-mixed}"

# ── 1. WSL check ──────────────────────────────────────────────────────────────
if ! grep -qi microsoft /proc/version 2>/dev/null; then
  echo "Warning: not running in WSL2. This script targets Ubuntu WSL2 + NVIDIA GPU."
fi

# ── 2. System packages ────────────────────────────────────────────────────────
MISSING_PKGS=()
for pkg in build-essential git cmake ninja-build gfortran gcc g++ make pkg-config ca-certificates curl wget libfftw3-dev; do
  dpkg -s "$pkg" >/dev/null 2>&1 || MISSING_PKGS+=("$pkg")
done
if [ ${#MISSING_PKGS[@]} -gt 0 ]; then
  echo "==> Installing missing system packages: ${MISSING_PKGS[*]}"
  sudo apt-get update -qq
  sudo apt-get install -y "${MISSING_PKGS[@]}"
else
  echo "==> All system build dependencies already installed"
fi

# ── 3. CUDA toolkit — auto-install if nvcc is missing ─────────────────────────
if ! command -v nvcc >/dev/null 2>&1; then
  echo "==> nvcc not found — attempting CUDA toolkit installation"
  UBUNTU_VER=$(. /etc/os-release && echo "${VERSION_ID/./}")   # e.g. 2204
  DEB_URL="https://developer.download.nvidia.com/compute/cuda/repos/ubuntu${UBUNTU_VER}/x86_64/cuda-keyring_1.1-1_all.deb"
  PIN_URL="https://developer.download.nvidia.com/compute/cuda/repos/ubuntu${UBUNTU_VER}/x86_64/cuda-ubuntu${UBUNTU_VER}.pin"
  wget -qO /tmp/cuda-pin "$PIN_URL" && sudo mv /tmp/cuda-pin /etc/apt/preferences.d/cuda-repository-pin-600
  wget -qO /tmp/cuda-keyring.deb "$DEB_URL" && sudo dpkg -i /tmp/cuda-keyring.deb
  sudo apt-get update -qq
  sudo apt-get install -y cuda-toolkit
  # Add CUDA to PATH for this session
  for d in /usr/local/cuda/bin /usr/local/cuda-*/bin; do
    [ -d "$d" ] && export PATH="$d:$PATH" && break
  done
fi

if ! command -v nvcc >/dev/null 2>&1; then
  echo "ERROR: nvcc still not found after installation attempt."
  echo "Install the NVIDIA CUDA Toolkit manually: https://developer.nvidia.com/cuda-downloads"
  exit 1
fi
echo "==> CUDA compiler: $(nvcc --version | head -1)"

# ── 4. nvidia-smi (WSL path fallback) ────────────────────────────────────────
# Add WSL lib path early so nvidia-smi is found even if not on default PATH
[ -d /usr/lib/wsl/lib ] && export PATH="/usr/lib/wsl/lib:$PATH"
if ! command -v nvidia-smi >/dev/null 2>&1; then
  echo "ERROR: nvidia-smi not found. Install/update the Windows NVIDIA driver with WSL2 CUDA support."
  exit 1
fi
echo "==> GPU:"
nvidia-smi --query-gpu=name,compute_cap --format=csv,noheader 2>/dev/null || nvidia-smi | head -12
echo "==> Compute capability target: ${COMPUTE_CAPABILITY} ($(echo "$COMPUTE_CAPABILITY" | sed 's/\(.\)$/.\1/'))"

# ── 5. Locate CUDA directory ──────────────────────────────────────────────────
CUDA_DIR="${CUDA_HOME:-}"
if [ -z "$CUDA_DIR" ]; then
  if [ -d /usr/local/cuda ]; then
    CUDA_DIR=/usr/local/cuda
  else
    CUDA_DIR="$(dirname "$(dirname "$(command -v nvcc)")")"
  fi
fi
export CUDACXX="${CUDACXX:-$(command -v nvcc)}"
echo "==> CUDA_DIR: $CUDA_DIR"

# ── 6. Detect NVHPC (optional — enables OpenACC path) ────────────────────────
if command -v nvc++ >/dev/null 2>&1; then
  GPU_LANG="OPENACC"
  export ACC="${ACC:-$(command -v nvc++)}"
  echo "==> NVHPC detected → OpenACC + CUDA build (ACC=$ACC)"
else
  GPU_LANG="CUDA"
  echo "==> NVHPC not found → CUDA-only build"
fi

# ── 7. Clone / update repository ─────────────────────────────────────────────
mkdir -p "$SRC_ROOT"
if [ -d "$REPO_DIR/.git" ]; then
  echo "==> Updating existing clone at $REPO_DIR"
  git -C "$REPO_DIR" checkout -- CMakeLists.txt src/cudart/gpucard.cpp src/hippo/edisp.cpp include/ff/spatial.h 2>/dev/null || true
  git -C "$REPO_DIR" pull --ff-only
else
  echo "==> Cloning tinker-gpu"
  git clone "$REPO_URL" "$REPO_DIR"
fi

echo "==> Initializing submodules (CPU Tinker source is required)"
git -C "$REPO_DIR" submodule update --init --recursive

# Patch out the test subdirectory — its Catch2 submodule doesn't compile on GCC 13+
# and the tests are not needed for the wdrop pipeline executables.
sed -i 's/^add_subdirectory (test)/# add_subdirectory (test)  # disabled: Catch2 GCC13 incompatibility/' \
  "$REPO_DIR/CMakeLists.txt"

# Patch gpucard.cpp for CUDA 13+ API compatibility:
# singleToDoublePrecisionPerfRatio, computeMode, and clockRate were removed from cudaDeviceProp.
export GPUCARD="$REPO_DIR/src/cudart/gpucard.cpp"
python3 - <<'PYEOF'
import os
path = os.environ['GPUCARD']
with open(path) as f:
    code = f.read()

# singleToDoublePrecisionPerfRatio removed in CUDA 13
code = code.replace(
    'a.single_double_ratio = prop.singleToDoublePrecisionPerfRatio;',
    '// singleToDoublePrecisionPerfRatio removed in CUDA 13\n   a.single_double_ratio = 0;'
)

# computeMode removed in CUDA 13 — replace the entire if/else chain
old_mode = (
    '   if (prop.computeMode == cudaComputeModeExclusive)\n'
    '      a.compute_mode_string = "exclusive thread";\n'
    '   else if (prop.computeMode == cudaComputeModeProhibited)\n'
    '      a.compute_mode_string = "prohibited";\n'
    '   else if (prop.computeMode == cudaComputeModeExclusiveProcess)\n'
    '      a.compute_mode_string = "exclusive process";\n'
    '   else\n'
    '      a.compute_mode_string = "default";'
)
new_mode = '   // computeMode removed in CUDA 13\n   a.compute_mode_string = "default";'
code = code.replace(old_mode, new_mode)

# clockRate removed in CUDA 13 — use cudaDeviceGetAttribute instead
code = code.replace(
    'a.clock_rate_kHz = prop.clockRate;',
    '{ int _clk = 0; cudaDeviceGetAttribute(&_clk, cudaDevAttrClockRate, device); a.clock_rate_kHz = _clk; }'
)

with open(path, 'w') as f:
    f.write(code)
print("gpucard.cpp patched OK")
PYEOF

# Patch hippo/edisp.cpp: 3 uses of "if CONSTEXPR" with runtime int values break in C++17.
# (CONSTEXPR macro becomes "constexpr" in C++17 but do_e/do_v are runtime, not constexpr.)
EDISP="$REPO_DIR/src/hippo/edisp.cpp"
sed -i \
  -e 's/if CONSTEXPR (do_e || do_v)/if (do_e || do_v)/g' \
  -e 's/if CONSTEXPR (do_e)/if (do_e)/g' \
  -e 's/if CONSTEXPR (do_v)/if (do_v)/g' \
  "$EDISP"
echo "edisp.cpp patched OK"

# Patch include/ff/spatial.h: increase LSTCAP from 48 to 128.
# The default 48 causes a runtime crash with proteins that have dense neighbor lists.
sed -i 's/static constexpr int LSTCAP = 48;/static constexpr int LSTCAP = 128;/' \
  "$REPO_DIR/include/ff/spatial.h"
echo "spatial.h LSTCAP patched OK (48 → 128)"

# ── 8. Configure CMake ────────────────────────────────────────────────────────
rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR"

CMAKE_ARGS=(
  -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
  -DCMAKE_INSTALL_PREFIX="$PREFIX"
  -DGPU_LANG="$GPU_LANG"
  -DCOMPUTE_CAPABILITY="$COMPUTE_CAPABILITY"
  -DPREC="$PREC"
  -DCUDA_DIR="$CUDA_DIR"
  -DBUILD_TESTING=OFF
  -DSTD=17
)

echo "==> Configuring CMake"
cmake -S "$REPO_DIR" -B "$BUILD_DIR" "${CMAKE_ARGS[@]}"

# ── 9. Build ──────────────────────────────────────────────────────────────────
echo "==> Building Tinker9-GPU with $JOBS parallel jobs"
echo "    (this may take 10–30 minutes)"
cmake --build "$BUILD_DIR" -j "$JOBS" --target tinker9

# ── 10. Install to PREFIX ─────────────────────────────────────────────────────
echo "==> Installing to $PREFIX"
mkdir -p "$PREFIX"
# Use sudo only if PREFIX is not writable by current user
if [ -w "$PREFIX" ]; then
  cmake --install "$BUILD_DIR"
else
  sudo cmake --install "$BUILD_DIR"
fi

# ── 11. Copy executables to ~/bin ─────────────────────────────────────────────
mkdir -p "$HOME_BIN"
# tinker-gpu installs to $PREFIX/gpu-m/ (not bin/).
# We copy them as real files (cp -L) so ~/bin is self-contained.
INSTALLED=()
for subdir in gpu-m bin .; do
  if [ -d "$PREFIX/$subdir" ]; then
    for exe in tinker9 analyze9 minimize9 dynamic9 bar9 testgrad9 info9; do
      src="$PREFIX/$subdir/$exe"
      if [ -e "$src" ]; then
        cp -fL "$src" "$HOME_BIN/$exe"
        chmod +x "$HOME_BIN/$exe"
        INSTALLED+=("$exe")
        echo "  ~/bin/$exe"
      fi
    done
    break
  fi
done
echo "==> Installed ${#INSTALLED[@]} executables to $HOME_BIN"

# ── 12. Persistent environment ────────────────────────────────────────────────
cat > "$PROFILE_FILE" <<EOF
export TINKER9_HOME="$PREFIX"
export PATH="\$PATH:$HOME/bin"
EOF

if ! grep -q 'source ~/.tinker9_env' "$HOME/.bashrc" 2>/dev/null; then
  printf '\n# Tinker9-GPU environment\nsource ~/.tinker9_env\n' >> "$HOME/.bashrc"
fi
# shellcheck disable=SC1090
source "$PROFILE_FILE"

# ── 13. Verify ────────────────────────────────────────────────────────────────
echo "==> Verifying installation"
FAIL=0
for exe in minimize9 dynamic9; do
  if command -v "$exe" >/dev/null 2>&1; then
    echo "  OK  $exe → $(command -v "$exe")"
  else
    echo "  ERR $exe — not found in PATH"
    FAIL=1
  fi
done
[ "$FAIL" -eq 1 ] && exit 1

echo
echo "Done! Tinker9-GPU installed at $PREFIX"
echo "  ~/bin/tinker9:   $(command -v tinker9   2>/dev/null || echo not-found)"
echo "  ~/bin/minimize9: $(command -v minimize9 2>/dev/null || echo not-found)"
echo "  ~/bin/dynamic9:  $(command -v dynamic9  2>/dev/null || echo not-found)"
echo "  ~/bin/analyze9:  $(command -v analyze9  2>/dev/null || echo not-found)"
echo
echo "To use GPU-accelerated Tinker in the wdrop pipeline:"
echo "  export TINKER_GPU=1"
echo "  ./pipeline/wdrop.sh input.pdb tinker SHORT"
