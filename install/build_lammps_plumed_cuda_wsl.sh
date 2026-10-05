#!/usr/bin/env bash
# LAMMPS patch_10Sep2025 + PLUMED 2.9.4 with the GPU package (CUDA, mixed precision) in WSL2,
# for the machine with the RTX 4050 Laptop (Ada, sm_89, 6 GB).
#
# Prerequisites on that machine (once, needs sudo):
#   - NVIDIA Windows driver with WSL support (nvidia-smi must work inside WSL)
#   - CUDA toolkit inside WSL:  sudo apt install nvidia-cuda-toolkit   (or the toolkit from developer.nvidia.com)
#   - build tools:              sudo apt install build-essential cmake git openmpi-bin libopenmpi-dev
#
#   run:  wsl -d Ubuntu-22.04 -- bash /mnt/c/Users/ileao/OneDrive/Documentos/Nucleation_MD_GPU/install/build_lammps_plumed_cuda_wsl.sh
#
# Machine with the RTX 4050 (2026-10-04): CUDA 12.6 already in /usr/local/cuda, g++ 11, OpenMPI 4.1.2;
# Ubuntu-22.04 has no cmake -> a portable cmake is unpacked in ~/opt/cmake (no sudo).
set -euo pipefail

export PATH=/usr/local/cuda/bin:/usr/lib/wsl/lib:$HOME/opt/cmake/bin:$PATH
CMAKE_VER=3.31.6
if ! command -v cmake >/dev/null; then
  mkdir -p "$HOME/opt/cmake"
  curl -sSL "https://github.com/Kitware/CMake/releases/download/v$CMAKE_VER/cmake-$CMAKE_VER-linux-x86_64.tar.gz" \
    | tar xz -C "$HOME/opt/cmake" --strip-components=1
fi
cmake --version | head -1

TAG=patch_10Sep2025
SRC=$HOME/opt/src/lammps-10Sep2025
BUILD=$SRC/build-cuda
PREFIX=$HOME/opt/lammps-10Sep2025-cuda
VENV=$HOME/opt/venvs/nucmd
NPROC=${NPROC:-8}   # WSL has ~7 GB RAM here; nvcc with -j20 risks OOM

command -v nvidia-smi >/dev/null || { echo "nvidia-smi not found in WSL: install/update the NVIDIA Windows driver"; exit 1; }
command -v nvcc >/dev/null || { echo "nvcc not found: install the CUDA toolkit inside WSL"; exit 1; }
nvidia-smi --query-gpu=name,memory.total,compute_cap --format=csv

[ -d "$SRC" ] || git clone --depth 1 --branch "$TAG" https://github.com/lammps/lammps.git "$SRC"

if [ ! -x "$VENV/bin/pip" ]; then
  python3 -m venv --without-pip "$VENV"
  curl -sSL https://bootstrap.pypa.io/get-pip.py | "$VENV/bin/python"
fi
"$VENV/bin/pip" install -q numpy scipy matplotlib pandas pytest

mkdir -p "$BUILD"
cd "$BUILD"
cmake ../cmake \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_INSTALL_PREFIX="$PREFIX" \
  -D BUILD_MPI=yes -D BUILD_OMP=yes -D BUILD_SHARED_LIBS=yes -D LAMMPS_EXCEPTIONS=yes \
  -D PKG_GPU=yes -D GPU_API=cuda -D GPU_PREC=mixed -D GPU_ARCH=sm_89 \
  -D PKG_MANYBODY=yes -D PKG_MEAM=yes \
  -D PKG_EXTRA-FIX=yes -D PKG_EXTRA-COMPUTE=yes -D PKG_EXTRA-DUMP=yes -D PKG_EXTRA-PAIR=yes \
  -D PKG_MISC=yes -D PKG_MC=yes -D PKG_REPLICA=yes -D PKG_MOLECULE=yes \
  -D PKG_COLVARS=yes -D PKG_OPENMP=yes \
  -D PKG_PLUMED=yes -D DOWNLOAD_PLUMED=yes -D PLUMED_MODE=runtime \
  -D Python_EXECUTABLE="$VENV/bin/python"
cmake --build . -j "$NPROC"
cmake --install .
"$VENV/bin/pip" install -q "$SRC/python"

PLUMED_PREFIX=$(find "$BUILD" -maxdepth 2 -type d -name 'plumed_build-prefix' | head -1)
cat > "$HOME/opt/nucmd_cuda_env.sh" <<EOF
# source ~/opt/nucmd_cuda_env.sh
export PATH=$PREFIX/bin:$PLUMED_PREFIX/bin:/usr/local/cuda/bin:\$PATH
export LD_LIBRARY_PATH=$PREFIX/lib:$PLUMED_PREFIX/lib:/usr/local/cuda/lib64:/usr/lib/wsl/lib:\${LD_LIBRARY_PATH:-}
export PLUMED_KERNEL=$PLUMED_PREFIX/lib/libplumedKernel.so
export OMP_NUM_THREADS=1
source $VENV/bin/activate
EOF
echo "BUILD OK -> source ~/opt/nucmd_cuda_env.sh ; lmp -sf gpu -pk gpu 1 -in ..."
