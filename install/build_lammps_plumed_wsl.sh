#!/usr/bin/env bash
# Build LAMMPS patch_10Sep2025 (same version as the Windows binary) with PLUMED 2.9.4
# (version pinned by this LAMMPS release, MD5-checked) inside WSL Ubuntu-22.04.
# Everything goes under ~/opt; no sudo, system untouched.
#   run:  wsl -d Ubuntu-22.04 -- bash /mnt/c/Users/ileao/OneDrive/Documentos/Nucleation_MD_GPU/install/build_lammps_plumed_wsl.sh
set -euo pipefail

TAG=patch_10Sep2025
SRC=$HOME/opt/src/lammps-10Sep2025
BUILD=$SRC/build-plumed
PREFIX=$HOME/opt/lammps-10Sep2025
VENV=$HOME/opt/venvs/nucmd
NPROC=$(nproc)

[ -d "$SRC" ] || git clone --depth 1 --branch "$TAG" https://github.com/lammps/lammps.git "$SRC"

# Ubuntu's python3 lacks ensurepip (python3.10-venv needs sudo): bootstrap pip by hand
if [ ! -x "$VENV/bin/pip" ]; then
  rm -rf "$VENV"
  python3 -m venv --without-pip "$VENV"
  curl -sSL https://bootstrap.pypa.io/get-pip.py | "$VENV/bin/python"
fi
"$VENV/bin/pip" install -q --upgrade pip
"$VENV/bin/pip" install -q numpy scipy matplotlib pandas pytest

mkdir -p "$BUILD"
cd "$BUILD"
cmake ../cmake \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_INSTALL_PREFIX="$PREFIX" \
  -D BUILD_MPI=yes -D BUILD_OMP=yes \
  -D BUILD_SHARED_LIBS=yes \
  -D LAMMPS_EXCEPTIONS=yes \
  -D PKG_MANYBODY=yes -D PKG_MEAM=yes \
  -D PKG_EXTRA-FIX=yes -D PKG_EXTRA-COMPUTE=yes -D PKG_EXTRA-DUMP=yes -D PKG_EXTRA-PAIR=yes \
  -D PKG_MISC=yes -D PKG_MC=yes -D PKG_REPLICA=yes -D PKG_MOLECULE=yes \
  -D PKG_COLVARS=yes -D PKG_OPENMP=yes \
  -D PKG_PLUMED=yes -D DOWNLOAD_PLUMED=yes -D PLUMED_MODE=runtime \
  -D PKG_PYTHON=no \
  -D Python_EXECUTABLE="$VENV/bin/python"

cmake --build . -j "$NPROC"
cmake --install .

# LAMMPS python module into the venv (uses the shared liblammps just installed)
"$VENV/bin/pip" install -q "$SRC/python"

# PLUMED installed by the ExternalProject (CLI: sum_hills, driver, ...)
PLUMED_PREFIX=$(find "$BUILD" -maxdepth 2 -type d -name 'plumed_build-prefix' | head -1)

cat > "$HOME/opt/nucmd_env.sh" <<EOF
# source ~/opt/nucmd_env.sh
export PATH=$PREFIX/bin:$PLUMED_PREFIX/bin:\$PATH
export LD_LIBRARY_PATH=$PREFIX/lib:$PLUMED_PREFIX/lib:\${LD_LIBRARY_PATH:-}
export PLUMED_KERNEL=$PLUMED_PREFIX/lib/libplumedKernel.so  # runtime mode: LAMMPS loads this kernel
export OMP_NUM_THREADS=1
source $VENV/bin/activate
EOF
echo "BUILD OK -> source ~/opt/nucmd_env.sh"
