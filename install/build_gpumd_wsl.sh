#!/usr/bin/env bash
# GPUMD (native NEP/NEP89 MD on the GPU) for the RTX 4050 (sm_89) in WSL Ubuntu-22.04. No sudo.
#   wsl -d Ubuntu-22.04 -- bash /mnt/c/Users/ileao/OneDrive/Documentos/Nucleation_MD_GPU/install/build_gpumd_wsl.sh
# Result: ~/opt/src/GPUMD/src/{gpumd,nep}; ~/opt/nucmd_cuda_env.sh gets GPUMD on the PATH.
set -euo pipefail
export PATH=/usr/local/cuda/bin:$PATH
SRC=$HOME/opt/src/GPUMD
NPROC=${NPROC:-8}   # WSL has ~7 GB RAM here

[ -d "$SRC" ] || git clone --depth 1 https://github.com/brucefan1983/GPUMD.git "$SRC"
git -C "$SRC" log -1 --format="GPUMD commit %h %cd"
make -C "$SRC/src" -j "$NPROC" CUDA_ARCH=-arch=sm_89
ENV=$HOME/opt/nucmd_cuda_env.sh
grep -q GPUMD "$ENV" 2>/dev/null || echo "export PATH=$SRC/src:\$PATH   # GPUMD (gpumd, nep)" >> "$ENV"
echo "GPUMD OK -> $SRC/src/gpumd"
