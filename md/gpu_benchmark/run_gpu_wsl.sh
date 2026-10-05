#!/usr/bin/env bash
# Benchmark on the RTX 4050 with the WSL CUDA build (Ubuntu-22.04) + CPU reference on the same machine.
#   wsl -d Ubuntu-22.04 -- bash /mnt/c/Users/ileao/OneDrive/Documentos/Nucleation_MD_GPU/md/gpu_benchmark/run_gpu_wsl.sh
source ~/opt/nucmd_cuda_env.sh
cd "$(dirname "$0")"
for n in 40 63 100; do
  lmp -sf gpu -pk gpu 1 -in in.bench -var n $n -log bench_cuda_n$n.log -screen none
done
mpirun -np 8 lmp -in in.bench -var n 40 -var steps 500 -log bench_cpu8_n40.log -screen none
grep -H Performance bench_*.log
