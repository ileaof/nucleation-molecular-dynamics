#!/usr/bin/env bash
# Small systems (E2 size: n=13 -> 8788 atoms, n=20 -> 32000 atoms): GPU vs CPU (8 MPI ranks) on the RTX 4050 machine.
#   wsl -d Ubuntu-22.04 -- bash /mnt/c/Users/ileao/OneDrive/Documentos/Nucleation_MD_GPU/md/gpu_benchmark/run_small_wsl.sh
source ~/opt/nucmd_cuda_env.sh
cd "$(dirname "$0")"
for n in 13 20; do
  lmp -sf gpu -pk gpu 1 -in in.bench -var n $n -var steps 5000 -log bench_cuda_n$n.log -screen none
  mpirun -np 8 lmp -in in.bench -var n $n -var steps 5000 -log bench_cpu8_n$n.log -screen none
done
grep -H "^Performance" bench_*_n13.log bench_*_n20.log
