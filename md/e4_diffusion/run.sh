#!/usr/bin/env bash
source ~/opt/nucmd_cuda_env.sh
cd "$(dirname "$0")"
for T in 813 851 936; do lmp -sf gpu -pk gpu 1 -in in.msd -var T $T -log log_$T.lammps > /dev/null 2>&1; echo "T=$T exit=$?"; done
