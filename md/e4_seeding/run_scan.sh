#!/usr/bin/env bash
# E4 seeding scan: seed radius R [Å] at several T; one replica per (R, T) unless SEEDS is set.
#   wsl -d Ubuntu-22.04 -- bash .../md/e4_seeding/run_scan.sh <R> <T1> <T2> ...      (SEEDS="1 2 3" for replicas)
source ~/opt/nucmd_cuda_env.sh
cd "$(dirname "$0")"
R=$1; shift
mkdir -p logs
for T in "$@"; do for s in ${SEEDS:-1}; do
  f=logs/seed_R${R}_T${T}_s${s}.log
  [ -s $f ] && grep -q "Total wall time" $f && continue
  lmp -sf gpu -pk gpu 1 -in in.seed -var R $R -var T $T -var seed $s -log $f > /dev/null 2>&1
  echo "R=$R T=$T s=$s exit=$?" >> status.txt
done; done
echo "DONE R=$R" >> status.txt
