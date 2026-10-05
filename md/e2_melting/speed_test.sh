#!/usr/bin/env bash
source ~/opt/nucmd_env.sh
cd "$(dirname "$0")"
for p in alnb jelinek; do
  for np in 8 16; do
    mpirun --use-hwthread-cpus -np $np lmp -in in.coexist -var pot $p -var T0 930 -var neq 500 -var nprod 1000 -log speed_${p}_np$np.log > /dev/null 2>&1
    echo "$p np=$np exit=$?  $(grep -E 'Performance' speed_${p}_np$np.log | tail -1)"
    grep ERROR speed_${p}_np$np.log
  done
done
rm -f solid_*.data coexist_*.data coexist_*.dump
