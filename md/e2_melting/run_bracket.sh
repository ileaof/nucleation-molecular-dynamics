#!/usr/bin/env bash
# Isothermal bracketing of T_m: NPT coexistence at fixed T, 2 jobs x 4 ranks at a time.
source ~/opt/nucmd_env.sh
cd "$(dirname "$0")"
run() { mpirun --use-hwthread-cpus -np 4 lmp -in in.coexist -var pot $1 -var T0 $2 -var ens npt -var nx 6 -var nz 24 \
        -var neq 5000 -var nprod 50000 -log bracket_$1_$2.log > /dev/null 2>&1; echo "$1 $2 exit=$?"; }
for T in 900 930 960 990; do
  run alnb $T & run jelinek $T & wait
done
