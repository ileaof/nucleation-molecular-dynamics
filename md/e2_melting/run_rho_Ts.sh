#!/usr/bin/env bash
# Densities of solid and (undercooled) liquid Al at the eutectic temperatures of the SDAS scripts.
source ~/opt/nucmd_env.sh
cd "$(dirname "$0")"
run() { mpirun --use-hwthread-cpus -np 2 lmp -in in.enthalpy -var pot borovikov -var T $1 -var phase $2 -var n 8 \
        -var neq 5000 -var nprod 15000 -log H_borovikov_$2_$1.log > /dev/null 2>&1; echo "$1 $2 exit=$?" >> rho_status.txt; }
rm -f rho_status.txt
run 817.75 solid & run 817.75 liquid & wait
run 798.15 solid & run 798.15 liquid & wait
echo DONE >> rho_status.txt
