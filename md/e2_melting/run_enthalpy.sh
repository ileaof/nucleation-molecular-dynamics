#!/usr/bin/env bash
# Latent heat at T_m for the Borovikov 2024 potential: <H>_liquid - <H>_solid (NPT, P = 0).
source ~/opt/nucmd_env.sh
cd "$(dirname "$0")"
run() { mpirun --use-hwthread-cpus -np 2 lmp -in in.enthalpy -var pot borovikov -var T 935 -var phase $1 -var n 8 \
        -var neq 10000 -var nprod 50000 -log H_borovikov_$1.log > /dev/null 2>&1; echo "$1 exit=$?" >> H_status.txt; }
rm -f H_status.txt
run solid & run liquid & wait
echo DONE >> H_status.txt
