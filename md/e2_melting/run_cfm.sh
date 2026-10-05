#!/usr/bin/env bash
# Capillary fluctuation method, Borovikov 2024, T_m = 936 K: (100), (110), (111) in parallel, 4 ranks each, 400 ps.
source ~/opt/nucmd_env.sh
cd "$(dirname "$0")"
run() { mpirun --use-hwthread-cpus -np 4 lmp -in in.cfm -var pot borovikov -var Tm 936 -var orient $1 \
        -var nprod 200000 -var ndump 2000 -log cfm_borovikov_$1.log > /dev/null 2>&1; echo "$1 exit=$?" >> cfm_status.txt; }
rm -f cfm_status.txt
run 100 & run 110 & run 111 & wait
echo DONE >> cfm_status.txt
