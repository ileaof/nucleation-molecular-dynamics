#!/usr/bin/env bash
source ~/opt/nucmd_env.sh
cd "$(dirname "$0")"
mpirun --use-hwthread-cpus -np 4 lmp -in in.cfm -var pot borovikov -var Tm 935 -var orient 100 \
   -var nstage 1000 -var neq 1000 -var nprod 6000 -var ndump 500 -log cfm_test.log > /dev/null 2>&1
echo exit=$?; grep -E "ERROR|Performance" cfm_test.log | tail -2
