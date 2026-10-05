#!/usr/bin/env bash
# T_m bracketing (pure Al), EAM potentials, 8192 atoms, NPT coexistence 100 ps, 3 jobs x 4 ranks.
source ~/opt/nucmd_env.sh
cd "$(dirname "$0")"
run() { mpirun --use-hwthread-cpus -np 4 lmp -in in.coexist -var pot $1 -var T0 $2 -var ens npt -var nx 8 -var nz 32 \
        -var neq 5000 -var nprod 50000 -log bracket_$1_$2.log > /dev/null 2>&1; echo "$1 $2 exit=$?" >> eam_status.txt; }
rm -f eam_status.txt
jobs_list=( "borovikov 900" "borovikov 930" "borovikov 960" "borovikov 990" "mendelev 900" "mendelev 920" "mendelev 940" "mendelev 960" )
i=0
for j in "${jobs_list[@]}"; do
  run $j &
  i=$((i+1)); if (( i % 3 == 0 )); then wait; fi
done
wait
echo DONE >> eam_status.txt
