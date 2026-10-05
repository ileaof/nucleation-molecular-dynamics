#!/usr/bin/env bash
# Refine T_m (pure Al, Al-Nb MEAM 2022): isothermal coexistence at 940/945/950 K, 3 jobs x 4 ranks.
source ~/opt/nucmd_env.sh
cd "$(dirname "$0")"
run() { mpirun --use-hwthread-cpus -np 4 lmp -in in.coexist -var pot alnb -var T0 $1 -var ens npt -var nx 6 -var nz 24 \
        -var neq 5000 -var nprod 50000 -log bracket_alnb_$1.log > /dev/null 2>&1; echo "alnb $1 exit=$?" >> refine_status.txt; }
rm -f refine_status.txt
run 940 & run 945 & run 950 & wait
echo DONE >> refine_status.txt
