#!/usr/bin/env bash
# CFM (111) rerun with a periodic ribbon (3 lattice cells along [11-2], W = 19.7 Å). 8 ranks once free.
source ~/opt/nucmd_env.sh
cd "$(dirname "$0")"
mpirun --use-hwthread-cpus -np 6 lmp -in in.cfm -var pot borovikov -var Tm 936 -var orient 111 \
       -var nprod 200000 -var ndump 2000 -log cfm_borovikov_111.log > /dev/null 2>&1
echo "111 exit=$?" > cfm111_status.txt
