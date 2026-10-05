#!/usr/bin/env bash
# T_m by solid-liquid coexistence (NPH) for both potentials, pure Al, 3456 atoms, 100 ps.
source ~/opt/nucmd_env.sh
cd "$(dirname "$0")"
for p in alnb jelinek; do
  mpirun --use-hwthread-cpus -np 8 lmp -in in.coexist -var pot $p -var T0 930 -var nx 6 -var nz 24 \
         -var neq 5000 -var nprod 50000 -log coexist_${p}_930.log > /dev/null 2>&1
  echo "$p exit=$?"
done
