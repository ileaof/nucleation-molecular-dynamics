#!/usr/bin/env bash
# interface stress series: "ex ey seed" per argument; T = T_m of Borovikov (936 K, E2)
#   wsl -d Ubuntu-22.04 -- bash .../md/e5_surface_stress/run_istress.sh "0 0 1" "0.01 0 1" ...
source ~/opt/nucmd_cuda_env.sh
cd "$(dirname "$0")"
T=${T:-936}; NPROD=${NPROD:-500000}
for spec in "$@"; do set -- $spec
  tag=T${T}_ex$1_ey$2_s$3
  [ -s glob_$tag.txt ] && grep -q "Total wall time" log_$tag.lammps 2>/dev/null && continue
  lmp -sf gpu -pk gpu 1 -in in.istress -var T $T -var ex $1 -var ey $2 -var seed $3 -var nprod $NPROD -var tag $tag -log log_$tag.lammps > /dev/null 2>&1
  echo "$tag exit=$?" >> status.txt
done
echo DONE >> status.txt
