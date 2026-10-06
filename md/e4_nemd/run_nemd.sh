#!/usr/bin/env bash
# E4 NEMD runs: "R T0 G seed" per argument (G in K/Å). Outputs in ~/runs/e4_nemd/<tag>/ (WSL disk; T3d/dumps are large);
# logs and the 1D profile are copied back to md/e4_nemd/runs/<tag>/.
#   wsl -d Ubuntu-22.04 -- bash .../md/e4_nemd/run_nemd.sh "20 809 0.0 1" "20 809 0.24 1"
source ~/opt/nucmd_cuda_env.sh
cd "$(dirname "$0")"
for spec in "$@"; do
  set -- $spec; R=$1; T0=$2; G=$3; s=$4
  tag=R${R}_T${T0}_G${G}_s${s}; O=~/runs/e4_nemd/$tag; mkdir -p $O runs/$tag
  if ! grep -q "Total wall time" $O/log.lammps 2>/dev/null; then
    lmp -sf gpu -pk gpu 1 -in in.nemd -var R $R -var T0 $T0 -var G $G -var seed $s -var out $O -log $O/log.lammps > /dev/null 2>&1
    echo "$tag exit=$?" >> status.txt
  fi
  cp $O/log.lammps $O/Tz.txt runs/$tag/
done
echo "DONE $*" >> status.txt
