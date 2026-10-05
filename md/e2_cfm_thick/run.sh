#!/usr/bin/env bash
# γ₀ robustness (HANDOFF §4, option 2): Borovikov 2024 CFM with thicker ribbons (W ≈ 28–39 Å vs 12–20 Å in E2)
# and 1 ns production (vs 400 ps). Same protocol as md/e2_melting/in.cfm; LAMMPS GPU (eam/fs/gpu).
# Dumps go to WSL native disk (~/runs/cfm_thick, ~0.5–1 GB each) to keep them out of OneDrive.
#   wsl -d Ubuntu-22.04 -- bash .../md/e2_cfm_thick/run.sh
source ~/opt/nucmd_cuda_env.sh
cd "$(dirname "$0")"
OUT=~/runs/cfm_thick; mkdir -p $OUT
rm -f status.txt
for spec in "100 8" "110 5" "111 6"; do
  set -- $spec
  lmp -sf gpu -pk gpu 1 -in in.cfm -var pot borovikov -var Tm 936 -var orient $1 -var ny $2 \
      -var nprod 500000 -var ndump 4000 -var dump $OUT/cfm_borovikov_$1.dump \
      -log cfm_borovikov_$1.log > /dev/null 2>&1
  echo "$1 ny=$2 exit=$?" >> status.txt
done
echo DONE >> status.txt
