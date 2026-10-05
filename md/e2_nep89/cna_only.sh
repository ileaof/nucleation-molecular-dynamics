#!/usr/bin/env bash
# Recompute the CNA fcc fraction of existing GPUMD runs (no new dynamics):  bash cna_only.sh coexist_980 ...
source ~/opt/nucmd_cuda_env.sh
cd "$(dirname "$0")"
for d in "$@"; do
  python nep_tools.py todump $d/prod.xyz $d/prod.dump $d/frame0.data > /dev/null
  lmp -in in.cna_rerun -var d $d -log $d/cna.log -screen none; echo "$d cna exit=$?"
  rm -f $d/prod.dump
done
