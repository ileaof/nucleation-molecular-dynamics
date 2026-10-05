#!/usr/bin/env bash
source ~/opt/nucmd_cuda_env.sh
cd "$(dirname "$0")"
for o in 100 110 111; do
  mkdir -p cfm_$o
  lmp -in in.build_cfm -var orient $o -log none -screen none
  python nep_tools.py fromdata cfm_$o
done
