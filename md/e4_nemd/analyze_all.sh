#!/usr/bin/env bash
# Analyse every finished NEMD run in ~/runs/e4_nemd (optionally a glob, e.g. "R20_T809_*").
source ~/opt/nucmd_cuda_env.sh
cd "$(dirname "$0")"
for d in ~/runs/e4_nemd/${1:-R*}; do
  grep -q "Total wall time" $d/log.lammps 2>/dev/null || continue
  python analyze_nemd.py $d > runs/$(basename $d)/analysis.md 2>&1
  grep RESUMO runs/$(basename $d)/analysis.md
done
