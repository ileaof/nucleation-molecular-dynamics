#!/usr/bin/env bash
source ~/opt/nucmd_env.sh
cd "$(dirname "$0")"
for p in alnb jelinek; do for dt in 0.002 0.001; do
  mpirun --use-hwthread-cpus -np 8 lmp -in in.nve_drift -var pot $p -var dt $dt -log drift_${p}_${dt}.log >/dev/null 2>&1
  awk -v p=$p -v dt=$dt '/Step/{f=1;next} /Loop/{f=0} f&&NF==4{t[++n]=$2;T[n]=$3;E[n]=$4} END{printf "%s dt=%s  T_end=%.0f K  Etot drift = %+.4f eV/atom/ns  (E0=%.3f)\n", p, dt, T[n], (E[n]-E[int(n/4)])/2048/((t[n]-t[int(n/4)])/1000), E[1]}' drift_${p}_${dt}.log
done; done
