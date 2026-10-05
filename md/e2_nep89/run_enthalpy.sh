#!/usr/bin/env bash
# Enthalpy and density of solid and liquid Al, NEP89 (GPUMD), protocol of md/e2_melting/in.enthalpy:
# 8 x 8 x 8 cells (2048 atoms), NPT iso P = 0, dt = 2 fs; liquid first melted at 1700 K (10000 steps);
# 10000 steps at T, then 50000 production steps (100 ps), thermo every 100 steps.
#   wsl -d Ubuntu-22.04 -- bash .../md/e2_nep89/run_enthalpy.sh 950 798.15 817.75
source ~/opt/nucmd_cuda_env.sh
cd "$(dirname "$0")"
POT=$(realpath ../../potentials/NEP89_20250409/nep89_20250409.txt)
for T in "$@"; do for ph in solid liquid; do
  d=H_${ph}_$T; mkdir -p $d; rm -f $d/thermo.out
  python nep_tools.py bulk $d 4.05 8 > /dev/null
  {
    echo "potential $POT"
    if [ $ph = liquid ]; then
      echo "velocity 1700"; echo "time_step 2"
      echo "ensemble npt_mttk temp 1700 1700 iso 0 0 tperiod 100 pperiod 1000"; echo "run 10000"
    else
      echo "velocity $T"; echo "time_step 2"
    fi
    echo "ensemble npt_mttk temp $T $T iso 0 0 tperiod 100 pperiod 1000"; echo "run 10000"
    echo "ensemble npt_mttk temp $T $T iso 0 0 tperiod 100 pperiod 1000"
    echo "dump_thermo 100"; echo "dump_xyz 50000 final.xyz"; echo "run 50000"
  } > $d/run.in
  (cd $d && gpumd > gpumd.out 2>&1); echo "$ph $T exit=$?"
done; done
python analyze_h_nep.py
# final frame of every run by CNA: solid must stay fcc, liquid must not crystallise
for d in H_*_*; do
  [ -f $d/final.xyz ] || continue
  cp $d/final.xyz $d/prod.xyz
  python nep_tools.py todump $d/prod.xyz $d/prod.dump $d/frame0.data > /dev/null
  lmp -in in.cna_rerun -var d $d -log $d/cna.log -screen none
  echo "$d fcc fraction (final frame): $(awk '/v_ffcc/{f=1;next} /Loop/{f=0} f&&NF==2{v=$2} END{print v}' $d/cna.log)"
  rm -f $d/prod.dump $d/prod.xyz
done
