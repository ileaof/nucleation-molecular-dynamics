#!/usr/bin/env bash
# Δμ(Cu−Al) ↔ x_Cu calibration with SGC MCMD (GPUMD + NEP89), bulk solid and liquid, NPT P = 0.
# GPUMD SGC: swapping Al→Cu costs ΔE + μ_Cu − μ_Al; here μ_Al = 0, μ_Cu = Δμ.
#   wsl -d Ubuntu-22.04 -- bash .../md/e4_mcmd/run_dmu.sh <T> <dmu1> <dmu2> ...
source ~/opt/nucmd_cuda_env.sh
cd "$(dirname "$0")"
POT=$(realpath ../../potentials/NEP89_20250409/nep89_20250409.txt)
T=$1; shift
for dmu in "$@"; do for ph in solid liquid; do
  d=runs/dmu_${ph}_T${T}_${dmu}; mkdir -p $d; rm -f $d/*.out
  python ../e2_nep89/nep_tools.py bulk $d 4.05 8 > /dev/null
  {
    echo "potential $POT"
    if [ $ph = liquid ]; then
      echo "velocity 1500"; echo "time_step 2"
      echo "ensemble npt_mttk temp 1500 1500 iso 0 0 tperiod 100 pperiod 1000"; echo "run 10000"
    else
      echo "velocity $T"; echo "time_step 2"
    fi
    echo "ensemble npt_mttk temp $T $T iso 0 0 tperiod 100 pperiod 1000"; echo "run 5000"
    echo "ensemble npt_mttk temp $T $T iso 0 0 tperiod 100 pperiod 1000"
    echo "mc sgc 100 1000 $T $T 2 Al 0 Cu $dmu"
    echo "dump_thermo 500"; echo "dump_xyz 50000 final.xyz"; echo "run 50000"
  } > $d/run.in
  t0=$(date +%s); (cd $d && gpumd > gpumd.out 2>&1); echo "$ph T=$T dmu=$dmu exit=$? $(( $(date +%s)-t0 )) s"
done; done
