#!/usr/bin/env bash
# Liquid diffusion of Cu and Al in Al-3 wt% Cu, NEP89 (GPUMD): melt 1700 K (10 ps), NPT at T (20 ps), production NPT
# 200 ps with compute_msd per species (sample every 100 steps = 0.2 ps, 200 correlation steps = 40 ps window; 10 x 2000 was ~100x slower).
#   wsl -d Ubuntu-22.04 -- bash .../md/e4_diffusion/nep89/run_dl.sh 900 936.58 1000 1100 1200
source ~/opt/nucmd_cuda_env.sh
cd "$(dirname "$0")"
POT=$(realpath ../../../potentials/NEP89_20250409/nep89_20250409.txt)
for T in "$@"; do
  d=dl_T$T; python mk_alcu.py $d 16 3.0 > $d.log 2>&1 || { mkdir -p $d; python mk_alcu.py $d 16 3.0 > $d.log; }
  cat > $d/run.in <<EOR
potential $POT
velocity 1700
time_step 2
ensemble npt_mttk temp 1700 1700 iso 0 0 tperiod 100 pperiod 1000
run 5000
ensemble npt_mttk temp $T $T iso 0 0 tperiod 100 pperiod 1000
run 10000
ensemble npt_mttk temp $T $T iso 0 0 tperiod 100 pperiod 1000
compute_msd 100 200 all_groups 0
dump_thermo 1000
run 100000
EOR
  t0=$(date +%s); (cd $d && rm -f msd.out thermo.out && gpumd > gpumd.out 2>&1); echo "T=$T exit=$? $(( $(date +%s)-t0 )) s $(cat $d.log)"
done
echo DONE
