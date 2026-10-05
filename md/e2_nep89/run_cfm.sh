#!/usr/bin/env bash
# Capillary fluctuation method, NEP89 (GPUMD), protocol of md/e2_melting/in.cfm (same ribbons, built by in.build_cfm):
#   NPT aniso at T_m (10000) -> upper half melted at 1500 K, lower half frozen (10000) -> liquid at T_m, frozen (10000)
#   -> NPT z only at T_m (25000) -> production NPH z only, 200000 steps = 400 ps, frames every 2000 steps.
# Heights then from the Morris order parameter (md/e2_melting/cfm_morris.py, pot = nep89).
#   wsl -d Ubuntu-22.04 -- bash .../md/e2_nep89/run_cfm.sh <T_m> 100 110 111
source ~/opt/nucmd_cuda_env.sh
cd "$(dirname "$0")"
POT=$(realpath ../../potentials/NEP89_20250409/nep89_20250409.txt)
T=$1; shift
for o in "$@"; do
  d=cfm_$o; rm -f $d/thermo.out $d/prod.xyz
  cat > $d/run.in <<EOR
potential $POT
velocity $T
time_step 2
ensemble npt_mttk temp $T $T aniso 0 0 tperiod 100 pperiod 1000
dump_thermo 1000
run 10000
fix 0
ensemble nvt_bdp 1500 1500 100
dump_thermo 1000
run 10000
fix 0
ensemble nvt_bdp $T $T 100
dump_thermo 1000
run 10000
ensemble npt_mttk temp $T $T z 0 0 tperiod 100 pperiod 1000
dump_thermo 1000
run 25000
ensemble nph_mttk z 0 0 pperiod 1000
dump_thermo 1000
dump_xyz 2000 prod.xyz
run 200000
EOR
  (cd $d && gpumd > gpumd.out 2>&1); echo "cfm $o exit=$?"
  grep "Speed of this run" $d/gpumd.out | tail -1
  python nep_tools.py todump $d/prod.xyz $d/cfm_nep89_$o.dump > /dev/null
done
