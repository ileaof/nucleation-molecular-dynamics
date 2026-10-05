#!/usr/bin/env bash
# T_m bracketing for NEP89 (GPUMD) with the protocol of md/e2_melting/in.coexist + run_bracket_eam.sh:
# fcc [001] || z, 8 x 8 x 32 cells (8192 atoms), dt = 2 fs, tau_T = 100 steps (0.2 ps), tau_P = 1000 steps (2 ps).
#   1) solid NPT aniso at T0 (5000)  2) upper half melted at 1700 K, lower half frozen (5000)
#   3) liquid at T0, lower half frozen (5000)  4) NPT z only, x-y at the solid lattice (2500)
#   5) production NPT z only, 50000 steps = 100 ps, frames every 1000 steps
# Then the fcc fraction comes from the same LAMMPS cna/atom as E2 (rerun, cutoff 0.854 a0, a0 = 3.9899 Å NEP89 0 K).
#   wsl -d Ubuntu-22.04 -- bash .../md/e2_nep89/run_coexist.sh 800 900 1000 1100
source ~/opt/nucmd_cuda_env.sh
cd "$(dirname "$0")"
POT=$(realpath ../../potentials/NEP89_20250409/nep89_20250409.txt)
for T in "$@"; do
  d=coexist_$T; mkdir -p $d
  python nep_tools.py coexist $d 4.05 8 32 > /dev/null
  cat > $d/run.in <<EOR
potential $POT
velocity $T
time_step 2
ensemble npt_mttk temp $T $T aniso 0 0 tperiod 100 pperiod 1000
dump_thermo 1000
run 5000
fix 0
ensemble nvt_bdp 1700 1700 100
dump_thermo 1000
run 5000
fix 0
ensemble nvt_bdp $T $T 100
dump_thermo 1000
run 5000
ensemble npt_mttk temp $T $T z 0 0 tperiod 100 pperiod 1000
dump_thermo 1000
run 2500
ensemble npt_mttk temp $T $T z 0 0 tperiod 100 pperiod 1000
dump_xyz 1000 prod.xyz
dump_thermo 1000
run 50000
EOR
  (cd $d && rm -f thermo.out prod.xyz && gpumd > gpumd.out 2>&1); echo "T=$T gpumd exit=$?"
  grep "Speed of this run" $d/gpumd.out | tail -1
  python nep_tools.py todump $d/prod.xyz $d/prod.dump $d/frame0.data > /dev/null
  lmp -in in.cna_rerun -var d $d -log $d/cna.log -screen none; echo "T=$T cna exit=$?"
  rm -f $d/prod.dump
done
