#!/usr/bin/env bash
# T_L(Δμ) of Al-Cu (NEP89) by isothermal solid–liquid coexistence with SGC MCMD at fixed Δμ = μ_Cu − μ_Al.
# Protocol of md/e2_nep89/run_coexist.sh (8192 atoms, fcc [001] || z); MC Al↔Cu swaps (2000 trials every 100 steps)
# from the quench of the liquid on, so both phases hold their equilibrium compositions (β = 1, infinite diffusion).
# Growth (T < T_L) or melting (T > T_L) from the CNA fcc fraction (LAMMPS rerun, as in E2).
#   wsl -d Ubuntu-22.04 -- bash .../md/e4_mcmd/run_coexist_sgc.sh <dmu> <T1> <T2> ...
source ~/opt/nucmd_cuda_env.sh
cd "$(dirname "$0")"
POT=$(realpath ../../potentials/NEP89_20250409/nep89_20250409.txt)
DMU=$1; shift
for T in "$@"; do
  d=runs/coexist_dmu${DMU}_T$T; mkdir -p $d
  python ../e2_nep89/nep_tools.py coexist $d 4.05 8 32 > /dev/null
  MC="mc sgc 100 2000 $T $T 2 Al 0 Cu $DMU"
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
$MC
dump_thermo 1000
run 5000
ensemble npt_mttk temp $T $T z 0 0 tperiod 100 pperiod 1000
$MC
dump_thermo 1000
run 2500
ensemble npt_mttk temp $T $T z 0 0 tperiod 100 pperiod 1000
$MC
dump_xyz 1000 prod.xyz
dump_thermo 1000
run 50000
EOR
  t0=$(date +%s); (cd $d && rm -f thermo.out prod.xyz mcmd.out && gpumd > gpumd.out 2>&1); echo "T=$T gpumd exit=$? $(( $(date +%s)-t0 )) s"
  python ../e2_nep89/nep_tools.py todump $d/prod.xyz $d/prod.dump $d/frame0.data > /dev/null
  lmp -in ../e2_nep89/in.cna_rerun -var d $d -log $d/cna.log -screen none; echo "T=$T cna exit=$?"
  rm -f $d/prod.dump
done
echo DONE
