#!/usr/bin/env bash
# Seeding with NEP89 (GPUMD): pure Al (mode pure) or Al-Cu at fixed Δμ with SGC MC (mode alloy, β = 1).
#   wsl -d Ubuntu-22.04 -- bash .../md/e4_mcmd/run_seed_nep.sh <mode> <dmu|0> <R> <T1> <T2> ...
# Per T: solid 8³ cells NPT (a(T)), liquid 24³ cells NPT (alloy: SGC MC in the liquid stage), seed inserted by
# seed_tools.py, 5 ps with the seed frozen (NVT), then 60 ps NPT iso (alloy: SGC MC in a shell R + 12 Å around the seed).
# Solid atoms by the ten Wolde criterion (LAMMPS rerun, in.q6_rerun, cutoff 0.854 a(T)), frames every 2 ps.
source ~/opt/nucmd_cuda_env.sh
cd "$(dirname "$0")"
POT=$(realpath ../../potentials/NEP89_20250409/nep89_20250409.txt)
MODE=$1; DMU=$2; R=$3; shift 3
for T in "$@"; do
  P=runs/seedprep_${MODE}_T$T; mkdir -p $P/solid $P/liquid
  if [ ! -s $P/solid/final.xyz ]; then
    python ../e2_nep89/nep_tools.py bulk $P/solid 4.05 8 > /dev/null
    printf "potential $POT\nvelocity $T\ntime_step 2\nensemble npt_mttk temp $T $T iso 0 0 tperiod 100 pperiod 1000\ndump_thermo 1000\ndump_xyz 10000 final.xyz\nrun 10000\n" > $P/solid/run.in
    (cd $P/solid && gpumd > gpumd.out 2>&1)
  fi
  if [ ! -s $P/liquid/final.xyz ]; then
    python ../e2_nep89/nep_tools.py bulk $P/liquid 4.05 24 > /dev/null
    MC=""; [ $MODE = alloy ] && MC="mc sgc 100 6000 $T $T 2 Al 0 Cu $DMU"
    printf "potential $POT\nvelocity 1700\ntime_step 2\nensemble npt_mttk temp 1700 1700 iso 0 0 tperiod 100 pperiod 1000\nrun 5000\nensemble npt_mttk temp $T $T iso 0 0 tperiod 100 pperiod 1000\n$MC\ndump_thermo 1000\ndump_xyz 15000 final.xyz\nrun 15000\n" > $P/liquid/run.in
    (cd $P/liquid && gpumd > gpumd.out 2>&1)
  fi
  aT=$(python seed_tools.py lattice $P/solid/final.xyz 8)
  d=runs/seed_${MODE}_R${R}_T$T; mkdir -p $d
  python seed_tools.py insert $P/liquid/final.xyz $P/solid/final.xyz $R $d > $d/insert.log
  MC=""; [ $MODE = alloy ] && MC="mc sgc 100 500 $T $T 2 Al 0 Cu $DMU group 1 0"
  printf "potential $POT\nvelocity $T\ntime_step 2\nfix 0\nensemble nvt_bdp $T $T 100\ndump_thermo 1000\nrun 2500\nensemble npt_mttk temp $T $T iso 0 0 tperiod 100 pperiod 1000\n$MC\ndump_thermo 1000\ndump_xyz 1000 prod.xyz\nrun 30000\n" > $d/run.in
  t0=$(date +%s); (cd $d && rm -f thermo.out prod.xyz mcmd.out && gpumd > gpumd.out 2>&1); echo "$MODE R=$R T=$T exit=$? $(( $(date +%s)-t0 )) s; $(cat $d/insert.log)"
  python ../e2_nep89/nep_tools.py todump $d/prod.xyz $d/prod.dump $d/frame0.data > /dev/null
  lmp -in in.q6_rerun -var d $d -var rc $(python -c "print(0.854*$aT)") -log $d/q6.log -screen none
  rm -f $d/prod.dump
done
echo "DONE $MODE R=$R"
