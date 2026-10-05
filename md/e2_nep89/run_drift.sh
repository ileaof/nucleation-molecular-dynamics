#!/usr/bin/env bash
# NVE energy conservation of NEP89 (GPUMD), 2048 Al atoms started at 2*950 K (as md/e2_melting/in.nve_drift), 20 ps.
#   wsl -d Ubuntu-22.04 -- bash .../md/e2_nep89/run_drift.sh
source ~/opt/nucmd_cuda_env.sh
cd "$(dirname "$0")"
POT=$(realpath ../../potentials/NEP89_20250409/nep89_20250409.txt)
for dt in 2 1; do
  d=drift_dt$dt; mkdir -p $d
  python nep_tools.py bulk $d 4.05 8 > /dev/null
  nsteps=$((20000 / dt)); every=$((1000 / dt))
  cat > $d/run.in <<EOR
potential $POT
velocity 1900
time_step $dt
ensemble nve
dump_thermo $every
run $nsteps
EOR
  (cd $d && gpumd > gpumd.out 2>&1)
  python - $d $dt <<'EOP'
import sys, numpy as np
d, dt = sys.argv[1], float(sys.argv[2])
t = np.loadtxt(f"{d}/thermo.out"); E = t[:, 1] + t[:, 2]; n = len(t); time = np.arange(1, n + 1) * 1.0   # ps
i = n // 4
print(f"NEP89 dt={dt:.0f} fs  T_end={t[-1,0]:.0f} K  Etot drift = {(E[-1]-E[i])/2048/((time[-1]-time[i])/1000):+.4f} eV/atom/ns  (E0={E[0]:.3f})")
EOP
done
