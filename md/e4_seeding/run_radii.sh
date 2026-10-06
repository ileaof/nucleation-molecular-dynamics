#!/usr/bin/env bash
# E4: remaining radii, T ranges centred on ΔT = 2Γ/R with Γ ≈ 1.2e-7 m·K (pilot R = 20 Å: T* = 809 K)
D="$(dirname "$0")"
bash "$D/run_scan.sh" 15 740 760 780 800
bash "$D/run_scan.sh" 25 810 830 850 870
bash "$D/run_scan.sh" 30 830 850 870 890
echo "DONE ALL" >> "$D/status.txt"
