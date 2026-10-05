#!/usr/bin/env bash
# Full height spectra (Morris) for the 9 ribbons, same analysis: Borovikov thin (E2), Borovikov thick (this dir), NEP89.
# Spectra -> ~/runs/spectra/<set>/spectrum_<pot>_<o>.npz ; copied to md/e2_cfm_thick/spectra/<set>/
source ~/opt/nucmd_cuda_env.sh
HERE="$(cd "$(dirname "$0")" && pwd)"; MD="$HERE/.."
export PYTHONPATH="$MD/e2_melting"
R=~/runs/spectra; mkdir -p $R/thin $R/thick $R/nep89
ln -sf "$MD"/e2_melting/cfm_borovikov_{100,110,111}.dump $R/thin/
ln -sf ~/runs/cfm_thick/cfm_borovikov_{100,110,111}.dump $R/thick/
for o in 100 110 111; do
  [ -s $R/nep89/cfm_nep89_$o.dump ] || (cd "$MD/e2_nep89" && python nep_tools.py todump cfm_$o/prod.xyz $R/nep89/cfm_nep89_$o.dump > /dev/null)
done
run() { (cd $R/$1 && python "$MD/e2_melting/cfm_morris.py" $2 $3 40 $4 > morris_$2.txt 2>&1); }
run thin 100 927.9 borovikov & run thin 110 931.9 borovikov & run thin 111 931.0 borovikov & wait
run thick 100 928.0 borovikov & run thick 110 928.1 borovikov & run thick 111 927.9 borovikov & wait
run nep89 100 925.5 nep89 & run nep89 110 926.0 nep89 & run nep89 111 924.2 nep89 & wait
mkdir -p "$HERE/spectra"; for s in thin thick nep89; do mkdir -p "$HERE/spectra/$s"; cp $R/$s/*.npz $R/$s/morris_*.txt "$HERE/spectra/$s/"; done
grep -h "window" "$HERE"/spectra/*/morris_*.txt
