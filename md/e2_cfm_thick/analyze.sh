#!/usr/bin/env bash
# γ₀ from the thick ribbons: Morris heights (../e2_melting/cfm_morris.py) on the dumps in ~/runs/cfm_thick.
#   wsl -d Ubuntu-22.04 -- bash .../md/e2_cfm_thick/analyze.sh T100 T110 T111
source ~/opt/nucmd_cuda_env.sh
HERE="$(cd "$(dirname "$0")" && pwd)"
cd ~/runs/cfm_thick
export PYTHONPATH="$HERE/../e2_melting"
i=0; for T in "$@"; do o=$(echo 100 110 111 | cut -d' ' -f$((i+1))); i=$((i+1))
  python "$HERE/../e2_melting/cfm_morris.py" $o $T 40 borovikov > "$HERE/morris_$o.txt" 2>&1 &
done; wait
cat "$HERE"/morris_*.txt
