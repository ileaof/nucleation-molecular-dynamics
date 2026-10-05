#!/usr/bin/env bash
source ~/opt/nucmd_env.sh
cd "$(dirname "$0")"
for p in alnb jelinek; do
  lmp -in in.cna_rerun -var pot $p -log cna_$p.log > /dev/null 2>&1
  echo "$p: $(awk '/v_ffcc/{f=1;next} f&&NF==2{printf "step %s fcc=%.3f; ", $1, $2}' cna_$p.log)"
done
