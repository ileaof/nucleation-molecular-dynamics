#!/usr/bin/env bash
D="$(dirname "$0")"
bash "$D/run_seed_nep.sh" pure 0 20 810 790 830
bash "$D/run_seed_nep.sh" alloy 1.013 20 800 760 720 840
echo "QUEUE DONE"
