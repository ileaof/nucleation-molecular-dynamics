#!/usr/bin/env bash
D="$(dirname "$0")"
bash "$D/run_seed_nep.sh" alloy 1.013 20 840 860 880
bash "$D/run_seed_nep.sh" pure 0 20 850 870 890
echo "QUEUE2 DONE"
