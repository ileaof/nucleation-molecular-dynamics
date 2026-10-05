#!/usr/bin/env bash
source ~/opt/nucmd_env.sh
cd "$(dirname "$0")"
for s in d022fixed l12; do lmp -in in.al3nb_check -var s $s -log chk_$s.log > /dev/null 2>&1; grep "^RESULT" chk_$s.log; done
