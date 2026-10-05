#!/usr/bin/env bash
source ~/opt/nucmd_env.sh
cd "$(dirname "$0")"
lmp -in in.al3nb_d022_bor -log al3nb_bor.log > /dev/null 2>&1; echo exit=$?; grep "^RESULT" al3nb_bor.log
for s in d022fixed l12; do lmp -in in.al3nb_check_bor -var s $s -log chk_bor_$s.log > /dev/null 2>&1; grep "^RESULT" chk_bor_$s.log; done
