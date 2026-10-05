#!/usr/bin/env bash
source ~/opt/nucmd_env.sh
cd "$(dirname "$0")"
lmp -in in.al3nb_d022 -log al3nb.log > /dev/null 2>&1; echo exit=$?
grep RESULT al3nb.log | grep -v print
