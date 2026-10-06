#!/usr/bin/env bash
source ~/opt/nucmd_cuda_env.sh
cd "$(dirname "$0")"
for T in 760 820 880 920; do for ph in liquid; do
  lmp -sf gpu -pk gpu 1 -in in.calib -var T $T -var phase $ph -log calib_${ph}_$T.log > /dev/null 2>&1 || echo "FAIL $ph $T"
done; done
