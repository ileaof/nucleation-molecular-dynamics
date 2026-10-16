# -*- coding: utf-8 -*-
"""
Random Al-Cu fcc block for the liquid diffusion runs (melted in GPUMD before production).
group (grouping method 0): 0 = Al, 1 = Cu  -> compute_msd ... all_groups 0 gives one MSD per species.

    python mk_alcu.py <dir> <n_cells> <wt%Cu> [seed]
"""
import sys
from pathlib import Path

import numpy as np
from ase.build import bulk
from ase.io import write

d, n, wt = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
x = (wt / 63.546) / (wt / 63.546 + (100 - wt) / 26.9815)
at = bulk("Al", "fcc", a=4.05, cubic=True).repeat((n, n, n))
rng = np.random.default_rng(seed)
ncu = int(round(x * len(at)))
idx = rng.choice(len(at), ncu, replace=False)
sym = np.array(at.get_chemical_symbols())
sym[idx] = "Cu"
at.set_chemical_symbols(sym.tolist())
at.new_array("group", (sym == "Cu").astype(int))
Path(d).mkdir(parents=True, exist_ok=True)
write(Path(d) / "model.xyz", at, format="extxyz")
print(f"{len(at)} atoms, {ncu} Cu (x_Cu = {ncu/len(at):.5f}, target {x:.5f})")
