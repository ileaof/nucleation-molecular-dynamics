# -*- coding: utf-8 -*-
"""
Seed insertion for GPUMD seeding (NEP89), avoiding the compressed-liquid problem of a fixed box:
the liquid and the solid are equilibrated separately at T (NPT), then a spherical cavity of radius R + gap is cut
in the liquid and filled with an fcc sphere of radius R built on the solid lattice a(T).

    python seed_tools.py insert <liquid.xyz> <solid.xyz> <R> <out_dir> [shell]
        liquid.xyz / solid.xyz : last frames of the NPT runs (GPUMD dump_xyz, extxyz); a(T) from the solid box
        writes <out_dir>/model.xyz with two grouping methods:
            group 0: 0 = seed, 1 = liquid            (seed frozen while the liquid relaxes around it)
            group 1: 0 = MC shell (|r − c| < R + shell, seed included), 1 = outside
        Cu is kept where it is in the liquid; the seed starts as pure Al (SGC MC sets its composition).
    python seed_tools.py lattice <solid.xyz> <ncell>   -> prints a(T)
"""
import sys
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.io import read, write


def lattice(solid_xyz, n):
    at = read(solid_xyz, index=-1)
    return float(np.mean(at.cell.lengths()) / n)


def insert(liq_xyz, sol_xyz, R, out, shell=12.0, gap=1.0, n_sol=8):
    liq = read(liq_xyz, index=-1)
    a = lattice(sol_xyz, n_sol)
    L = liq.cell.lengths()
    c = L / 2
    d = liq.positions - c
    d -= L * np.round(d / L)
    keep = np.linalg.norm(d, axis=1) > R + gap
    liq = liq[keep]
    m = int(np.ceil(2 * R / a)) + 2
    blk = bulk("Al", "fcc", a=a, cubic=True).repeat((m, m, m))
    p = blk.positions - blk.positions.mean(0)
    p = p[np.linalg.norm(p, axis=1) <= R]
    seed = Atoms("Al" * len(p), positions=p + c, cell=liq.cell, pbc=True)
    at = seed + liq
    at.wrap()
    r = at.positions - c
    r -= L * np.round(r / L)
    dist = np.linalg.norm(r, axis=1)
    g0 = np.r_[np.zeros(len(seed), int), np.ones(len(liq), int)]
    g1 = (dist >= R + shell).astype(int)
    at.new_array("group", np.column_stack([g0, g1]))
    for k in list(at.arrays):
        if k not in ("numbers", "positions", "group"):
            del at.arrays[k]
    Path(out).mkdir(parents=True, exist_ok=True)
    write(Path(out) / "model.xyz", at, format="extxyz")
    print(f"a(T) = {a:.4f} Å; seed {len(seed)} at.; liquid {len(liq)} at. (Cu {sum(at.symbols == 'Cu')}); "
          f"MC shell {int((g1 == 0).sum())} at.; box {L.round(2)}")


if __name__ == "__main__":
    if sys.argv[1] == "insert":
        insert(sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5], float(sys.argv[6]) if len(sys.argv) > 6 else 12.0)
    elif sys.argv[1] == "lattice":
        print(f"{lattice(sys.argv[2], int(sys.argv[3])):.5f}")
