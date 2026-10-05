# -*- coding: utf-8 -*-
"""
Helpers to run the E2 protocols (md/e2_melting, LAMMPS) with GPUMD + NEP89.

    python nep_tools.py coexist <dir> <a> <nx> <nz>      -> <dir>/model.xyz, fcc [001] || z,
                                                           group 0 = lower half (frozen while the upper melts)
    python nep_tools.py bulk <dir> <a> <n>               -> <dir>/model.xyz, n x n x n fcc cubic cells
    python nep_tools.py fromdata <dir>                   -> <dir>/model.xyz from <dir>/ribbon.data (LAMMPS), group 0 = lower half
    python nep_tools.py todump <xyz> <dump> [<data>]     -> GPUMD dump_xyz -> LAMMPS text dump (id x y z)
                                                           (+ LAMMPS data file of frame 0 for the CNA rerun)
    python nep_tools.py thermo <dir>                     -> per-atom averages of <dir>/thermo.out (last run)

GPUMD units: time_step fs, pressure GPa, thermo.out columns T K U Pxx Pyy Pzz Pyz Pxz Pxy then the box
(3 numbers for an orthogonal box, 9 for a triclinic one).
"""
import sys
from pathlib import Path

import numpy as np
from ase.build import bulk
from ase.io import read, write

AL_MASS = 26.982


def coexist(d, a, nx, nz):
    at = bulk("Al", "fcc", a=a, cubic=True).repeat((nx, nx, nz))
    z = at.positions[:, 2]
    at.new_array("group", (z >= at.cell[2, 2] / 2).astype(int))     # 0 = lower half, 1 = upper half
    Path(d).mkdir(parents=True, exist_ok=True)
    write(Path(d) / "model.xyz", at, format="extxyz")
    print(f"{len(at)} atoms, box {at.cell.lengths()}")


def bulk_model(d, a, n):
    at = bulk("Al", "fcc", a=a, cubic=True).repeat((n, n, n))
    Path(d).mkdir(parents=True, exist_ok=True)
    write(Path(d) / "model.xyz", at, format="extxyz")
    print(f"{len(at)} atoms")


def fromdata(d):
    from ase.io.lammpsdata import read_lammps_data
    from ase import Atoms
    src = read_lammps_data(Path(d) / "ribbon.data", atom_style="atomic", units="metal")
    at = Atoms("Al" * len(src), positions=src.positions, cell=src.cell, pbc=True)   # only species, pos, group
    at.wrap()
    at.new_array("group", (at.positions[:, 2] >= at.cell[2, 2] / 2).astype(int))
    write(Path(d) / "model.xyz", at, format="extxyz")
    print(f"{len(at)} atoms, box {at.cell.lengths()}")


def todump(xyz, dump, data=None):
    frames = read(xyz, index=":")
    with open(dump, "w") as fh:
        for k, at in enumerate(frames):
            L = np.diag(at.cell)
            p = at.positions % L
            fh.write(f"ITEM: TIMESTEP\n{k}\nITEM: NUMBER OF ATOMS\n{len(at)}\n"
                     f"ITEM: BOX BOUNDS pp pp pp\n0 {L[0]:.8f}\n0 {L[1]:.8f}\n0 {L[2]:.8f}\n"
                     "ITEM: ATOMS id x y z\n")
            np.savetxt(fh, np.column_stack([np.arange(1, len(at) + 1), p]),
                       fmt=["%d", "%.6f", "%.6f", "%.6f"])
    if data:
        at = frames[0]
        L = np.diag(at.cell)
        with open(data, "w") as fh:
            fh.write(f"frame 0 of {xyz}\n\n{len(at)} atoms\n1 atom types\n\n"
                     f"0 {L[0]:.8f} xlo xhi\n0 {L[1]:.8f} ylo yhi\n0 {L[2]:.8f} zlo zhi\n\n"
                     f"Masses\n\n1 {AL_MASS}\n\nAtoms # atomic\n\n")
            for i, r in enumerate(at.positions % L, 1):
                fh.write(f"{i} 1 {r[0]:.6f} {r[1]:.6f} {r[2]:.6f}\n")
    print(f"{len(frames)} frames -> {dump}")


def thermo(d, skip=0.2):
    t = np.loadtxt(Path(d) / "thermo.out")
    n = len(read(Path(d) / "model.xyz"))
    t = t[int(len(t) * skip):]
    T, K, U = t[:, 0], t[:, 1], t[:, 2]
    box = t[:, 9:12] if t.shape[1] == 12 else t[:, [9, 13, 17]]
    V = box.prod(1) / n
    P = t[:, 3:6].mean(1)
    print(f"T = {T.mean():.2f} K, P = {P.mean():+.4f} GPa, U = {U.mean()/n:.5f} eV/at, "
          f"H = {(K + U).mean()/n:.5f} eV/at, V = {V.mean():.4f} Å3/at")
    return T.mean(), (K + U).mean() / n, V.mean()


if __name__ == "__main__":
    cmd, *a = sys.argv[1:]
    if cmd == "coexist":
        coexist(a[0], float(a[1]), int(a[2]), int(a[3]))
    elif cmd == "bulk":
        bulk_model(a[0], float(a[1]), int(a[2]))
    elif cmd == "fromdata":
        fromdata(a[0])
    elif cmd == "todump":
        todump(*a)
    elif cmd == "thermo":
        thermo(a[0])
