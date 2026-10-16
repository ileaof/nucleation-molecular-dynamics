# -*- coding: utf-8 -*-
"""
Static screening of NEP89 for Al-Cu (calorine, relaxed at 0 K):
  fcc Al, fcc Cu lattice and cohesive energies; θ-Al2Cu (C16, I4/mcm) formation energy;
  dilute Cu substitution in fcc Al (energy per Cu relative to fcc Al + fcc Cu).
Reference values are printed by ../../scripts/calphad_alcu_ref.py (COST507) for comparison.

    python screen_alcu.py      (WSL venv ~/opt/venvs/nucmd)
"""
from pathlib import Path

import numpy as np
from ase.build import bulk
from ase.filters import FrechetCellFilter as ExpCellFilter
from ase.optimize import BFGS
from ase.spacegroup import crystal
from calorine.calculators import CPUNEP

POT = str(Path(__file__).resolve().parents[2] / "potentials" / "NEP89_20250409" / "nep89_20250409.txt")


def relax(at, cell=True, fmax=1e-3):
    at.calc = CPUNEP(POT)
    opt = BFGS(ExpCellFilter(at) if cell else at, logfile=None)
    opt.run(fmax=fmax, steps=500)
    return at.get_potential_energy() / len(at), at


eAl, al = relax(bulk("Al", "fcc", a=4.0, cubic=True))
eCu, cu = relax(bulk("Cu", "fcc", a=3.6, cubic=True))
print(f"fcc Al: a = {al.cell[0, 0]:.4f} Å, E = {eAl:.4f} eV/at")
print(f"fcc Cu: a = {cu.cell[0, 0]:.4f} Å, E = {eCu:.4f} eV/at")

th = crystal(["Al", "Cu"], basis=[(0.1581, 0.6581, 0.0), (0.0, 0.0, 0.25)], spacegroup=140,
             cellpar=[6.067, 6.067, 4.877, 90, 90, 90])
eth, th = relax(th)
xcu = sum(1 for s in th.get_chemical_symbols() if s == "Cu") / len(th)
dH = eth - ((1 - xcu) * eAl + xcu * eCu)
c = th.cell.cellpar()
print(f"θ-Al2Cu (C16): {len(th)} át., x_Cu = {xcu:.3f}, a = {c[0]:.3f}, c = {c[2]:.3f} Å, ΔH_f = {dH:+.4f} eV/at = {dH*96.485:+.2f} kJ/mol-át")

sc = bulk("Al", "fcc", a=al.cell[0, 0] / 1.0, cubic=True).repeat(4)       # 256 atoms
e0 = relax(sc.copy(), cell=False)[0] * len(sc)
sc[0].symbol = "Cu"
e1 = relax(sc, cell=True)[0] * len(sc)
Esub = e1 - e0 - eCu + eAl
print(f"Cu substitucional em Al (256 át.): ΔE = {Esub:+.4f} eV por Cu (relativo a Al cfc + Cu cfc) = {Esub*96.485:+.2f} kJ/mol Cu")
