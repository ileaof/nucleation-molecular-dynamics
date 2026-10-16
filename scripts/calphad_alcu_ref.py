# -*- coding: utf-8 -*-
"""
CALPHAD reference for the NEP89 Al-Cu screening (COST507-modified.tdb, pycalphad):
  ΔH_f of θ-Al2Cu (ALCU_THETA) at 298.15 K relative to fcc Al + fcc Cu;
  enthalpy of the dilute fcc solution per Cu (∂H/∂x_Cu at x → 0);
  Al-3 wt% Cu: liquidus T_L(FCC_A1), solidus, partition k = x_Cu(FCC)/x_Cu(LIQ) at T_L.

    python scripts/calphad_alcu_ref.py
"""
import warnings
from pathlib import Path

import numpy as np
from pycalphad import Database, calculate, equilibrium, variables as v

warnings.filterwarnings("ignore")
TDB = Path(r"C:\Users\ileao\OneDrive\Documentos\OpenCalphad\OC6\COST507-modified.tdb")
db = Database(str(TDB))
comps = ["AL", "CU", "VA"]
T0 = 298.15


def hm(phase, x, T=T0):
    eq = equilibrium(db, comps, [phase], {v.T: T, v.P: 101325, v.N: 1, v.X("CU"): x}, output="HM")
    return float(eq.HM.values.squeeze())


hAl = float(equilibrium(db, ["AL", "VA"], ["FCC_A1"], {v.T: T0, v.P: 101325, v.N: 1}, output="HM").HM.values.squeeze())
hCu = float(equilibrium(db, ["CU", "VA"], ["FCC_A1"], {v.T: T0, v.P: 101325, v.N: 1}, output="HM").HM.values.squeeze())
x = 1 / 3
hth = hm("ALCU_THETA", x)
print(f"θ-Al2Cu (COST507): ΔH_f(298 K) = {(hth - (1-x)*hAl - x*hCu)/1000:+.2f} kJ/mol-át = {(hth - (1-x)*hAl - x*hCu)/96485:+.4f} eV/at")
for xs in (0.002, 0.005):
    hs = hm("FCC_A1", xs)
    print(f"cfc Al-{xs*100:.1f} at% Cu: ΔH_mix/x_Cu = {(hs - (1-xs)*hAl - xs*hCu)/xs/1000:+.2f} kJ/mol Cu = {(hs - (1-xs)*hAl - xs*hCu)/xs/96485:+.4f} eV por Cu")

wt = 3.0
xcu = (wt / 63.546) / (wt / 63.546 + (100 - wt) / 26.9815)


def phases_at(T):
    eq = equilibrium(db, comps, ["LIQUID", "FCC_A1"], {v.T: T, v.P: 101325, v.N: 1, v.X("CU"): xcu})
    ph = [str(p) for p in eq.Phase.values.squeeze() if p]
    return eq, ph


lo, hi = 850.0, 950.0
while hi - lo > 0.01:
    m = 0.5 * (lo + hi)
    lo, hi = (m, hi) if "FCC_A1" in phases_at(m)[1] else (lo, m)
TL = 0.5 * (lo + hi)
eq, ph = phases_at(TL - 0.05)
X = eq.X.sel(component="CU").values.squeeze()
pl = list(eq.Phase.values.squeeze())
k = float(X[pl.index("FCC_A1")] / X[pl.index("LIQUID")])
lo, hi = 700.0, TL
while hi - lo > 0.01:
    m = 0.5 * (lo + hi)
    lo, hi = (m, hi) if "LIQUID" not in phases_at(m)[1] else (lo, m)
print(f"Al-3 wt% Cu (x_Cu = {xcu:.5f}), só LIQUID + FCC_A1: T_L = {TL:.2f} K, k = {k:.3f}, T_S (FCC) = {0.5*(lo+hi):.2f} K")
