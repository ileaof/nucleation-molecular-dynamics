# -*- coding: utf-8 -*-
"""
Independent recomputation of the continuum table from the script's own equations.

Each column of tablefull_paper_B_MRS_2026.csv used by nucleation_md is recomputed
step by step, following the loop of the continuum script (L870-L1022), and compared
with the printed value (6 significant digits). Nothing is fitted.

    python scripts/verify_continuum.py
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from nucleation_md.sources.continuum import CONTINUUM_PARAMS as P, read_tablefull  # noqa: E402
from nucleation_md.surface import gurtin_murdoch_legacy, gm_surface_stress          # noqa: E402
from nucleation_md.geometry import f_theta, df_legacy                               # noqa: E402

t = read_tablefull()
t = t[np.isclose(t["i"], np.round(t["i"]))].reset_index(drop=True)
i = t["i"].to_numpy().astype(int)

req = P["req"]
dreq = 0.005 * req                                    # L751
DSv0 = -P["DeltaH"] * P["rhos"] / P["TL"]             # L222
gam = gurtin_murdoch_legacy(P["gamma_0"], req, P["sigma_0"], P["lambda_"], P["mu"], P["lambda0"], P["mu0"])

# full loop i = 0..71 (row 0 is not printed)
ii = np.arange(0, 72)
r = req - ii * dreq                                   # L871
g, dg = gam(r), gam.d(r)                              # L879, L883
DT = -2.0 * g / (r * DSv0)                            # value returned by brentq(f_DT) (L890), see note 1
dDT = np.r_[np.nan, np.diff(DT) / np.diff(r)]         # L901
dSvdr = -3 * (DSv0 * DT + dg) / (2 * r * DT) - DSv0 / DT * dDT    # L902
DSv = np.empty_like(r); DSv[0] = DSv0                 # L900
for k in range(1, len(r)):
    DSv[k] = DSv[k - 1] + dSvdr[k] * dreq / DT[k]     # L937
dDSv = np.r_[np.nan, np.diff(DSv) / np.diff(r)]       # L939
th2 = t["theta_2nd"].to_numpy()                       # from brentq(f_het_2nd), L920 (taken from table)
sel = i                                               # printed rows 1..71

checks = {
    "gam_hom (L879)":        (g[sel], t["gam_hom"]),
    "dfgamdr_ana (L883)":    (dg[sel], t["dfgamdr_ana"]),
    "DT (L890)":             (DT[sel], t["DT"]),
    "DSv_hom (L937)":        (DSv[sel], t["DSv_hom"]),
    "dDSv_homdr (L939)":     (dDSv[sel], t["dDSv_homdr"]),
    "dTdr = DT/(8 pi r) (L923)": (DT[sel] / (8 * np.pi * r[sel]), t["dTdr"]),
    "GT_het_2nd = 4 pi r^2 dTdr (L930 hom.)": (DT[sel] * r[sel] / 2, t["GT_het_2nd"]),
    "GB = DSv_hom*DT (L977)":  (DSv[sel] * DT[sel], t["GB"]),
    "GS = gamma (L978)":       (g[sel], t["GS"]),
    "GC = gam/f*dfthetadr (L979)": (g[sel] / f_theta(th2) * df_legacy(th2), t["GC"]),
    "DSs_hom = dgam/DT (L973)":    (dg[sel] / DT[sel], t["DSs_hom"]),
    "DSc (L975)":              (g[sel] / DT[sel] / f_theta(th2) * df_legacy(th2), t["DSc"]),
    "GT_het = -DT*gam*f/b, theta=pi (1st order)": (-DT[sel] * g[sel] / (DSv0 * DT[sel] + dg[sel]), t["GT_het"]),
}

def main():
    print(f"{'column':<45}{'max rel. deviation':>20}")
    for name, (calc, tab) in checks.items():
        tab = np.asarray(tab, dtype=float)
        dev = np.max(np.abs(calc - tab) / np.abs(tab))
        print(f"{name:<45}{dev:>20.2e}")

    # Note 1: what brentq(f_DT, 1e-4, 1) returns
    r1, g1, dg1 = r[1], g[1], dg[1]
    f_DT = lambda x: dg1 / (2 * g1 / (r1 * x) + DSv0) + x
    x0 = -2 * g1 / (r1 * DSv0)
    print(f"\nNote 1 (i=1): f_DT changes sign at DT = {x0:.6f} K "
          f"(f_DT(DT-1e-9) = {f_DT(x0 - 1e-9):.3e}, f_DT(DT+1e-9) = {f_DT(x0 + 1e-9):.3e}); "
          f"zero of f_DT = {-(dg1 + 2 * g1 / r1) / DSv0:.6f} K")
    print(f"Gibbs-Thomson form: DT = 2*Gamma/r with Gamma = -gamma/DSv = {-g1 / DSv0:.6e} m K")


if __name__ == "__main__":
    main()
