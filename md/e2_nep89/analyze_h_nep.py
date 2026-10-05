# -*- coding: utf-8 -*-
"""
ΔH_m, densities and ΔS_V of NEP89 from H_<phase>_<T>/thermo.out (production run only: dump_thermo is in that block).
H = K + U (P ≈ 0, PV < 1e-5 eV/atom). The final frame of each run is checked by CNA elsewhere (solid stays fcc,
liquid does not crystallise) — see analyze notes in E2_nep89_results.md.

    python analyze_h_nep.py [T_m]
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
N, M_AL, NA, EV = 2048, 26.982e-3, 6.02214076e23, 1.602176634e-19


def stats(d):
    t = np.loadtxt(d / "thermo.out")
    box = t[:, 9:12] if t.shape[1] == 12 else t[:, [9, 13, 17]]
    v = box.prod(1) / N                                   # Å3/atom
    h = (t[:, 1] + t[:, 2]) / N                           # eV/atom
    nb = 10                                               # block averages for the standard error
    hb = h[: len(h) // nb * nb].reshape(nb, -1).mean(1)
    return t[:, 0].mean(), h.mean(), hb.std(ddof=1) / np.sqrt(nb), v.mean()


def main():
    Ts = sorted({d.name.split("_", 2)[2] for d in HERE.glob("H_solid_*")}, key=float)
    Tm = float(sys.argv[1]) if len(sys.argv) > 1 else None
    for T in Ts:
        s, l = HERE / f"H_solid_{T}", HERE / f"H_liquid_{T}"
        if not (s / "thermo.out").exists() or not (l / "thermo.out").exists():
            continue
        Ts_, hs, es, vs = stats(s)
        Tl_, hl, el, vl = stats(l)
        rho = lambda v: M_AL / NA / (v * 1e-30)           # kg/m3
        dH = hl - hs                                      # eV/atom
        dH_kJkg = dH * EV * NA / M_AL / 1e3
        print(f"T = {T} K  (<T> solid {Ts_:.1f}, liquid {Tl_:.1f})")
        print(f"  H_s = {hs:.5f} ± {es:.5f}, H_l = {hl:.5f} ± {el:.5f} eV/atom;  ΔH = {dH*1e3:.2f} meV/atom "
              f"= {dH*EV*NA/1e3:.3f} kJ/mol = {dH_kJkg:.1f} kJ/kg")
        print(f"  ρ_s = {rho(vs):.1f}, ρ_l = {rho(vl):.1f} kg/m3,  ΔV/V_s = {(vl/vs-1)*100:.2f} %")
        if Tm:
            dSv = -dH_kJkg * 1e3 * rho(vs) / Tm
            print(f"  ΔS_V = -ΔH·ρ_s/T_m = {dSv:.4e} J m-3 K-1  (T_m = {Tm} K)")


if __name__ == "__main__":
    main()
