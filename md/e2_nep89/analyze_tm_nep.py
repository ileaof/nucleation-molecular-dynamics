# -*- coding: utf-8 -*-
"""
T_m of NEP89 from isothermal coexistence (coexist_<T>/): fcc fraction f(t) by CNA (cna.log, frames every 2 ps),
df/dt > 0 below T_m, < 0 above; T_m = zero of a linear fit of df/dt vs T (same as md/e2_melting/analyze_tm.py).
Runs where the solid vanished or the liquid crystallised completely (saturation) are listed but excluded from the fit,
and so are runs farther than WINDOW from the first sign change (growth kinetics is not linear far from T_m:
at 800 K the solid grows slower than at 900 K).

    python analyze_tm_nep.py  ->  tm_bracket_nep89.csv
"""
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
FRAME_PS = 2.0          # dump_xyz every 1000 steps x 2 fs
NPROD = 50              # production thermo rows (last run block)
WINDOW = 50.0           # K around the sign change of df/dt used in the linear fit


def series(d):
    txt = (d / "cna.log").read_text()
    blk = txt[txt.rindex("v_ffcc"):]
    rows = [l.split() for l in blk.splitlines()[1:] if re.match(r"^\s+\d+\s+[\d.]+\s*$", l)]
    f = np.array([float(r[1]) for r in rows])
    T = np.loadtxt(d / "thermo.out")[-NPROD:, 0]
    return f, T


def main():
    pts, out = [], []
    for d in sorted(HERE.glob("coexist_*"), key=lambda p: float(p.name.split("_")[1])):
        T0 = float(d.name.split("_")[1])
        try:
            f, T = series(d)
        except (ValueError, OSError):
            continue
        t = np.arange(len(f)) * FRAME_PS
        slope = np.polyfit(t, f, 1)[0]
        sat = f.min() < 0.02 or f.max() > 0.95
        print(f"T0 = {T0:6.0f} K  <T> = {T.mean():7.1f}  f_fcc {f[0]:.3f} -> {f[-1]:.3f}  "
              f"df/dt = {slope*1e3:+.3f} /ns{'  (saturated: excluded)' if sat else ''}")
        out.append(f"{T0},{T.mean():.2f},{f[0]:.4f},{f[-1]:.4f},{slope:.6e},{int(sat)}")
        if not sat:
            pts.append((T0, slope))
    pts = np.array(pts)
    if len(pts) >= 2 and pts[:, 1].min() < 0 < pts[:, 1].max():
        i = np.where(np.diff(np.sign(pts[:, 1])) < 0)[0][0]
        T_lo, T_hi = pts[i, 0], pts[i + 1, 0]
        T_int = T_lo + (T_hi - T_lo) * pts[i, 1] / (pts[i, 1] - pts[i + 1, 1])
        w = np.abs(pts[:, 0] - T_int) <= WINDOW
        (b, a), cov = np.polyfit(pts[w, 0], pts[w, 1], 1, cov=True) if w.sum() > 3 else (np.polyfit(pts[w, 0], pts[w, 1], 1), None)
        Tm = -a / b
        print(f"sign change between {T_lo:.0f} and {T_hi:.0f} K; linear interpolation T_m = {T_int:.1f} K")
        msg = f"linear fit of {w.sum()} runs within ±{WINDOW:.0f} K ({', '.join(f'{t:.0f}' for t in pts[w, 0])}): T_m = {Tm:.1f} K"
        if cov is not None:
            J = np.array([a / b**2, -1 / b])                     # d(-a/b)/d(b, a)
            msg += f" ± {np.sqrt(J @ cov @ J):.1f} K (fit)"
        print(msg)
    else:
        print("T_m not bracketed")
    (HERE / "tm_bracket_nep89.csv").write_text("T0,T_mean,f0,f1,slope_per_ps,saturated\n" + "\n".join(out) + "\n")


if __name__ == "__main__":
    main()
