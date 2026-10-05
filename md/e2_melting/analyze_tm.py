# -*- coding: utf-8 -*-
"""
T_m from isothermal coexistence (bracket_<pot>_<T>.log):
fcc fraction f(t) from CNA; its slope df/dt > 0 below T_m (solid grows), < 0 above.
T_m = zero of a linear fit of df/dt versus T.

    python analyze_tm.py  ->  tm_bracket.csv, prints T_m per potential
"""
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DT_PS = 0.002


def production(log):
    txt = log.read_text()
    blk = txt[txt.rindex("v_fsol"):]
    blk = blk[:blk.index("Loop time")] if "Loop time" in blk else blk
    rows = [l.split() for l in blk.splitlines()[1:] if re.match(r"^\s+\d+\s", l)]
    return np.array([[float(x) for x in r] for r in rows])   # step temp press pzz pe enthalpy lz fsol


def main():
    out = []
    for pot in ("alnb", "jelinek", "borovikov", "mendelev"):
        pts = []
        for log in sorted(HERE.glob(f"bracket_{pot}_*.log")):
            T = float(log.stem.split("_")[-1])
            try:
                d = production(log)
            except ValueError:
                continue
            if len(d) < 5:
                continue
            t = (d[:, 0] - d[0, 0]) * DT_PS
            slope = np.polyfit(t, d[:, 7], 1)[0]                     # 1/ps
            pts.append((T, d[:, 1].mean(), d[0, 7], d[-1, 7], slope))
            out.append(f"{pot},{T},{d[:,1].mean():.2f},{d[0,7]:.4f},{d[-1,7]:.4f},{slope:.6e}")
        if not pts:
            continue
        pts = np.array(pts)
        print(f"\n{pot}")
        for T, Tm_run, f0, f1, s in pts:
            print(f"  T = {T:6.0f} K  <T> = {Tm_run:7.1f}  f_fcc {f0:.3f} -> {f1:.3f}  df/dt = {s*1e3:+.3f} /ns")
        if len(pts) >= 2 and pts[:, 4].min() < 0 < pts[:, 4].max():
            b, a = np.polyfit(pts[:, 0], pts[:, 4], 1)
            print(f"  T_m (df/dt = 0, linear fit) = {-a/b:.1f} K")
        else:
            print("  T_m not bracketed by these temperatures")
    (HERE / "tm_bracket.csv").write_text("pot,T,T_mean,f0,f1,slope_per_ps\n" + "\n".join(out) + "\n")


if __name__ == "__main__":
    main()
