# -*- coding: utf-8 -*-
"""
T_L(Δμ) of Al-Cu (NEP89) from the SGC coexistence runs (runs/coexist_dmu<Δμ>_T<T>/): slope of the CNA fcc fraction
(cna.log, frames every 2 ps) and the mean Cu fraction of the production (mcmd.out, second half);
T_L = zero of df/dt by linear interpolation between the bracketing temperatures.

    python analyze_coexist_sgc.py [dmu]
"""
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
dmu = sys.argv[1] if len(sys.argv) > 1 else "1.013"
pts = []
for d in sorted(HERE.glob(f"runs/coexist_dmu{dmu}_T*"), key=lambda p: float(p.name.split("_T")[1])):
    T = float(d.name.split("_T")[1])
    if not (d / "cna.log").exists():
        continue
    txt = (d / "cna.log").read_text(errors="ignore")
    if "v_ffcc" not in txt:
        continue
    blk = txt[txt.rindex("v_ffcc"):]
    f = np.array([float(x) for x in re.findall(r"^\s+\d+\s+([\d.]+)\s*$", blk, re.M)])
    m = np.loadtxt(d / "mcmd.out", comments="#")
    x = m[len(m) // 2:, 3].mean()
    s = np.polyfit(np.arange(len(f)) * 2.0, f, 1)[0] * 1000
    pts.append((T, s))
    print(f"T = {T:.0f} K: df/dt = {s:+.3f} /ns  (f {f[0]:.3f} → {f[-1]:.3f}), x_Cu médio (sistema) = {x:.4f}")
for (T1, s1), (T2, s2) in zip(pts, pts[1:]):
    if s1 > 0 >= s2:
        print(f"T_L(Δμ = {dmu} eV) = {T1 + (T2 - T1) * s1 / (s1 - s2):.1f} K (colchete {T1:.0f}–{T2:.0f} K)")
