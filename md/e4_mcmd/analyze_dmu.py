# -*- coding: utf-8 -*-
"""
x_Cu(Δμ) in bulk solid and liquid from the SGC MCMD calibration runs (runs/dmu_<phase>_T<T>_<dmu>/):
mean Cu fraction over the second half of mcmd.out, mean potential energy per atom over the second half of the
last run in thermo.out (solid and liquid must stay apart: a solid that melts or a liquid that crystallises shows here),
and the partition k = x_s/x_l at equal Δμ (= equal chemical potential: the two-phase equilibrium condition at T_L).

    python analyze_dmu.py [T]      ->  dmu_calibration.md
"""
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
T = sys.argv[1] if len(sys.argv) > 1 else "905"
rows = defaultdict(dict)
for d in sorted((HERE / "runs").glob(f"dmu_*_T{T}_*")):
    ph, dmu = re.match(rf"dmu_(solid|liquid)_T{T}_([\d.]+)", d.name).groups()
    f = d / "mcmd.out"
    if not f.exists():
        continue
    m = np.loadtxt(f, comments="#")
    if m.ndim < 2 or len(m) < 20:
        continue
    h = m[len(m) // 2:]
    xcu, xe = h[:, 3].mean(), h[:, 3].std() / np.sqrt(len(h) / 10)       # ~10-sample correlation
    th = np.loadtxt(d / "thermo.out")
    n = 2048
    U = th[-len(th) // 2:, 2].mean() / n if (d / "thermo.out").exists() else np.nan
    rows[float(dmu)][ph] = (xcu, xe, U)

lines = [f"# Calibração Δμ ↔ x_Cu, SGC MCMD, NEP89, T = {T} K (2048 át., NPT P = 0, 100 ps)\n",
         "| Δμ = μ_Cu − μ_Al [eV] | x_Cu sólido | U sólido [eV/át] | x_Cu líquido | U líquido [eV/át] | k = x_s/x_l |",
         "|---|---|---|---|---|---|"]
for dmu in sorted(rows):
    r = rows[dmu]
    s = r.get("solid", (np.nan,) * 3); l = r.get("liquid", (np.nan,) * 3)
    lines.append(f"| {dmu:.3f} | {s[0]:.4f} ± {s[1]:.4f} | {s[2]:.4f} | {l[0]:.4f} ± {l[1]:.4f} | {l[2]:.4f} | {s[0]/l[0]:.3f} |")
text = "\n".join(lines) + "\n"
(HERE / f"dmu_calibration_T{T}.md").write_text(text, encoding="utf-8")
print(text)
