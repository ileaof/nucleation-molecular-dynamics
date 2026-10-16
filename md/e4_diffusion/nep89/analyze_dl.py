# -*- coding: utf-8 -*-
"""
Liquid diffusion coefficients of Cu and Al in Al-3 wt% Cu (NEP89) from dl_T<T>/msd.out (compute_msd all_groups 0:
group 0 = Al, group 1 = Cu): D = slope of MSD(t) = msdx + msdy + msdz over 10–40 ps / 6. Arrhenius fit D = D0·exp(−Q/RT)
and comparison with the Arrhenius of the SDAS script (Cu in liquid Al: D0L = 1.06e-7 m²/s, QL = 24000 J/mol,
scriptsdasmodelAl3Cu5Nb01Fe_zoom_horizontal_B_MRS_Meeting_2026.py), then the SDAS of Al-3Cu-5Nb-0,1Fe recomputed with
the MD D_L of Cu (everything else unchanged).

    python analyze_dl.py      (Windows Python, repo root on sys.path)  ->  dl_nep89.md, dl_arrhenius.png
"""
import copy
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
R = 8.31451
D0_SCRIPT, Q_SCRIPT = 1.06e-7, 24000.0


def read_msd(f):
    d = np.loadtxt(f, comments="#")
    t = d[:, 0]
    out = {}
    for g, name in ((0, "Al"), (1, "Cu")):
        c = 1 + 6 * g
        out[name] = d[:, c] + d[:, c + 1] + d[:, c + 2]
    return t, out


rows = []
for d in sorted((p for p in HERE.glob("dl_T*") if p.is_dir()), key=lambda p: float(p.name[4:])):
    if not d.is_dir() or not (d / "msd.out").exists():
        continue
    T = float(d.name[4:])
    t, m = read_msd(d / "msd.out")
    w = (t >= 10) & (t <= 40)
    D = {k: np.polyfit(t[w], v[w], 1)[0] / 6 * 1e-8 for k, v in m.items()}     # Å²/ps → m²/s
    # uncertainty: slopes over 10–25 and 25–40 ps
    w1, w2 = (t >= 10) & (t < 25), (t >= 25) & (t <= 40)
    Derr = {k: abs(np.polyfit(t[w1], v[w1], 1)[0] - np.polyfit(t[w2], v[w2], 1)[0]) / 6 * 1e-8 / 2 for k, v in m.items()}
    rows.append((T, D["Cu"], Derr["Cu"], D["Al"], Derr["Al"]))
rows = np.array(rows)
lines = ["# Difusão no líquido Al-3 wt% Cu, NEP89 (GPUMD, 16 384 át., 212 Cu, NPT P = 0, 200 ps)\n",
         "| T [K] | D_Cu MD [m²/s] | D_Cu script [m²/s] | MD/script | D_Al MD [m²/s] |", "|---|---|---|---|---|"]
for T, dc, ec, da, ea in rows:
    ds = D0_SCRIPT * np.exp(-Q_SCRIPT / (R * T))
    lines.append(f"| {T:.2f} | {dc:.3e} ± {ec:.1e} | {ds:.3e} | {dc/ds:.2f} | {da:.3e} ± {ea:.1e} |")
fit = {}
for k, col in (("Cu", 1), ("Al", 3)):
    b, a = np.polyfit(1 / rows[:, 0], np.log(rows[:, col]), 1)
    fit[k] = (np.exp(a), -b * R)
lines += ["", f"Arrhenius MD: Cu D0 = {fit['Cu'][0]:.3e} m²/s, Q = {fit['Cu'][1]:.0f} J/mol; "
          f"Al D0 = {fit['Al'][0]:.3e} m²/s, Q = {fit['Al'][1]:.0f} J/mol",
          f"Script (Cu): D0 = {D0_SCRIPT:.3e} m²/s, Q = {Q_SCRIPT:.0f} J/mol"]

# SDAS of Al-3Cu-5Nb-0,1Fe with the MD D_L of Cu
import pandas as pd                                                          # noqa: E402
from nucleation_md import FerreiraModel                                      # noqa: E402
from nucleation_md.sdas import ALLOY_SDAS, sdas_from_gamma, power_law_fit    # noqa: E402
from nucleation_md.sources.continuum import ContinuumSource                  # noqa: E402
alloy = "Al3Cu5Nb01Fe"
sy0 = ALLOY_SDAS[alloy]
sy1 = copy.deepcopy(sy0)
cu = [s for s in sy1.solutes if s.name == "Cu"][0]
cu.D0L, cu.QL = fit["Cu"]
src = ContinuumSource(alloy=alloy, dSv_mode="closure", theta_mode="homogeneous")
st = src.state_at_gradient(src.grad_T_exp)
g1, g2 = FerreiraModel().gamma_tensor(st, 1), FerreiraModel().gamma_tensor(st, 2)
exp = pd.read_csv(ROOT / "data" / "experimental" / f"sdas_{alloy}.csv", comment="#")
lines += ["", f"## SDAS Al-3Cu-5Nb-0,1Fe: D_L do Cu do script × da MD (T_L = {sy0.TL} K: "
          f"script {sy0.DL()[0]:.3e}, MD {sy1.DL()[0]:.3e} m²/s)\n",
          "| curva | D_L do Cu | SDAS [μm] | dentro das barras | bias_log |", "|---|---|---|---|---|"]
for lab, G, model in (("MRB Γ¹ˢᵗ contínuo", g1, "MRB"), ("MRB Γ²ⁿᵈ contínuo", g2, "MRB"), ("RB Γ²ⁿᵈ contínuo", g2, "RB")):
    for dl, sy in (("script", sy0), ("MD", sy1)):
        t, S = sdas_from_gamma(G, sy, model)
        a, b = power_law_fit(t, S)
        p = a * exp.tSL_s ** b
        inside = ((p <= exp.SDAS_um + exp.err_plus_um) & (p >= exp.SDAS_um - exp.err_minus_um)).mean()
        lines.append(f"| {lab} | {dl} | {a:.2f}·t^{b:.3f} | {inside:.0%} | {np.mean(np.log(p/exp.SDAS_um)):+.3f} |")
text = "\n".join(lines) + "\n"
(HERE / "dl_nep89.md").write_text(text, encoding="utf-8")
print(text)

import matplotlib                                                            # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                              # noqa: E402
fig, ax = plt.subplots(figsize=(6, 4.2))
Ti = np.linspace(850, 1250, 50)
ax.semilogy(1e3 / rows[:, 0], rows[:, 1], "o", color="#2a78d6", ms=8, label="Cu, MD (NEP89)")
ax.semilogy(1e3 / rows[:, 0], rows[:, 3], "s", color="#1baf7a", ms=7, label="Al, MD (NEP89)")
ax.semilogy(1e3 / Ti, D0_SCRIPT * np.exp(-Q_SCRIPT / (R * Ti)), "k-", lw=1.8, label="Cu, script (D0 = 1,06e-7, Q = 24 kJ/mol)")
ax.semilogy(1e3 / Ti, fit["Cu"][0] * np.exp(-fit["Cu"][1] / (R * Ti)), color="#2a78d6", lw=1.5, ls="--", label="Cu, Arrhenius MD")
ax.set_xlabel("1000/T [1/K]"); ax.set_ylabel("D [m²/s]"); ax.legend(fontsize=8, frameon=False)
ax.set_title("Difusão no líquido Al-3 wt% Cu", loc="left")
fig.tight_layout(); fig.savefig(HERE / "dl_arrhenius.png", dpi=140)
