# -*- coding: utf-8 -*-
"""
Collects the RESUMO lines of runs/<tag>/analysis.md and, per imposed gradient G, finds T*(G): the temperature where the
mean growth rate dN/dt over the replicas changes sign (linear interpolation between the bracketing T0).
For each G: ΔT* = T_m − T*, mean local (∇T)_z at T0 nearest T*, and the two tests of the author:
    (a) ΔT*·r_c/2  vs  A·(∇T)_z           (b) ΔT*·r_c/2 − Γ(G = 0)  vs  A·(∇T)_z
r_c is taken as the mean equivalent radius of the runs at the T0 nearest T* (seed R = 20 Å).

    python summarize_nemd.py [T_m]     ->  nemd_summary.md
"""
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TM = float(sys.argv[1]) if len(sys.argv) > 1 else 936.0
num = r"([+-]?[\d.]+(?:e[+-]?\d+)?)"

rows = []
for f in sorted((HERE / "runs").glob("*/analysis.md")):
    line = next((l for l in f.read_text(encoding="utf-8").splitlines() if l.startswith("RESUMO")), None)
    if not line:
        continue
    R, T0, G, s = re.search(r"R(\d+)_T(\d+)_G([\d.]+)_s(\d+)", line).groups()
    gz = re.search(r"<∇T local> = \(.*?,.*?, " + num + r"±" + num + r"\)", line)
    rows.append(dict(R=int(R), T0=float(T0), G=float(G) * 1e10, s=int(s),
                     dNdt=float(re.search(r"dN/dt = " + num, line).group(1)),
                     dzdt=float(re.search(r"dz_c/dt = " + num, line).group(1)),
                     gz=float(gz.group(1)) * 1e9, gz_e=float(gz.group(2)) * 1e9,
                     gam=float(re.search(r"<ΔT·r/2> = " + num, line).group(1)),
                     Agz=float(re.search(r"<A·\(∇T\)_z> = " + num, line).group(1)),
                     side=float(re.search(r"T_quente − T_fria = " + num, line).group(1)),
                     nout=float(re.search(r"máx\) = " + num, line).group(1)) if "máx" in line else np.nan))

by = defaultdict(lambda: defaultdict(list))
for r in rows:
    by[r["G"]][r["T0"]].append(r)

out = [f"# E4 NEMD — resumo por gradiente imposto (R = 20 Å, T_m = {TM} K)\n",
       "| G [K/m] | T₀ [K] | réplicas | dN/dt por réplica [át/ps] | média | (∇T)_z local [10⁹ K/m] | A·(∇T)_z [10⁻⁷ m·K] | T_q − T_f [K] | sólido fora (máx) |",
       "|---|---|---|---|---|---|---|---|---|"]
star = {}
for G in sorted(by):
    Ts = sorted(by[G])
    means = []
    for T in Ts:
        rr = by[G][T]
        m = np.mean([r["dNdt"] for r in rr]); means.append(m)
        reps = ", ".join("%+.1f" % r["dNdt"] for r in rr)
        out.append(f"| {G:.1e} | {T:.0f} | {len(rr)} | {reps} | {m:+.1f} | "
                   f"{np.mean([r['gz'] for r in rr])/1e9:.2f} | {np.mean([r['Agz'] for r in rr])*1e7:.2f} | "
                   f"{np.mean([r['side'] for r in rr]):+.1f} | {np.nanmax([r['nout'] for r in rr]):.0f} |")
    for i in range(len(Ts) - 1):
        if means[i] > 0 >= means[i + 1]:
            Tst = Ts[i] + (Ts[i + 1] - Ts[i]) * means[i] / (means[i] - means[i + 1])
            near = Ts[i] if abs(Ts[i] - Tst) <= abs(Ts[i + 1] - Tst) else Ts[i + 1]
            rr = by[G][near]
            star[G] = dict(Tst=Tst, lo=Ts[i], hi=Ts[i + 1], gz=np.mean([r["gz"] for r in rr]),
                           Agz=np.mean([r["Agz"] for r in rr]),
                           r=np.mean([2 * r["gam"] / (TM - near) for r in rr]))   # r_eq from ΔT·r/2 at that T0
            break

out += ["", "| G [K/m] | colchete T₀ [K] | T* [K] | ΔT* [K] | r_c [nm] | ΔT*·r_c/2 [10⁻⁷ m·K] | (a) A·(∇T)_z [10⁻⁷ m·K] | "
        "(b) ΔT*·r_c/2 − Γ(0) [10⁻⁷ m·K] |", "|---|---|---|---|---|---|---|---|"]
G0 = None
if 0.0 in star:
    s0 = star[0.0]; G0 = (TM - s0["Tst"]) * s0["r"] / 2
for G in sorted(star):
    s_ = star[G]; dT = TM - s_["Tst"]; gam = dT * s_["r"] / 2
    b = f"{(gam - G0)*1e7:+.2f}" if G0 is not None else "—"
    out.append(f"| {G:.1e} | {s_['lo']:.0f}–{s_['hi']:.0f} | {s_['Tst']:.1f} | {dT:.1f} | {s_['r']*1e9:.3f} | {gam*1e7:.2f} | "
               f"{s_['Agz']*1e7:.2f} | {b} |")
out.append("\nT* sem colchete (todas as T crescem ou todas fundem) não aparece na 2ª tabela.")
text = "\n".join(out) + "\n"
(HERE / "nemd_summary.md").write_text(text, encoding="utf-8")
print(text)
