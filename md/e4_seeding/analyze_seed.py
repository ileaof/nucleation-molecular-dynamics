# -*- coding: utf-8 -*-
"""
E4 seeding analysis: growth or melting of a spherical seed at T, and the critical condition r_c(T).

For each logs/seed_R<R>_T<T>_s<s>.log (production stage, column v_Nsol):
    N_cr(t) = N_sol(t) / f_s(T)        f_s = fraction of a perfect crystal labelled solid by the criterion
                                       (calibration in_calib, q6·q6 > 0.5 and >= 7 bonds; liquid false positives ~0)
    dN/dt from a linear fit over the production;  sign > 0 grows, < 0 melts
    r = (3 N_cr v_at / 4π)^(1/3), v_at = a(T)³/4, with N_cr averaged over the first 10 ps
For each R, T* = temperature where dN/dt changes sign (linear interpolation between the bracketing T);
r_c = r at T*, ΔT = T_m − T*, Γ_GT = ΔT·r_c/2 (∇T = 0).

    python analyze_seed.py [T_m]      ->  seeding_results.md, seeding_Nt.png
"""
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TM = float(sys.argv[1]) if len(sys.argv) > 1 else 936.0       # Borovikov 2024, md/e2_melting/E2_results.md
# calibration (in.calib, 4000 atoms, 20 ps): fraction of the bulk solid labelled solid, q6·q6 > 0.5, >= 7 bonds
CAL_T = np.array([760.0, 820.0, 880.0, 920.0])
CAL_F = np.array([0.9285, 0.8948, 0.8728, 0.8472])


def f_s(T):
    return float(np.interp(T, CAL_T, CAL_F))


def read(path):
    txt = path.read_text(errors="ignore")
    m = re.search(r"lattice at T: a = ([\d.]+)", txt)
    a = float(m.group(1))
    blocks = txt.split("v_Nsol")
    prod = blocks[-1]                                           # last thermo block = production
    rows = [l.split() for l in prod.splitlines() if re.match(r"^\s+\d+\s+[\d.]", l)]
    rows = [r for r in rows if len(r) == 6]
    d = np.array(rows, float)
    t = (d[:, 0] - d[0, 0]) * 0.002                              # ps
    return a, t, d[:, 5], "Total wall time" in txt


def main():
    data = defaultdict(dict)
    for p in sorted((HERE / "logs").glob("seed_R*_T*_s*.log")):
        R, T, s = map(int, re.findall(r"R(\d+)_T(\d+)_s(\d+)", p.name)[0])
        a, t, Ns, done = read(p)
        if len(t) < 5:
            continue
        Ncr = Ns / f_s(T)
        slope = np.polyfit(t, Ncr, 1)[0]                          # atoms/ps
        N0 = Ncr[t <= 10].mean()
        r0 = (3 * N0 * a**3 / 4 / (4 * np.pi)) ** (1 / 3) / 10    # nm
        data[R].setdefault(T, []).append(dict(s=s, t=t, Ncr=Ncr, slope=slope, N0=N0, r0=r0, a=a, done=done))

    lines = [f"# E4 — seeding a ∇T = 0 (Borovikov 2024, T_m = {TM:.1f} K)\n",
             "| R [Å] | T [K] | réplica | N_cr inicial | r inicial [nm] | dN_cr/dt [át/ps] | resultado | completa |",
             "|---|---|---|---|---|---|---|---|"]
    crit = []
    for R in sorted(data):
        Ts = sorted(data[R])
        for T in Ts:
            for d in data[R][T]:
                lines.append(f"| {R} | {T} | {d['s']} | {d['N0']:.0f} | {d['r0']:.3f} | {d['slope']:+.1f} | "
                             f"{'cresce' if d['slope'] > 0 else 'funde'} | {'sim' if d['done'] else 'não'} |")
        sl = np.array([np.mean([d["slope"] for d in data[R][T]]) for T in Ts])
        for i in range(len(Ts) - 1):
            if sl[i] > 0 >= sl[i + 1]:
                Tstar = Ts[i] + (Ts[i + 1] - Ts[i]) * sl[i] / (sl[i] - sl[i + 1])
                r = np.interp(Tstar, Ts, [np.mean([d["r0"] for d in data[R][T]]) for T in Ts])
                dT = TM - Tstar
                crit.append((R, Ts[i], Ts[i + 1], Tstar, r, dT, dT * r * 1e-9 / 2))
    lines += ["", "| R [Å] | T que cresce | T que funde | T* [K] | r_c [nm] | ΔT = T_m − T* [K] | Γ = ΔT·r_c/2 [m·K] |",
              "|---|---|---|---|---|---|---|"]
    for R, Tg, Tm_, Ts_, r, dT, G in crit:
        lines.append(f"| {R} | {Tg} | {Tm_} | {Ts_:.1f} | {r:.3f} | {dT:.1f} | {G:.3e} |")
    if crit:
        G = np.array([c[6] for c in crit])
        lines += ["", f"Γ médio = {G.mean():.3e} ± {G.std(ddof=1) if len(G) > 1 else 0:.1e} m·K (desvio entre raios; "
                  "1 réplica por T, colchete de 20 K → T* ± ~5 K)"]
    text = "\n".join(lines) + "\n"
    (HERE / "seeding_results.md").write_text(text, encoding="utf-8")
    print(text)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(data), figsize=(5 * len(data), 4), squeeze=False)
    for ax, R in zip(axes[0], sorted(data)):
        for T in sorted(data[R]):
            for d in data[R][T]:
                ax.plot(d["t"], d["Ncr"], label=f"{T} K" if d["s"] == data[R][T][0]["s"] else None)
        ax.set_title(f"R = {R} Å"); ax.set_xlabel("t [ps]"); ax.set_ylabel("N cristalino"); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(HERE / "seeding_Nt.png", dpi=130)
    if crit:
        fig, ax = plt.subplots(figsize=(5.5, 4))
        r = np.array([c[4] for c in crit]); dT = np.array([c[5] for c in crit]); G = np.array([c[6] for c in crit])
        rr = np.linspace(1.2, 3.2, 100)
        ax.plot(rr, 2 * G.mean() / (rr * 1e-9), "k--", label=f"ΔT = 2Γ/r_c, Γ = {G.mean():.2e} m·K")
        ax.plot(r, dT, "o", ms=8, label="seeding (Borovikov 2024)")
        ax.set_xlabel("r_c [nm]"); ax.set_ylabel("ΔT = T_m − T* [K]"); ax.legend(); ax.set_title("E4 — Gibbs–Thomson a ∇T = 0")
        fig.tight_layout(); fig.savefig(HERE / "seeding_DT_rc.png", dpi=130)


if __name__ == "__main__":
    main()
