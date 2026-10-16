# -*- coding: utf-8 -*-
"""
NEP89 seeding (runs/seed_<mode>_R<R>_T<T>/q6.log): N_sol(t) by the ten Wolde criterion, every 2 ps.
f_s(T) = fraction of the bulk solid labelled solid, from the last frame of runs/seedprep_<mode>_T<T>/solid (computed
here with the same LAMMPS rerun); N_cr = N_sol / f_s; r = (3 N_cr a³/4 / 4π)^(1/3) from the first 10 ps.
T* by the sign change of dN_cr/dt; ΔT = T_ref − T* with T_ref = T_m (pure, 919 K, E2 NEP89) or T_L(Δμ) (alloy,
from the SGC coexistence); Γ = ΔT·r_c/2.

    python analyze_seed_nep.py [T_L_alloy]     (WSL venv; calls lmp for f_s)  ->  seeding_nep.md
"""
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TREF = {"pure": 919.0, "alloy": float(sys.argv[1]) if len(sys.argv) > 1 else np.nan}


def series(log):
    rows = re.findall(r"^\s+(\d+)\s+(\d+)\s*$", log.read_text(errors="ignore"), re.M)
    return np.array([[float(a), float(b)] for a, b in rows])


def lattice(xyz, n=8):
    from ase.io import read
    return float(np.mean(read(xyz, index=-1).cell.lengths()) / n)


def f_s(prep):
    """Fraction of the bulk solid labelled solid (last frame of the solid prep run)."""
    sol = prep / "solid"
    a = lattice(sol / "final.xyz")
    if not (sol / "q6.log").exists():
        subprocess.run(["python", str(HERE.parent / "e2_nep89" / "nep_tools.py"), "todump", str(sol / "final.xyz"),
                        str(sol / "prod.dump"), str(sol / "frame0.data")], check=True, capture_output=True)
        subprocess.run(["lmp", "-in", str(HERE / "in.q6_rerun"), "-var", "d", str(sol), "-var", "rc", str(0.854 * a),
                        "-log", str(sol / "q6.log"), "-screen", "none"], check=True)
        (sol / "prod.dump").unlink()
    s = series(sol / "q6.log")
    return s[-1, 1] / 2048.0, a


def main():
    data = defaultdict(dict)
    for d in sorted(HERE.glob("runs/seed_*_R*_T*")):
        mode, R, T = re.match(r"seed_(\w+)_R(\d+)_T(\d+)", d.name).groups()
        if not (d / "q6.log").exists() or "Total wall time" not in (d / "q6.log").read_text(errors="ignore"):
            continue
        fs, a = f_s(HERE / "runs" / f"seedprep_{mode}_T{T}")
        s = series(d / "q6.log")
        t, N = s[:, 0] * 2.0, s[:, 1] / fs
        slope = np.polyfit(t, N, 1)[0]
        r0 = (3 * N[t <= 10].mean() * a**3 / 4 / (4 * np.pi)) ** (1 / 3) / 10
        data[(mode, int(R))][int(T)] = dict(slope=slope, r0=r0, fs=fs, N0=N[0], N1=N[-1])
    out = ["# Seeding NEP89 (GPUMD): Al puro × Al-Cu com MC SGC (β = 1)\n",
           "| modo | R [Å] | T [K] | f_s | N_cr inicial → final | dN/dt [át/ps] | r inicial [nm] |", "|---|---|---|---|---|---|---|"]
    res = []
    for (mode, R), runs in sorted(data.items()):
        Ts = sorted(runs)
        for T in Ts:
            r = runs[T]
            out.append(f"| {mode} | {R} | {T} | {r['fs']:.3f} | {r['N0']:.0f} → {r['N1']:.0f} | {r['slope']:+.1f} | {r['r0']:.3f} |")
        for i in range(len(Ts) - 1):
            a, b = runs[Ts[i]]["slope"], runs[Ts[i + 1]]["slope"]
            if a > 0 >= b:
                Tst = Ts[i] + (Ts[i + 1] - Ts[i]) * a / (a - b)
                rc = np.interp(Tst, Ts, [runs[T]["r0"] for T in Ts])
                dT = TREF[mode] - Tst
                res.append((mode, R, Ts[i], Ts[i + 1], Tst, rc, dT, dT * rc * 1e-9 / 2))
    out += ["", "| modo | R [Å] | colchete [K] | T* [K] | r_c [nm] | T_ref [K] | ΔT [K] | Γ = ΔT·r_c/2 [m·K] |", "|---|---|---|---|---|---|---|---|"]
    for mode, R, lo, hi, Tst, rc, dT, G in res:
        out.append(f"| {mode} | {R} | {lo}–{hi} | {Tst:.1f} | {rc:.3f} | {TREF[mode]:.1f} | {dT:.1f} | {G:.3e} |")
    text = "\n".join(out) + "\n"
    (HERE / "seeding_nep.md").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
