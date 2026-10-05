# -*- coding: utf-8 -*-
"""
CFM stiffness with a white-noise floor in the height spectrum:

    <|A(k)|²> = k_B T / (W L γ̃ k²) + c          (c: height-resolution noise, independent of k)

weighted linear least squares in (1/k², 1) over modes kmin..kmax, then γ₀, ε₁, ε₂ by the cubic expansion
(Hoyt, Asta, Karma 2001) with a Monte Carlo over the stiffness errors — same combine() as E2.
Without the floor (c = 0) the window fit of cfm_morris.py is biased when the amplitudes approach the floor,
which happens for thick ribbons (amplitudes ∝ 1/W).

    python fit_floor.py [kmin] [kmax]          (reads spectra/<set>/spectrum_<pot>_<o>.npz)
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "e2_melting"))
from analyze_cfm import combine, KB   # noqa: E402

SETS = {"thin": "borovikov", "thick": "borovikov", "nep89": "nep89"}
LABEL = {"thin": "Borovikov, fitas da E2", "thick": "Borovikov, fitas espessas (1 ns)", "nep89": "NEP89"}


def fit(d, kmin, kmax, floor=True):
    s = slice(kmin - 1, kmax)
    k, m, se = d["k"][s], d["m"][s], d["se"][s]
    X = np.column_stack([1 / k**2, np.ones_like(k)]) if floor else (1 / k**2)[:, None]
    w = 1 / se
    coef, *_ = np.linalg.lstsq(X * w[:, None], m * w, rcond=None)
    cov = np.linalg.inv((X * w[:, None]).T @ (X * w[:, None]))
    r = (m - X @ coef) / se
    cov *= max(1.0, (r @ r) / max(len(k) - X.shape[1], 1))          # inflate by reduced χ² if > 1
    a, sa = coef[0], np.sqrt(cov[0, 0])
    g = KB * float(d["T"]) / (float(d["W"]) * float(d["L"]) * a) * 1e20
    return g, g * sa / a, (coef[1] if floor else 0.0), (r @ r) / max(len(k) - X.shape[1], 1)


def main(kmin=2, kmax=15):
    rng = np.random.default_rng(1)
    out = [f"# γ₀ com piso de ruído — modos {kmin}–{kmax}\n",
           "| conjunto | orient. | W [Å] | T [K] | γ̃ sem piso [J/m²] | γ̃ com piso [J/m²] | piso c [Å²] | χ²_red |",
           "|---|---|---|---|---|---|---|---|"]
    summary = []
    for st, pot in SETS.items():
        g, e = {}, {}
        for o in ("100", "110", "111"):
            f = HERE / "spectra" / st / f"spectrum_{pot}_{o}.npz"
            if not f.exists():
                continue
            d = np.load(f)
            g0, _, _, _ = fit(d, kmin, kmax, floor=False)
            g[o], e[o], c, chi = fit(d, kmin, kmax, floor=True)
            out.append(f"| {LABEL[st]} | ({o}) | {float(d['W']):.1f} | {float(d['T']):.1f} | {g0:.4f} | "
                       f"{g[o]:.4f} ± {e[o]:.4f} | {c:.4f} | {chi:.2f} |")
        if len(g) == 3:
            gg, e1, e2 = combine(g)
            mc = np.array([combine({o: rng.normal(g[o], e[o]) for o in g}) for _ in range(20000)])
            sd = mc.std(0)
            summary.append(f"| {LABEL[st]} | **{gg:.3f} ± {sd[0]:.3f}** | {e1:.3f} ± {sd[1]:.3f} | {e2:.3f} ± {sd[2]:.3f} |")
    out += ["", "| conjunto | γ₀ [J/m²] | ε₁ | ε₂ |", "|---|---|---|---|"] + summary
    text = "\n".join(out) + "\n"
    (HERE / f"gamma0_floor_k{kmin}-{kmax}.md").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main(*map(int, sys.argv[1:3]))
