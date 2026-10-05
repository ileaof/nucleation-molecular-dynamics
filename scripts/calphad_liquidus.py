# -*- coding: utf-8 -*-
"""
CALPHAD (pycalphad) for the MD alloy: Al-Nb binary and the nominal Al-Cu-Nb-Fe alloy.

    python scripts/calphad_liquidus.py  ->  results/calphad/

Databases (read-only, the author's OpenCalphad folder):
    ALFENB-B2-2SL.TDB      Al-Fe-Nb assessment (LIQUID, FCC_A1, AL3NB, ...)  -> Al-Nb binary
    COST507-modified.tdb   COST 507 light alloys (has Al-Nb, Al-Cu, Al-Fe in LIQUID/FCC_A1,
                           no Nb intermetallics) -> metastable FCC_A1 liquidus of the alloy

For each composition:
    T_L stable       highest T at which a solid appears, all phases of the database
    T_L(FCC_A1)      same with only LIQUID + FCC_A1 (metastable; what an FCC seed sees)
    k = x_Nb(FCC)/x_Nb(LIQ) at T_L(FCC_A1), and ΔH_m = H_LIQ - H_FCC at T_L(FCC_A1)
"""
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402
from pycalphad import Database, binplot, calculate, equilibrium, variables as v   # noqa: E402
from pycalphad.core.utils import filter_phases                                    # noqa: E402

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OCDIR = Path(r"C:\Users\ileao\OneDrive\Documentos\OpenCalphad\OC6")
OUT = ROOT / "results" / "calphad"
M = {"AL": 26.9815, "CU": 63.546, "NB": 92.90638, "FE": 55.845}


def wt_to_x(wt):
    n = {k: w / M[k] for k, w in wt.items()}
    s = sum(n.values())
    return {k: n[k] / s for k in n}


def phases_for(db, comps):
    return filter_phases(db, [v.Species(c) for c in comps])


def solid_present(db, comps, phases, X, T):
    conds = {v.T: T, v.P: 101325, v.N: 1, **{v.X(k): x for k, x in X.items()}}
    eq = equilibrium(db, comps, phases, conds)
    names = eq.Phase.values.squeeze()
    fr = eq.NP.values.squeeze()
    solids = {str(p): float(f) for p, f in zip(names, fr) if p and p != "LIQUID" and f > 1e-8}
    return solids, eq


def liquidus(db, comps, phases, X, lo=850.0, hi=2200.0, tol=0.05):
    """Bisection on 'a solid phase is present' (true below T_L)."""
    s_hi, _ = solid_present(db, comps, phases, X, hi)
    if s_hi:
        raise ValueError(f"solid present at {hi} K: {s_hi}")
    s_lo, _ = solid_present(db, comps, phases, X, lo)
    if not s_lo:
        raise ValueError(f"no solid at {lo} K")
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        s, _ = solid_present(db, comps, phases, X, mid)
        lo, hi = (mid, hi) if s else (lo, mid)
    first, eq = solid_present(db, comps, phases, X, lo)
    return 0.5 * (lo + hi), first, eq


def fcc_partition_and_latent_heat(db, comps, X, TL):
    """x_Nb in FCC and LIQ just below T_L(FCC); ΔH_m = H_LIQ(x_l) - H_FCC(x_s) [J/mol]."""
    _, eq = solid_present(db, comps, ["LIQUID", "FCC_A1"], X, TL - 0.5)
    ph = list(eq.Phase.values.squeeze())
    xnb = eq.X.sel(component="NB").values.squeeze()
    xl, xs = float(xnb[ph.index("LIQUID")]), float(xnb[ph.index("FCC_A1")])
    out = {}
    for name, x in (("LIQUID", xl), ("FCC_A1", xs)):
        sub = [c for c in comps if c != "VA"]
        # site fractions: LIQUID (AL,..,NB); FCC_A1 (AL,..,NB : VA) — binary Al-Nb only
        pts = np.array([[1 - x, x] + ([1.0] if name == "FCC_A1" else [])])
        r = calculate(db, ["AL", "NB", "VA"], name, T=TL, P=101325, N=1, points=pts, output="HM")
        out[name] = float(r.HM.values.squeeze())
        del sub
    return xl, xs, out["LIQUID"] - out["FCC_A1"]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    lines = ["# CALPHAD — liquidus para a MD\n"]

    # ------------------------------------------------------------ Al-Nb binary
    db = Database(str(OCDIR / "ALFENB-B2-2SL.TDB"))
    comps = ["AL", "NB", "VA"]
    allph = phases_for(db, comps)
    cases = {
        "Al puro": {},
        "Al-5.00wt%Nb (x=0.01506)": {"NB": wt_to_x({"AL": 95, "NB": 5})["NB"]},
        "Al-5.16wt%Nb (x=0.01556)": {"NB": wt_to_x({"AL": 91.9, "NB": 5})["NB"]},
        "x_Nb = 0.03468 (script L723)": {"NB": 0.0346836},
    }
    lines.append("## Binário Al-Nb — ALFENB-B2-2SL.TDB\n")
    lines.append(f"Fases da base para Al-Nb: {', '.join(sorted(allph))}\n")
    lines.append("| composição | T_L estável [K] | 1º sólido | T_L(FCC_A1) metaestável [K] | x_Nb(FCC) | k | ΔH_m [J/mol] | ΔH_m [J/kg] |\n|---|---|---|---|---|---|---|---|")
    for label, X in cases.items():
        if not X:
            TLf, _, _ = liquidus(db, ["AL", "VA"], ["LIQUID", "FCC_A1"], {}, lo=900, hi=1000, tol=0.01)
            lines.append(f"| {label} | {TLf:.2f} | FCC_A1 | {TLf:.2f} | — | — | — | — |")
            continue
        TLs, first, _ = liquidus(db, comps, allph, X)
        TLf, _, _ = liquidus(db, comps, ["LIQUID", "FCC_A1"], X, lo=850, hi=2200)
        xl, xs, dH = fcc_partition_and_latent_heat(db, comps, X, TLf)
        M_mix = (1 - xl) * M["AL"] + xl * M["NB"]
        lines.append(f"| {label} | {TLs:.1f} | {', '.join(first)} | {TLf:.2f} | {xs:.5f} | {xs/xl:.3f} | "
                     f"{dH:.0f} | {dH/M_mix*1e3:.0f} |")

    # phase diagram
    for name, xmax, tmin, tmax in (("AlNb_full", 1.0, 700, 3000), ("AlNb_Al_rich", 0.06, 850, 1900)):
        fig = plt.figure(figsize=(7, 5.5))
        ax = fig.gca()
        binplot(db, comps, allph, {v.X("NB"): (0, xmax, xmax / 120), v.T: (tmin, tmax, 10), v.P: 101325, v.N: 1},
                plot_kwargs={"ax": ax})
        ax.set_title(f"Al-Nb ({'completo' if xmax == 1 else 'canto rico em Al'}) — ALFENB-B2-2SL.TDB")
        if xmax < 1:
            for x in (0.01506, 0.0346836):
                ax.axvline(x, color="k", ls=":", lw=0.8)
        fig.savefig(OUT / f"{name}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    # ---------------------------------------- nominal alloy, metastable FCC liquidus
    db2 = Database(str(OCDIR / "COST507-modified.tdb"))
    comps2 = ["AL", "CU", "FE", "NB", "VA"]
    lines.append("\n## Liga nominal Al-3Cu-5Nb-0,1Fe — COST507-modified.tdb (só LIQUID + FCC_A1)\n")
    lines.append("A COST507 não tem intermetálicos de Nb; o resultado é o liquidus metaestável da FCC_A1.\n")
    lines.append("| composição | T_L(FCC_A1) [K] |\n|---|---|")
    x4 = wt_to_x({"AL": 91.9, "CU": 3.0, "NB": 5.0, "FE": 0.1})
    for label, X, cps in (("Al-3Cu-5Nb-0,1Fe", {k: x4[k] for k in ("CU", "FE", "NB")}, comps2),
                          ("Al-5Nb (binário, COST507)", {"NB": wt_to_x({"AL": 95, "NB": 5})["NB"]}, ["AL", "NB", "VA"]),
                          ("Al-3Cu (binário, COST507)", {"CU": wt_to_x({"AL": 97, "CU": 3})["CU"]}, ["AL", "CU", "VA"])):
        TLf, _, _ = liquidus(db2, cps, ["LIQUID", "FCC_A1"], X, lo=850, hi=2200)
        lines.append(f"| {label} | {TLf:.2f} |")
    lines.append("\nReferência DTA (script de nucleação, L218): T_L = 936,58 K.")
    (OUT / "liquidus.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    sys.exit(main())
