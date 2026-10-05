# -*- coding: utf-8 -*-
"""
CALPHAD (pycalphad, COST507-modified.tdb) for the two alloys of the continuum scripts:
liquidus, solidus and enthalpy of fusion, as a reference for the solute effect that the pure-Al MD lacks.

    python scripts/calphad_alloys.py  ->  results/calphad/alloys_TL_dH.md, alloys_H_T.png, alloys_path_<liga>.csv

Database (read-only): C:\\Users\\ileao\\OneDrive\\Documentos\\OpenCalphad\\OC6\\COST507-modified.tdb

For each composition (wt.%), equilibrium with all phases of the database for its elements:
    T_L          highest T with a solid; first solid phase
    T_L(FCC_A1)  T at which FCC_A1 (α-Al) appears
    T_S          T at which the liquid disappears (equilibrium solidification)
Enthalpy of fusion, reported in three definitions (the choice is the author's):
    ΔH_a = H_LIQ(x0) − H_FCC(x0) at T_L(FCC_A1): single liquid vs single FCC of the nominal composition
           (the analogue of ΔH_m of a pure metal; what the MD of pure Al measures)
    ΔH_b = H_eq(T_L) − H_eq(T_S): total heat released along equilibrium solidification (latent + sensible)
    ΔH_c = ΔH_b − ∫_{T_S}^{T_L} c_p,solid dT: latent part only, c_p of the equilibrium solid assemblage
    ΔH_b, ΔH_c are taken over the α-Al interval [T_S, T_L(FCC_A1)] (for Al-Cu-Nb-Fe Al3Nb forms far above it);
    ΔH_a is also given with the composition of the liquid at T_L(FCC_A1) (= x0 when FCC_A1 is the first solid).
Phase excluded: ALTI (ordered L1_0 γ-TiAl; its Al-rich end is FCC Al + 2 J/mol but it carries
L(Al,Nb) = −80800 + 30T, so without Ti it swallows the Nb and the Al3Nb — artifact for these alloys).
"""
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402
from pycalphad import Database, equilibrium, variables as v   # noqa: E402
from pycalphad.core.utils import filter_phases                # noqa: E402

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
TDB = Path(r"C:\Users\ileao\OneDrive\Documentos\OpenCalphad\OC6\COST507-modified.tdb")
OUT = ROOT / "results" / "calphad"
M = {"AL": 26.9815, "SI": 28.0855, "MG": 24.305, "FE": 55.845, "CU": 63.546, "NB": 92.90638}
ALLOYS = {
    "Al puro": {"AL": 100.0},
    "Al-0,8Si-0,6Mg-0,2Fe": {"AL": 98.4, "SI": 0.8, "MG": 0.6, "FE": 0.2},
    "Al-3Cu-5Nb-0,1Fe": {"AL": 91.9, "CU": 3.0, "NB": 5.0, "FE": 0.1},
}
EXCLUDE = {"ALTI"}
SCRIPT = {"Al-0,8Si-0,6Mg-0,2Fe": (925.32, 335300.0), "Al-3Cu-5Nb-0,1Fe": (936.58, 335300.0)}


def wt_to_x(wt):
    n = {k: w / M[k] for k, w in wt.items()}
    s = sum(n.values())
    return {k: n[k] / s for k in n}


class Alloy:
    def __init__(self, db, wt):
        self.db, self.wt = db, wt
        self.x = wt_to_x(wt)
        self.comps = sorted(wt) + ["VA"]
        self.phases = sorted(set(filter_phases(db, [v.Species(c) for c in self.comps])) - EXCLUDE)
        self.Mmol = sum(self.x[k] * M[k] for k in self.x) * 1e-3          # kg/mol of atoms

    def eq(self, T, phases=None, x=None):
        x = x or self.x
        conds = {v.T: T, v.P: 101325, v.N: 1, **{v.X(k): x[k] for k in x if k != "AL"}}
        return equilibrium(self.db, self.comps, phases or self.phases, conds, output="HM")

    def state(self, T):
        """{phase: mole fraction} and HM [J/mol] at a single T."""
        e = self.eq(float(T))
        ph = e.Phase.values.squeeze()
        f = e.NP.values.squeeze()
        out = {}
        for p, n in zip(np.atleast_1d(ph), np.atleast_1d(f)):
            if p and n > 1e-9:
                out[str(p)] = out.get(str(p), 0.0) + float(n)
        return out, float(e.HM.values.squeeze())

    def bisect(self, pred, lo, hi, tol=0.02):
        """pred(T) True at lo, False at hi -> boundary."""
        while hi - lo > tol:
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if pred(mid) else (lo, mid)
        return 0.5 * (lo + hi)

    def single_phase_H(self, phase, T, x=None):
        return float(self.eq(float(T), [phase], x).HM.values.squeeze())

    def liquid_x(self, T):
        if len(self.x) == 1:
            return dict(self.x)
        e = self.eq(float(T))
        ph = list(e.Phase.values.squeeze())
        X = e.X.values.squeeze()[ph.index("LIQUID")]
        comps = [str(c) for c in e.component.values]
        return {c: float(X[i]) for i, c in enumerate(comps)}


def path(al, T):
    rows = []
    for t in T:
        st, h = al.state(t)
        rows.append((t, h, st))
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    db = Database(str(TDB))
    lines = ["# CALPHAD — T_L, T_S e entalpia de fusão das ligas (COST507-modified.tdb, pycalphad)\n",
             f"Base: `{TDB}` (só leitura). Equilíbrio com todas as fases da base para os elementos de cada liga.\n",
             "Definições de ΔH (a escolha é do autor):",
             "- **ΔH_a** = H_LIQ(x₀) − H_FCC(x₀) a T_L(FCC_A1): líquido único × CFC única da composição nominal "
             "(análogo ao ΔH_m de metal puro; é o que a MD de Al puro mede).",
             "- **ΔH_b** = H_eq(T_L) − H_eq(T_S): calor total liberado na solidificação de equilíbrio (latente + sensível).",
             "- **ΔH_c** = ΔH_b − ∫ c_p,sólido dT entre T_S e T_L: só a parte latente "
             "(c_p do conjunto sólido de equilíbrio logo abaixo de T_S).\n"]
    table = ["| liga | T_L [K] | 1º sólido | T_L(FCC_A1) [K] | T_S [K] | ΔH_a (x₀ / x_L) [kJ/kg] | ΔH_b [kJ/kg] | ΔH_c [kJ/kg] "
             "| T_L / ΔH do script |", "|---|---|---|---|---|---|---|---|---|"]
    fig, ax = plt.subplots(figsize=(7, 5))
    notes = []
    for name, wt in ALLOYS.items():
        al = Alloy(db, wt)
        has_solid = lambda T: any(p != "LIQUID" for p in al.state(T)[0])        # noqa: E731
        has_fcc = lambda T: "FCC_A1" in al.state(T)[0]                          # noqa: E731
        has_liq = lambda T: "LIQUID" in al.state(T)[0]                          # noqa: E731
        assert has_solid(700.0) and has_fcc(700.0) and not has_liq(500.0), name
        TL = al.bisect(has_solid, 700.0, 2000.0)
        first = sorted(p for p in al.state(TL - 0.05)[0] if p != "LIQUID")
        TLf = al.bisect(has_fcc, 700.0, TL + 0.05)
        TS = al.bisect(lambda T: not has_liq(T), 500.0, TLf)
        kg = al.Mmol
        dHa = (al.single_phase_H("LIQUID", TLf) - al.single_phase_H("FCC_A1", TLf)) / kg
        xl = al.liquid_x(TLf + 0.05)
        kgl = sum(xl[k] * M[k] for k in xl) * 1e-3
        dHa_l = (al.single_phase_H("LIQUID", TLf, xl) - al.single_phase_H("FCC_A1", TLf, xl)) / kgl
        T_top = TLf
        Hl = al.state(TLf + 0.05)[1]
        Hs = al.state(TS - 0.05)[1]
        dHb = (Hl - Hs) / kg
        cps = (al.state(TS - 0.05)[1] - al.state(TS - 5.05)[1]) / 5.0          # J/mol/K, solid assemblage
        dHc = dHb - cps * (T_top - TS) / kg
        ref = SCRIPT.get(name)
        refs = f"{ref[0]:.2f} / {ref[1]/1e3:.1f}" if ref else "—"
        table.append(f"| {name} | {TL:.2f} | {', '.join(first)} | {TLf:.2f} | {TS:.2f} | {dHa/1e3:.1f} ({dHa_l/1e3:.1f}) | "
                     f"{dHb/1e3:.1f} | {dHc/1e3:.1f} | {refs} |")
        # equilibrium path below T_L(FCC): fractions and H
        T = np.arange(np.floor(TLf) + 5, max(TS - 30, 500), -2.0)
        rows = path(al, T)
        allph = sorted({p for _, _, st in rows for p in st})
        with open(OUT / f"alloys_path_{name.replace(',', '').replace(' ', '_')}.csv", "w", encoding="utf-8") as fh:
            fh.write("T_K,H_J_per_kg," + ",".join(allph) + "\n")
            for t, h, st in rows:
                fh.write(f"{t:.2f},{h/kg:.2f}," + ",".join(f"{st.get(p, 0):.6f}" for p in allph) + "\n")
        ax.plot(T, [(h - Hs) / kg / 1e3 for _, h, _ in rows], label=name)
        st_ts = al.state(TS - 0.05)[0]
        notes.append(f"- {name}: M = {kg*1e3:.4f} g/mol; x = " + ", ".join(f"{k} {x:.5f}" for k, x in al.x.items())
                     + "; líquido a T_L(FCC_A1): " + ", ".join(f"{k} {wl:.4f} wt.%" for k, wl in
                                                         ((k, xl[k] * M[k] / (kgl * 1e3) * 100) for k in xl if k != "AL"))
                     + f"; sólidos logo abaixo de T_S: " + ", ".join(f"{p} {f:.4f}" for p, f in st_ts.items()))
        print(table[-1], flush=True)
    ax.set_xlabel("T [K]"); ax.set_ylabel("H − H(T_S⁻) [kJ/kg]"); ax.legend(); ax.set_title("COST507-modified, equilíbrio")
    fig.savefig(OUT / "alloys_H_T.png", dpi=150, bbox_inches="tight")
    lines += table + ["", "Composições e fases:"] + notes
    (OUT / "alloys_TL_dH.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    sys.exit(main())
