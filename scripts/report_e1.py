# -*- coding: utf-8 -*-
"""
Stage E1 report on the continuum solution (no MD yet).

    python scripts/report_e1.py   ->  results/e1/*.csv and results/e1/report_e1.md

1. exact roots of a r² + 3b r + 6γf = 0 versus r_c² = -3b/(2a)
2. CNT x Ferreira (1st, 2nd order) critical radii and barriers
3. decomposition of the r_c¹ denominator: ΔG_V, ∂ΔG_S/∂r, ΔG_C
4. sign of ΔG* and trend of I with ΔT and ∇T (literal exponent +ΔG*/ΔG*_Eq)
for each continuum ΔS_V mode (closure / integrated / constant).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from nucleation_md import CNT, FerreiraModel                       # noqa: E402
from nucleation_md.geometry import cap_area                        # noqa: E402
from nucleation_md.sources.continuum import ContinuumSource, dself_legacy   # noqa: E402

OUT = ROOT / "results" / "e1"


def table_for(dSv_mode="closure", theta_mode="homogeneous"):
    src = ContinuumSource(dSv_mode=dSv_mode, theta_mode=theta_mode)
    ref = src.reference_state()
    fm, cnt = FerreiraModel(reference=ref, order=2), CNT()
    lam = src.params["ad"]
    rows = []
    for i in src.rows:
        s = src.state(i)
        k = fm.coefficients(s)
        r_m, r_p, disc = fm.exact_roots(s)
        rc1, rc2 = fm.r_critical(s, 1), fm.r_critical(s, 2)
        T = src.params["TL"] - k.DT
        D = dself_legacy(T, src.params["TL"])
        # prefactor with N/V = 1 (relative trend only; N/V is not varied along the table)
        I_rel = fm.rate(s, D=D, lam=lam, N_V=1.0)
        rows.append(dict(
            i=i, r=s.r_eval, DT=k.DT, gradT_table=float(src.row(i)["dTdr"]),
            a=k.a, b=k.b, disc=disc,
            rc2_vertex=rc2,
            root_minus=np.real(r_m), root_plus=np.real(r_p),
            root_imag=abs(np.imag(r_p)),
            rc1=rc1, rc_CNT=cnt.r_critical(s),
            dG_V=k.dG_V, ddG_S_dr=k.ddG_S_dr, dG_conf=k.dG_conf, denominator=k.denominator,
            G_star_rc2=fm.delta_G(rc2, s), G_star_rc1=fm.delta_G(rc1, s), G_star_CNT=cnt.barrier(s),
            G_eq=fm.dG_eq, Gamma1=fm.gamma_tensor(s, 1), Gamma2=fm.gamma_tensor(s, 2),
            gradT_model=fm.grad_T(s, 2),
            I_rel=I_rel,
        ))
    return pd.DataFrame(rows)


def trend(x, y):
    """Spearman-like sign of monotonic trend: fraction of positive consecutive slopes."""
    order = np.argsort(x)
    dy = np.diff(np.asarray(y)[order])
    return float(np.mean(dy > 0))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    md = ["# E1 – relatório sobre a solução do contínuo\n",
          "Gerado por `scripts/report_e1.py`. Fonte: `tablefull_paper_B_MRS_2026.csv` (71 linhas do laço).\n"]
    for mode in ("closure",):              # the continuum's own construction
        df = table_for(mode)
        df.to_csv(OUT / f"e1_continuum_{mode}.csv", index=False)
        sel = df[df.i.isin([1, 10, 20, 38, 50, 71])]
        md.append(f"\n## ΔS_V modo `{mode}`\n")
        md.append("### Raízes exatas × vértice r_c² = −3b/(2a)\n")
        md.append("| i | r [μm] | r_c² vértice [μm] | r₋ [μm] | r₊ [μm] | Im [μm] | disc |\n|---|---|---|---|---|---|---|")
        for _, x in sel.iterrows():
            md.append(f"| {x.i} | {x.r*1e6:.4f} | {x.rc2_vertex*1e6:.4f} | {x.root_minus*1e6:.4f} | "
                      f"{x.root_plus*1e6:.4f} | {x.root_imag*1e6:.3g} | {x.disc:.3e} |")
        n_complex = int((df.disc < 0).sum())
        rel = np.abs(df[["root_minus", "root_plus"]].sub(df.rc2_vertex, axis=0)).min(axis=1) / df.rc2_vertex
        md.append(f"\nDiscriminante < 0 (raízes complexas) em {n_complex}/{len(df)} linhas. "
                  f"Distância relativa da raiz mais próxima ao vértice: mín {rel.min():.3f}, máx {rel.max():.3f}.\n")
        md.append("### Raios e barreiras\n")
        md.append("| i | r_CNT [μm] | r_c¹ [μm] | r_c² [μm] | ΔG*_CNT [J] | ΔG*(r_c¹) [J] | ΔG*(r_c²) [J] |\n|---|---|---|---|---|---|---|")
        for _, x in sel.iterrows():
            md.append(f"| {x.i} | {x.rc_CNT*1e6:.4f} | {x.rc1*1e6:.4f} | {x.rc2_vertex*1e6:.4f} | "
                      f"{x.G_star_CNT:.3e} | {x.G_star_rc1:.3e} | {x.G_star_rc2:.3e} |")
        neg = df[df.G_star_rc2 <= 0].i.tolist()
        md.append(f"\nΔG*(r_c²) ≤ 0 em {len(neg)} linhas" + (f" (a partir de i={neg[0]})." if neg else "."))
        md.append("\n### Decomposição do denominador de r_c¹ (J/m³)\n")
        md.append("| i | ΔG_V | ∂ΔG_S/∂r | ΔG_C | soma |\n|---|---|---|---|---|")
        for _, x in sel.iterrows():
            md.append(f"| {x.i} | {x.dG_V:.4e} | {x.ddG_S_dr:.4e} | {x.dG_conf:.4e} | {x.denominator:.4e} |")
        md.append("\n### Tensor de campo térmico Γ = A·∇T\n")
        md.append("| i | ∇T tabela [K/m] | ∇T = Γ²ⁿᵈ/A [K/m] | Γ¹ˢᵗ [m K] | Γ²ⁿᵈ [m K] |\n|---|---|---|---|---|")
        for _, x in sel.iterrows():
            md.append(f"| {x.i} | {x.gradT_table:.2f} | {x.gradT_model:.2f} | {x.Gamma1:.4e} | {x.Gamma2:.4e} |")
        md.append(f"\n### Tendência de I (expoente literal +ΔG*/ΔG*_Eq, referência i=1)\n")
        md.append(f"- fração de inclinações positivas de I contra ΔT: {trend(df.DT, df.I_rel):.2f}")
        md.append(f"- fração de inclinações positivas de I contra ∇T (coluna dTdr): {trend(df.gradT_table, df.I_rel):.2f}")
        md.append(f"- I(i=71)/I(i=1) = {df.I_rel.iloc[-1] / df.I_rel.iloc[0]:.3e}\n")
    (OUT / "report_e1.md").write_text("\n".join(md), encoding="utf-8")
    print("\n".join(md))


if __name__ == "__main__":
    main()
