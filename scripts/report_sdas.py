# -*- coding: utf-8 -*-
"""
End-to-end validation chain (continuum), the chain the MD Γ will be plugged into:

    ∇T_exp  ->  FerreiraModel Γ¹ˢᵗ, Γ²ⁿᵈ (= A·∇T)  ->  SDAS(t_SL) (MRB)  ->  experimental SDAS

    python scripts/report_sdas.py  ->  results/sdas/report_sdas.md, sdas_<alloy>.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402
import pandas as pd               # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from nucleation_md import FerreiraModel                                   # noqa: E402
from nucleation_md.sdas import ALLOY_SDAS, sdas_from_gamma, power_law_fit  # noqa: E402
from nucleation_md.sources.continuum import ContinuumSource, ALLOYS       # noqa: E402

OUT = ROOT / "results" / "sdas"


def experimental(alloy):
    return pd.read_csv(ROOT / "data" / "experimental" / f"sdas_{alloy}.csv", comment="#")


def compare(alloy):
    src = ContinuumSource(alloy=alloy, dSv_mode="closure", theta_mode="homogeneous")
    s = src.state_at_gradient(src.grad_T_exp)
    fm = FerreiraModel()
    gammas = {"1st": fm.gamma_tensor(s, 1), "2nd": fm.gamma_tensor(s, 2)}
    system = ALLOY_SDAS[alloy]
    exp = experimental(alloy)
    rows, curves = [], {}
    for label, model, G in (("MRB Γ¹ˢᵗ", "MRB", gammas["1st"]), ("MRB Γ²ⁿᵈ", "MRB", gammas["2nd"]),
                            ("RB Γ²ⁿᵈ", "RB", gammas["2nd"])):
        tSL, S = sdas_from_gamma(G, system, model)
        a, b = power_law_fit(tSL, S)
        pred = a * exp.tSL_s ** b
        ratio = np.log(pred / exp.SDAS_um)
        inside = ((pred <= exp.SDAS_um + exp.err_plus_um) & (pred >= exp.SDAS_um - exp.err_minus_um)).mean()
        rows.append(dict(model=label, Gamma=G, a=a, b=b, rms_log=np.sqrt(np.mean(ratio ** 2)),
                         bias_log=ratio.mean(), within_error_bars=inside))
        curves[label] = (tSL, S * 1e6)
    return src, s, gammas, pd.DataFrame(rows), exp, curves


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    md = ["# Cadeia de validação contínuo → SDAS\n",
          "∇T experimental → Γ (FerreiraModel, Γ = A·∇T) → SDAS (MRB, Ferreira 2024) → pontos experimentais.\n",
          "Métricas: rms_log = √⟨ln²(pred/exp)⟩; bias_log = ⟨ln(pred/exp)⟩; within = fração dos pontos cuja barra de erro contém a previsão.\n"]
    for alloy in ALLOYS:
        src, s, gammas, df, exp, curves = compare(alloy)
        df.to_csv(OUT / f"sdas_{alloy}.csv", index=False)
        md.append(f"\n## {alloy} — ∇T = {src.grad_T_exp} K/m ({ALLOYS[alloy]['experiment']})\n")
        md.append(f"r_c² = {s.r_eval*1e6:.4f} μm, ΔT = {float(s.DT(s.r_eval)):.5f} K, "
                  f"Γ¹ˢᵗ = {gammas['1st']:.5e} m·K, Γ²ⁿᵈ = {gammas['2nd']:.5e} m·K\n")
        md.append("| modelo | Γ [m K] | SDAS = a·t^b | rms_log | bias_log | within |\n|---|---|---|---|---|---|")
        for _, x in df.iterrows():
            md.append(f"| {x.model} | {x.Gamma:.4e} | {x.a:.2f}·t^{x.b:.3f} | {x.rms_log:.3f} | {x.bias_log:+.3f} | {x.within_error_bars:.0%} |")
        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        styles = {"MRB Γ¹ˢᵗ": "--", "MRB Γ²ⁿᵈ": "-.", "RB Γ²ⁿᵈ": "-"}
        for label, (t, S) in curves.items():
            ax.plot(t, S, styles[label], color="k", lw=1.3, label=label)
        ax.errorbar(exp.tSL_s, exp.SDAS_um, yerr=[exp.err_minus_um, exp.err_plus_um], fmt="sk", ms=4,
                    label=ALLOYS[alloy]["experiment"])
        ax.set(xscale="log", yscale="log", xlabel="t_SL [s]", ylabel="SDAS [μm]", title=alloy)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(OUT / f"sdas_{alloy}.png", dpi=130)
        plt.close(fig)
    (OUT / "report_sdas.md").write_text("\n".join(md), encoding="utf-8")
    print("\n".join(md))


if __name__ == "__main__":
    main()
