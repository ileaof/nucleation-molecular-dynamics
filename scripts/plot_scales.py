# -*- coding: utf-8 -*-
"""
Figures for evaluating the continuum and MD solutions of both alloys.

    python scripts/plot_scales.py  ->  results/figures/
        fig_gamma_gradT.png   (a) Γ vs ∇T   (b) ΔT vs r_c      continuum + atomistic scale
        fig_sdas.png          SDAS vs t_SL: experiment, MRB(Γ¹ˢᵗ), MRB(Γ²ⁿᵈ), RB, MRB(Γ_MD)

Continuum curves come from the continuum tables (FerreiraModel, closure, homogeneous).
Atomistic-scale curves are PREDICTIONS of Γ = A·∇T = ΔT·r_c/2 with the same Γ as the
continuum at the experimental gradient (Γ valid at both scales), r_c = 1–20 nm.
MD results, when they exist, are read from results/md/<alloy>/gamma_md.csv
(columns: r_c_m, DT_K, gradT_K_per_m, Gamma_mK) and drawn as black markers.
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

from nucleation_md import FerreiraModel                                    # noqa: E402
from nucleation_md.sdas import ALLOY_SDAS, sdas_from_gamma                  # noqa: E402
from nucleation_md.sources.continuum import ContinuumSource, ALLOYS        # noqa: E402

OUT = ROOT / "results" / "figures"
C1, C2, GRAY, INK, MUTED = "#2a78d6", "#eb6834", "#9a9a94", "#1f1f1e", "#6b6b66"
LABEL = {"Al3Cu5Nb01Fe": "Al-3Cu-5Nb-0.1Fe", "Al08Si06Mg02Fe": "Al-0.8Si-0.6Mg-0.2Fe"}
R_ATOM = np.geomspace(1e-9, 20e-9, 60)
R_MARK = np.array([2e-9, 5e-9, 10e-9, 20e-9])

plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlesize": 10,
                     "axes.titleweight": "bold", "legend.frameon": False})


def continuum(alloy):
    src = ContinuumSource(alloy=alloy, dSv_mode="closure", theta_mode="homogeneous")
    fm = FerreiraModel()
    rows = []
    for i in src.rows:
        s = src.state(i)
        rows.append(dict(r=s.r_eval, DT=float(s.DT(s.r_eval)), gradT=fm.grad_T(s, 2),
                         G1=fm.gamma_tensor(s, 1), G2=fm.gamma_tensor(s, 2)))
    s_exp = src.state_at_gradient(src.grad_T_exp)
    exp = dict(gradT=src.grad_T_exp, r=s_exp.r_eval, DT=float(s_exp.DT(s_exp.r_eval)),
               G1=fm.gamma_tensor(s_exp, 1), G2=fm.gamma_tensor(s_exp, 2), TL=src.params["TL"])
    return pd.DataFrame(rows), exp


def md_results(alloy):
    path = ROOT / "results" / "md" / alloy / "gamma_md.csv"
    return pd.read_csv(path) if path.exists() else None


def grid(ax):
    ax.grid(True, which="major", color="#e4e4df", lw=0.6)
    ax.set_axisbelow(True)


def fig_gamma_gradT():
    fig, axes = plt.subplots(2, 2, figsize=(10, 8.2))
    for row, alloy in enumerate(ALLOYS):
        df, exp = continuum(alloy)
        md = md_results(alloy)
        ax = axes[row, 0]
        # iso-radius guides Γ = 4πr²∇T
        g = np.geomspace(1e2, 1e11, 50)
        for r, lab in ((exp["r"], f"r = {exp['r']*1e6:.2f} μm"), (10e-9, "r = 10 nm"), (2e-9, "r = 2 nm")):
            ax.plot(g, 4 * np.pi * r ** 2 * g, ":", color=GRAY, lw=0.9)
            gx = 2e-7 / (4 * np.pi * r ** 2)
            if 1e2 < gx < 1e11:
                ax.annotate(lab, (gx, 2e-7), color=MUTED, fontsize=7, rotation=0,
                            xytext=(6, -2), textcoords="offset points",
                            bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none"))
        ax.plot(df.gradT, df.G1, "--", color=C1, lw=2, label="Γ¹ˢᵗ contínuo")
        ax.plot(df.gradT, df.G2, "-", color=C2, lw=2, label="Γ²ⁿᵈ contínuo")
        for G, c, key in ((exp["G1"], C1, "1st"), (exp["G2"], C2, "2nd")):
            ax.axhline(G, color=c, lw=0.7, alpha=0.6)
            ga = G / (4 * np.pi * R_MARK ** 2)
            ax.plot(ga, np.full_like(ga, G), "o", mfc="white", mec=c, mew=1.6, ms=8,
                    label=f"previsão atomística Γ{'¹ˢᵗ' if key == '1st' else '²ⁿᵈ'} (r_c = 2, 5, 10, 20 nm)")
        for r, ga in zip(R_MARK, exp["G2"] / (4 * np.pi * R_MARK ** 2)):
            ax.annotate(f"{r*1e9:.0f} nm", (ga, exp["G2"]), xytext=(2, 9), textcoords="offset points",
                        ha="left", rotation=45, fontsize=7, color=INK)
        ax.axvline(exp["gradT"], color=INK, lw=0.8, ls="-.")
        ax.annotate(f"∇T exp = {exp['gradT']:g} K/m", (exp["gradT"], 3e-9), rotation=90,
                    xytext=(-10, 0), textcoords="offset points", fontsize=7, color=INK)
        if md is not None:
            ax.plot(md.gradT_K_per_m, md.Gamma_mK, "*", color=INK, ms=11, label="MD (seeding + NEMD)")
        ax.set(xscale="log", yscale="log", xlim=(1e2, 1e11), ylim=(1e-9, 1e-5),
               xlabel="Gradiente térmico ∇T [K/m]", ylabel="Gibbs–Thomson Γ = A·∇T [m·K]",
               title=f"({'ac'[row]}) {LABEL[alloy]}: Γ × ∇T")
        grid(ax)
        ax.legend(fontsize=7, loc="lower right")

        ax = axes[row, 1]
        ax.axhspan(exp["TL"], 1e5, color="#f2f2ee")
        ax.annotate("ΔT > T_L: inacessível", (2.5e-8, exp["TL"] * 3), fontsize=7, color=MUTED)
        ax.plot(df.r, df.DT, "-", color=C2, lw=2, label="contínuo (ΔT = 2Γ²ⁿᵈ/r_c)")
        ax.plot(exp["r"], exp["DT"], "s", color=C2, ms=8, mec="white", mew=1.2, label="ponto experimental")
        for G, c, ls, lab in ((exp["G1"], C1, "--", "Γ¹ˢᵗ"), (exp["G2"], C2, "-", "Γ²ⁿᵈ")):
            ax.plot(R_ATOM, 2 * G / R_ATOM, ls, color=c, lw=1.6, alpha=0.9,
                    label=f"previsão atomística, {lab} = {G:.3g} m·K")
            ax.plot(R_MARK, 2 * G / R_MARK, "o", mfc="white", mec=c, mew=1.6, ms=7)
        if md is not None:
            ax.plot(md.r_c_m, md.DT_K, "*", color=INK, ms=11, label="MD (seeding)")
        ax.set(xscale="log", yscale="log", xlim=(1e-9, 1e-5), ylim=(1e-2, 1e5),
               xlabel="Raio crítico r_c [m]", ylabel="Super-resfriamento ΔT [K]",
               title=f"({'bd'[row]}) {LABEL[alloy]}: ΔT × r_c")
        grid(ax)
        ax.legend(fontsize=7, loc="lower left")
    fig.suptitle("Γ = A·∇T = ΔT·r_c/2 — contínuo (μm) e escala atomística (nm), mesmo Γ", fontsize=11, color=INK)
    fig.text(0.5, 0.005, "Marcadores vazados: previsão com Γ invariante entre escalas (sem MD ainda). "
             "Estrelas pretas aparecem quando results/md/<liga>/gamma_md.csv existir.",
             ha="center", fontsize=7.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    fig.savefig(OUT / "fig_gamma_gradT.png", dpi=160)
    plt.close(fig)


def fig_sdas():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4))
    for ax, alloy in zip(axes, ALLOYS):
        _, exp = continuum(alloy)
        system = ALLOY_SDAS[alloy]
        data = pd.read_csv(ROOT / "data" / "experimental" / f"sdas_{alloy}.csv", comment="#")
        t, S2 = sdas_from_gamma(exp["G2"], system)
        _, S1 = sdas_from_gamma(exp["G1"], system)
        _, SR = sdas_from_gamma(exp["G2"], system, "RB")
        ax.plot(t, SR * 1e6, "-", color=GRAY, lw=1.6, label="RB (Rappaz–Boettinger)")
        ax.plot(t, S1 * 1e6, "--", color=C1, lw=2, label=f"MRB contínuo Γ¹ˢᵗ = {exp['G1']:.3g}")
        ax.plot(t, S2 * 1e6, "-", color=C2, lw=2, label=f"MRB contínuo Γ²ⁿᵈ = {exp['G2']:.3g}")
        md = md_results(alloy)
        if md is not None:
            G_md = float(np.mean(md.Gamma_mK))
            _, Sm = sdas_from_gamma(G_md, system)
            ax.plot(t, Sm * 1e6, "-", color=INK, lw=2, label=f"MRB com Γ_MD = {G_md:.3g}")
        else:
            ax.annotate("Γ_MD: aguardando simulações\n(com Γ invariante, coincide com a curva Γ²ⁿᵈ)",
                        (0.03, 0.97), xycoords="axes fraction", va="top", fontsize=7.5, color=MUTED)
        ax.errorbar(data.tSL_s, data.SDAS_um, yerr=[data.err_minus_um, data.err_plus_um],
                    fmt="s", color=INK, ms=4.5, elinewidth=0.9, capsize=2,
                    label=ALLOYS[alloy]["experiment"])
        ax.set(xscale="log", yscale="log", ylim=(5, 400), xlabel="Tempo local de solidificação t_SL [s]",
               ylabel="SDAS [μm]", title=f"{LABEL[alloy]} (∇T = {exp['gradT']:g} K/m)")
        grid(ax)
        ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT / "fig_sdas.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    fig_gamma_gradT()
    fig_sdas()
    print(OUT)
