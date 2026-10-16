# -*- coding: utf-8 -*-
"""
Continuum × molecular dynamics figures.

1. SDAS × t_SL for both alloys: MRB/RB with the continuum Γ (Γ¹ˢᵗ, Γ²ⁿᵈ at the experimental ∇T), MRB with the MD Γ
   (seeding at ∇T = 0, md/e4_seeding: Γ = ΔT·r_c/2 = 1.219e-7 m·K, pure Al, Borovikov 2024), experimental points.
2. Solid–liquid interfacial energy γ_SL × radius: continuum γ(r) (column gam_hom of the continuum tables, μm radii),
   MD γ(r_c) = Γ(r_c)·|ΔS_V| from each seed radius (nm), and the planar γ₀ from capillary fluctuations (r → ∞).

    python scripts/plot_md_vs_continuum.py  ->  results/figures/fig_sdas_md_continuum.png, fig_gamma_md_continuum.png
"""
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402
import pandas as pd               # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from nucleation_md import FerreiraModel                                       # noqa: E402
from nucleation_md.sdas import ALLOY_SDAS, power_law_fit, sdas_from_gamma      # noqa: E402
from nucleation_md.sources.continuum import ContinuumSource                    # noqa: E402

OUT = ROOT / "results" / "figures"
ALLOYS = {"Al3Cu5Nb01Fe": ("Al-3Cu-5Nb-0,1Fe", "Mendes et al. (2023)"),
          "Al08Si06Mg02Fe": ("Al-0,8Si-0,6Mg-0,2Fe", "Marques et al. (2025)")}
# reference categorical palette (dataviz skill, light mode, fixed order) + text inks
C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK2, GRID = "#1f1f1e", "#5f5e58", "#e4e3dc"
# MD inputs (E2 and E4): ΔS_V at T_m (md/e2_melting/E2_results.md); γ₀ ranges of md/e2_cfm_thick/E2_thick_results.md
DSV_MD = -8.878e5                 # J m⁻³ K⁻¹, Borovikov 2024
GAMMA0_BORO = (0.105, 0.120)      # J/m², fitas finas e espessas, janelas 2–20
GAMMA0_NEP = (0.115, 0.125)


def style(ax):
    ax.grid(True, which="major", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelcolor=INK)


def md_seeding():
    txt = (ROOT / "md" / "e4_seeding" / "seeding_results.md").read_text(encoding="utf-8")
    rows = re.findall(r"^\| (\d+) \| \d+ \| \d+ \| ([\d.]+) \| ([\d.]+) \| ([\d.]+) \| ([\d.e+-]+) \|$", txt, re.M)
    d = pd.DataFrame([[float(x) for x in r] for r in rows], columns=["R_A", "Tstar", "rc_nm", "dT", "Gamma"])
    return d


def fig_sdas(G_md):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    lines = []
    for ax, (alloy, (name, ref)) in zip(axes, ALLOYS.items()):
        src = ContinuumSource(alloy=alloy, dSv_mode="closure", theta_mode="homogeneous")
        s = src.state_at_gradient(src.grad_T_exp)
        fm = FerreiraModel()
        g1, g2 = fm.gamma_tensor(s, 1), fm.gamma_tensor(s, 2)
        system = ALLOY_SDAS[alloy]
        exp = pd.read_csv(ROOT / "data" / "experimental" / f"sdas_{alloy}.csv", comment="#")
        curves = [("RB, Γ²ⁿᵈ contínuo", "RB", g2, C3, (0, (1, 1.5))),
                  ("MRB, Γ²ⁿᵈ contínuo", "MRB", g2, C2, (0, (6, 2))),
                  ("MRB, Γ¹ˢᵗ contínuo", "MRB", g1, C1, (0, (3, 1.5, 1, 1.5))),
                  ("MRB, Γ da MD", "MRB", G_md, C4, "-")]
        for label, model, G, c, ls in curves:
            t, S = sdas_from_gamma(G, system, model)
            a, b = power_law_fit(t, S)
            pred = a * exp.tSL_s ** b
            inside = ((pred <= exp.SDAS_um + exp.err_plus_um) & (pred >= exp.SDAS_um - exp.err_minus_um)).mean()
            ax.plot(t, S * 1e6, color=c, lw=2, ls=ls, label=f"{label} ({G:.2e} m·K)")
            lines.append(f"{name} | {label} | Γ = {G:.3e} m·K | SDAS = {a:.2f}·t^{b:.3f} μm | dentro das barras: {inside:.0%}")
        ax.errorbar(exp.tSL_s, exp.SDAS_um, yerr=[exp.err_minus_um, exp.err_plus_um], fmt="s", color=INK, ms=6,
                    mfc=INK, mec="white", mew=1, elinewidth=1.2, capsize=0, label=f"experimental, {ref}", zorder=5)
        ax.set(xscale="log", yscale="log")
        ax.set_xlabel("tempo local de solidificação t_SL [s]", color=INK)
        ax.set_title(name, color=INK, fontsize=12, loc="left")
        style(ax)
        ax.legend(fontsize=8, frameon=False, loc="upper left", labelcolor=INK)
    axes[0].set_ylabel("SDAS [μm]", color=INK)
    fig.suptitle("SDAS: Γ do contínuo (μm) × Γ da dinâmica molecular (nm, Al puro, ∇T = 0)", color=INK, fontsize=12)
    fig.text(0.01, 0.005, "Γ da MD: seeding com núcleos de 1,4–2,9 nm, Γ = ΔT·r_c/2 (md/e4_seeding). "
             "Mesma cinética experimental (V_L, t_SL) e mesmos solutos em todas as curvas.", fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    fig.savefig(OUT / "fig_sdas_md_continuum.png", dpi=150)
    plt.close(fig)
    return lines


def fig_gamma(md):
    fig, ax = plt.subplots(figsize=(8, 5))
    for (alloy, (name, _)), c, ls in zip(ALLOYS.items(), (C1, C2), ("-", (0, (6, 2)))):
        src = ContinuumSource(alloy=alloy)
        t = src.table
        r, g = t["r_hom_2nd"].to_numpy(), t["gam_hom"].to_numpy()
        o = np.argsort(r)
        ax.plot(r[o], g[o], color=c, lw=2, ls=ls, label=f"contínuo, {name} (γ de Gurtin–Murdoch)")
        st = src.state_at_gradient(src.grad_T_exp)
        ax.plot(st.r_eval, float(src.gamma(st.r_eval)), "o", ms=9, color=c, mec="white", mew=2, zorder=5)
        ax.annotate(f"∇T exp. = {src.grad_T_exp:g} K/m", (st.r_eval, float(src.gamma(st.r_eval))),
                    textcoords="offset points", xytext=(-12, 8), ha="right", fontsize=8, color=INK2)
    gmd = md.Gamma * abs(DSV_MD)
    ax.plot(md.rc_nm * 1e-9, gmd, "o", ms=9, color=C3, mec="white", mew=2, zorder=5,
            label="MD, seeding: γ(r_c) = Γ·|ΔS_V| (Al puro, Borovikov)")
    from matplotlib.patches import Patch
    ax.axhspan(*GAMMA0_BORO, color=C3, alpha=0.25, lw=0)
    ax.axhspan(*GAMMA0_NEP, color=C4, alpha=0.25, lw=0)
    h, l = ax.get_legend_handles_labels()
    h += [Patch(color=C3, alpha=0.25), Patch(color=C4, alpha=0.25)]
    l += ["MD, interface plana (flutuação capilar), Borovikov: 0,105–0,120 J/m²",
          "MD, interface plana (flutuação capilar), NEP89: 0,115–0,125 J/m²"]
    ax.set(xscale="log", yscale="log", xlim=(1e-9, 1e-5), ylim=(0.05, 3))
    ax.set_xlabel("raio do núcleo r [m]", color=INK)
    ax.set_ylabel("energia de superfície sólido–líquido γ_SL [J/m²]", color=INK)
    ax.set_title("γ_SL: contínuo (μm) × dinâmica molecular (nm)", color=INK, fontsize=12, loc="left")
    style(ax)
    ax.legend(h, l, fontsize=8, frameon=False, loc="upper left", labelcolor=INK)
    fig.text(0.01, 0.005, f"MD: ΔS_V = {DSV_MD:.3e} J m⁻³K⁻¹ (E2, a T_m) suposto constante até T*; "
             "contínuo: colunas r_hom_2nd e gam_hom das tabelas dos scripts.", fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(OUT / "fig_gamma_md_continuum.png", dpi=150)
    plt.close(fig)
    return [f"MD seeding r_c = {r:.3f} nm: γ = {g:.4f} J/m²" for r, g in zip(md.rc_nm, gmd)]


def fig_confronto(G_md, G_lo, G_hi):
    """SDAS continuum × MD: MRB with the continuum Γ¹ˢᵗ/Γ²ⁿᵈ and with the MD Γ (band = range of all MD estimates),
    experimental points; lower row: SDAS_MD / SDAS_continuum."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.2), sharex="col", gridspec_kw=dict(height_ratios=(3, 1.3)))
    lines = []
    for j, (alloy, (name, ref)) in enumerate(ALLOYS.items()):
        ax, axr = axes[0, j], axes[1, j]
        src = ContinuumSource(alloy=alloy, dSv_mode="closure", theta_mode="homogeneous")
        s = src.state_at_gradient(src.grad_T_exp)
        fm = FerreiraModel()
        g1, g2 = fm.gamma_tensor(s, 1), fm.gamma_tensor(s, 2)
        system = ALLOY_SDAS[alloy]
        exp = pd.read_csv(ROOT / "data" / "experimental" / f"sdas_{alloy}.csv", comment="#")
        t, S1 = sdas_from_gamma(g1, system, "MRB")
        _, S2 = sdas_from_gamma(g2, system, "MRB")
        _, Sm = sdas_from_gamma(G_md, system, "MRB")
        _, Slo = sdas_from_gamma(G_lo, system, "MRB")
        _, Shi = sdas_from_gamma(G_hi, system, "MRB")
        _, R2 = sdas_from_gamma(g2, system, "RB")      # Rappaz–Boettinger (equilibrium), as in the author's SDAS plots
        _, Rm = sdas_from_gamma(G_md, system, "RB")
        ax.plot(t, R2 * 1e6, color=INK, lw=1.8, label=f"RB (equilíbrio), contínuo Γ²ⁿᵈ = {g2:.2e} m·K")
        ax.annotate("RB, Γ²ⁿᵈ contínuo", (t[-1], R2[-1] * 1e6), textcoords="offset points", xytext=(-4, 6), ha="right", fontsize=8, color=INK)
        ax.plot(t, Rm * 1e6, color=INK, lw=1.8, ls=(0, (6, 3)), label=f"RB (equilíbrio), dinâmica molecular Γ = {G_md:.2e} m·K")
        ax.annotate("RB, Γ da MD", (t[-1], Rm[-1] * 1e6), textcoords="offset points", xytext=(-4, 6), ha="right", fontsize=8, color=INK)
        ax.plot(t, S2 * 1e6, color=C2, lw=2, label=f"contínuo Γ²ⁿᵈ = {g2:.2e} m·K, MRB")
        ax.plot(t, S1 * 1e6, color=C1, lw=2, label=f"contínuo Γ¹ˢᵗ = {g1:.2e} m·K, MRB")
        ax.fill_between(t, Slo * 1e6, Shi * 1e6, color=C4, alpha=0.3, lw=0)
        ax.plot(t, Sm * 1e6, color=C4, lw=2.4, label=f"dinâmica molecular Γ = {G_lo:.2e}–{G_hi:.2e} m·K, MRB")
        ax.errorbar(exp.tSL_s, exp.SDAS_um, yerr=[exp.err_minus_um, exp.err_plus_um], fmt="s", color=INK, ms=6,
                    mfc=INK, mec="white", mew=1, elinewidth=1.2, capsize=0, label=f"experimental, {ref}", zorder=5)
        ax.set(xscale="log", yscale="log")
        ax.set_title(name, color=INK, fontsize=12, loc="left")
        style(ax)
        ax.legend(fontsize=8, frameon=False, loc="upper left", labelcolor=INK)
        axr.axhline(1, color=INK2, lw=1)
        axr.plot(t, Sm / S1, color=C1, lw=2, label="MRB: MD / contínuo Γ¹ˢᵗ")
        axr.plot(t, Sm / S2, color=C2, lw=2, label="MRB: MD / contínuo Γ²ⁿᵈ")
        axr.plot(t, Rm / R2, color=INK, lw=1.8, ls=(0, (6, 3)), label="RB: MD / contínuo Γ²ⁿᵈ")
        axr.set(xscale="log", ylim=(0.3, 1.05))
        axr.set_xlabel("tempo local de solidificação t_SL [s]", color=INK)
        style(axr)
        for y, lab in ((Sm / S1, "MRB, MD / Γ¹ˢᵗ"), (Sm / S2, "MRB, MD / Γ²ⁿᵈ"), (Rm / R2, "RB, MD / Γ²ⁿᵈ")):
            axr.annotate(f"{lab}: {np.mean(y):.2f}", (t[0], y[0]), textcoords="offset points", xytext=(4, 4),
                         fontsize=8, color=INK)
        win = (t >= exp.tSL_s.min()) & (t <= exp.tSL_s.max())
        for lab, S in (("MRB Γ¹ˢᵗ contínuo", S1), ("MRB Γ²ⁿᵈ contínuo", S2), ("MRB Γ MD", Sm),
                       ("RB Γ²ⁿᵈ contínuo", R2), ("RB Γ MD", Rm)):
            a, b = power_law_fit(t, S)
            pred = a * exp.tSL_s ** b
            inside = ((pred <= exp.SDAS_um + exp.err_plus_um) & (pred >= exp.SDAS_um - exp.err_minus_um)).mean()
            rms = np.sqrt(np.mean(np.log(pred / exp.SDAS_um) ** 2)); bias = np.mean(np.log(pred / exp.SDAS_um))
            lines.append(f"| {name} | {lab} | {a:.2f}·t^{b:.3f} | {inside:.0%} | {rms:.3f} | {bias:+.3f} |")
        lines.append(f"| {name} | razão SDAS_MD/contínuo na janela experimental | Γ¹ˢᵗ: {np.mean(Sm[win]/S1[win]):.2f} | "
                     f"Γ²ⁿᵈ: {np.mean(Sm[win]/S2[win]):.2f} | RB Γ²ⁿᵈ: {np.mean(Rm[win]/R2[win]):.2f} | |")
    axes[0, 0].set_ylabel("SDAS [μm]", color=INK)
    axes[1, 0].set_ylabel("SDAS_MD / SDAS_contínuo", color=INK)
    fig.suptitle("Confronto SDAS: contínuo (Γ na escala μm) × dinâmica molecular (Γ na escala nm)", color=INK, fontsize=12)
    fig.text(0.01, 0.005, "MD: Al puro, Borovikov 2024; Γ = ΔT·r_c/2 do seeding (∇T = 0, r_c = 1,4–2,9 nm) e da NEMD (∇T até 4,8e9 K/m). "
             "Mesma cinética experimental e solutos em todas as curvas.\n"
             "Linhas coloridas = MRB (cor = origem de Γ); pretas = RB (Rappaz–Boettinger, equilíbrio): contínua = Γ²ⁿᵈ do contínuo, tracejada = Γ da MD.", fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0.045, 1, 0.96))
    fig.savefig(OUT / "fig_sdas_confronto_md_continuo.png", dpi=150)
    plt.close(fig)
    return lines


def fig_gap(G_md, G_lo, G_hi):
    """Γ × ∇T: continuum Γ¹ˢᵗ and Γ²ⁿᵈ of both alloys (table rows) and the MD equilibrium Γ (∇T = 0, horizontal band).
    The 'gap' ∇T_gap is where the continuum Γ¹ˢᵗ(∇T) reaches the MD equilibrium value."""
    fig, ax = plt.subplots(figsize=(8, 5.2))
    lines = []
    for (alloy, (name, _)), c in zip(ALLOYS.items(), (C1, C2)):
        src = ContinuumSource(alloy=alloy)
        t = src.table[src.table.dTdr > 0].sort_values("dTdr")
        g, G1, G2 = t.dTdr.values, t.GT_het.values, t.GT_het_2nd.values
        ax.plot(g, G1, color=c, lw=2, label=f"contínuo Γ¹ˢᵗ, {name}")
        ax.plot(g, G2, color=c, lw=2, ls=(0, (6, 3)), label=f"contínuo Γ²ⁿᵈ, {name}")
        o = np.argsort(G1)
        ggap = float(np.interp(G_md, G1[o], g[o]))
        ax.plot(ggap, G_md, "o", ms=9, color=c, mec="white", mew=2, zorder=5)
        ax.annotate(f"∇T_gap = {ggap:.0f} K/m", (ggap, G_md), textcoords="offset points",
                    xytext=(8, -14 if c == C1 else 6), fontsize=8, color=INK)
        ax.axvline(src.grad_T_exp, color=c, lw=1, ls=":")
        ax.annotate(f"∇T exp. {name}", (src.grad_T_exp, 1.6e-8), rotation=90, fontsize=7, color=INK2,
                    textcoords="offset points", xytext=(-9, 0))
        a, b = np.polyfit(g, G1, 1)
        lines.append(f"{name}: Γ¹ˢᵗ ≈ {b:.3e} + {a:.3e}·∇T; ∇T_gap (interp.) = {ggap:.0f} K/m")
    ax.axhspan(G_lo, G_hi, color=C4, alpha=0.35, lw=0)
    ax.axhline(G_md, color=C4, lw=2.4, label=f"dinâmica molecular, equilíbrio (∇T = 0): Γ = {G_md:.2e} m·K")
    ax.set(xscale="log", yscale="log", xlim=(300, 1.2e4), ylim=(1e-8, 3e-6))
    ax.set_xlabel("gradiente do campo térmico ∇T [K/m]", color=INK)
    ax.set_ylabel("tensor de campo térmico Γ [m·K]", color=INK)
    ax.set_title("Gap entre o equilíbrio da MD e o contínuo: Γ¹ˢᵗ(∇T_gap) = Γ_MD", color=INK, fontsize=12, loc="left")
    style(ax)
    ax.legend(fontsize=8, frameon=False, loc="upper left", labelcolor=INK)
    fig.tight_layout()
    fig.savefig(OUT / "fig_gap_gamma_gradT.png", dpi=150)
    plt.close(fig)
    return lines


def fig_gap_sdas(G_md):
    """Gap defined by the author (2026-10-07): the difference between the continuum RB (Γ²ⁿᵈ) and the MD RB curves.
    In log–log it is a constant vertical shift, k = SDAS_RB,cont / SDAS_RB,MD (independent of t_SL). The MD RB
    (equilibrium Gibbs–Thomson of the MD) is moved up by k onto the continuum RB, and the MD MRB is moved up by the
    same k; the unshifted MD MRB is kept as a thin line for reference."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    lines = []
    for ax, (alloy, (name, ref)) in zip(axes, ALLOYS.items()):
        src = ContinuumSource(alloy=alloy, dSv_mode="closure", theta_mode="homogeneous")
        s = src.state_at_gradient(src.grad_T_exp)
        fm = FerreiraModel()
        g1, g2 = fm.gamma_tensor(s, 1), fm.gamma_tensor(s, 2)
        system = ALLOY_SDAS[alloy]
        exp = pd.read_csv(ROOT / "data" / "experimental" / f"sdas_{alloy}.csv", comment="#")
        t, S1 = sdas_from_gamma(g1, system, "MRB")
        _, S2 = sdas_from_gamma(g2, system, "MRB")
        _, Sm = sdas_from_gamma(G_md, system, "MRB")
        _, R2 = sdas_from_gamma(g2, system, "RB")
        _, Rm = sdas_from_gamma(G_md, system, "RB")
        k = float(np.mean(R2 / Rm))                     # constant: (Γ²ⁿᵈ/Γ_MD)^(1/3)
        ax.plot(t, R2 * 1e6, color=INK2, lw=4, alpha=0.45, label=f"RB (equilíbrio), contínuo Γ²ⁿᵈ = {g2:.2e} m·K")
        ax.plot(t, k * Rm * 1e6, color=INK, lw=1.8, ls=(0, (6, 3)), zorder=4,
                label=f"RB (equilíbrio), MD Γ = {G_md:.2e} m·K, deslocada × {k:.2f}")
        ax.plot(t, S2 * 1e6, color=C2, lw=2, label=f"contínuo Γ²ⁿᵈ, MRB")
        ax.plot(t, S1 * 1e6, color=C1, lw=2, label=f"contínuo Γ¹ˢᵗ = {g1:.2e} m·K, MRB")
        ax.plot(t, Sm * 1e6, color=C4, lw=1.2, ls=(0, (2, 2)), label="MD, MRB, sem deslocamento")
        ax.plot(t, k * Sm * 1e6, color=C4, lw=2.6, label=f"MD, MRB, deslocada × {k:.2f} (gap)")
        ax.errorbar(exp.tSL_s, exp.SDAS_um, yerr=[exp.err_minus_um, exp.err_plus_um], fmt="s", color=INK, ms=6,
                    mfc=INK, mec="white", mew=1, elinewidth=1.2, capsize=0, label=f"experimental, {ref}", zorder=6)
        ax.set(xscale="log", yscale="log")
        ax.set_xlabel("tempo local de solidificação t_SL [s]", color=INK)
        ax.set_title(f"{name} — gap: × {k:.2f} (Δlog₁₀ = {np.log10(k):.3f})", color=INK, fontsize=11, loc="left")
        style(ax)
        ax.legend(fontsize=7.5, frameon=False, loc="upper left", labelcolor=INK)
        for lab, S in (("MRB Γ¹ˢᵗ contínuo", S1), ("MRB Γ²ⁿᵈ contínuo", S2), ("MRB MD sem gap", Sm),
                       ("MRB MD com gap", k * Sm)):
            a, b = power_law_fit(t, S)
            pred = a * exp.tSL_s ** b
            inside = ((pred <= exp.SDAS_um + exp.err_plus_um) & (pred >= exp.SDAS_um - exp.err_minus_um)).mean()
            rms = np.sqrt(np.mean(np.log(pred / exp.SDAS_um) ** 2)); bias = np.mean(np.log(pred / exp.SDAS_um))
            lines.append(f"| {name} | {lab} | {a:.2f}·t^{b:.3f} | {inside:.0%} | {rms:.3f} | {bias:+.3f} |")
        lines.append(f"| {name} | gap RB contínuo/RB MD | × {k:.3f} (Γ²ⁿᵈ/Γ_MD = {g2/G_md:.2f}; ^(1/3) = {(g2/G_md)**(1/3):.3f}) | | | |")
    axes[0].set_ylabel("SDAS [μm]", color=INK)
    fig.suptitle("Gap = RB contínuo − RB MD (constante em log): equilíbrio da MD levado ao do contínuo, demais curvas MD deslocadas",
                 color=INK, fontsize=11)
    fig.text(0.01, 0.005, "RB (Rappaz–Boettinger, equilíbrio): cinza largo = contínuo, preta tracejada = MD deslocada pelo gap (coincidem). "
             "Amarelo grosso = MRB da MD deslocada pelo mesmo gap; amarelo fino = sem deslocamento.", fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    fig.savefig(OUT / "fig_sdas_gap_md_continuo.png", dpi=150)
    plt.close(fig)
    return lines


def gamma_beta_corrected(G_md, system, n=0.5, tol=1e-10, itmax=200):
    """Author's hypothesis (2026-10-07): the MD Γ (pure Al) ignores the solute redistribution β of the MRB.
    Candidate correction Γ_eff = Γ_MD · β̄(Γ_eff)^(−n), solved self-consistently point by point in t_SL
    (β̄ = weighted β of the MRB, which depends on SDAS through the Fourier numbers). Returns t, Γ_eff, SDAS, β̄."""
    from nucleation_md.sdas import diff_length_scale, sdas_mrb
    VL, t = system.VL(system.P), system.tSL(system.P)
    _, kef, w = diff_length_scale(system, VL)
    G = np.full_like(t, G_md, dtype=float)
    for _ in range(itmax):
        S, b = sdas_mrb(G, system, kef, w, t)
        bw = np.asarray(sum(bi * wi for bi, wi in zip(b, w)))
        Gn = G_md * bw ** (-n)
        if np.max(np.abs(Gn - G) / G) < tol:
            break
        G = Gn
    S, b = sdas_mrb(G, system, kef, w, t)
    return t, G, S, np.asarray(sum(bi * wi for bi, wi in zip(b, w)))


def fig_beta(G_md, n=1/3, n_fit=0.5):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    lines = []
    for ax, (alloy, (name, ref)) in zip(axes, ALLOYS.items()):
        src = ContinuumSource(alloy=alloy, dSv_mode="closure", theta_mode="homogeneous")
        s = src.state_at_gradient(src.grad_T_exp)
        fm = FerreiraModel()
        g1, g2 = fm.gamma_tensor(s, 1), fm.gamma_tensor(s, 2)
        system = ALLOY_SDAS[alloy]
        exp = pd.read_csv(ROOT / "data" / "experimental" / f"sdas_{alloy}.csv", comment="#")
        t, S1 = sdas_from_gamma(g1, system, "MRB")
        _, S2 = sdas_from_gamma(g2, system, "MRB")
        _, Sm = sdas_from_gamma(G_md, system, "MRB")
        _, R2 = sdas_from_gamma(g2, system, "RB")
        _, Rm = sdas_from_gamma(G_md, system, "RB")
        tb, Geff, Sb, bw = gamma_beta_corrected(G_md, system, n)
        tf, Gf, Sf, _ = gamma_beta_corrected(G_md, system, n_fit)
        win = (tb >= exp.tSL_s.min()) & (tb <= exp.tSL_s.max())
        ax.plot(t, R2 * 1e6, color=INK, lw=1.8, label=f"RB (equilíbrio, β = 1), contínuo Γ²ⁿᵈ")
        ax.plot(t, Rm * 1e6, color=INK, lw=1.8, ls=(0, (6, 3)), label=f"RB (equilíbrio, β = 1), MD Γ = {G_md:.2e} m·K")
        ax.plot(t, S2 * 1e6, color=C2, lw=2, label="contínuo Γ²ⁿᵈ, MRB")
        ax.plot(t, S1 * 1e6, color=C1, lw=2, label=f"contínuo Γ¹ˢᵗ = {g1:.2e} m·K, MRB")
        ax.plot(t, Sm * 1e6, color=C4, lw=1.2, ls=(0, (2, 2)), label="MD, MRB, Γ_MD sem β")
        ax.plot(tf, Sf * 1e6, color=C4, lw=1.8, ls=(0, (6, 3)), label="MD, MRB, Γ_MD·β̄^(−1/2) (expoente ajustado)")
        ax.plot(tb, Sb * 1e6, color=C4, lw=2.8,
                label=f"MD, MRB, Γ_MD·β̄^(−1/3) (teoria) = {Geff[win].min():.2e}–{Geff[win].max():.2e} m·K")
        ax.errorbar(exp.tSL_s, exp.SDAS_um, yerr=[exp.err_minus_um, exp.err_plus_um], fmt="s", color=INK, ms=6,
                    mfc=INK, mec="white", mew=1, elinewidth=1.2, capsize=0, label=f"experimental, {ref}", zorder=6)
        ax.set(xscale="log", yscale="log")
        ax.set_xlabel("tempo local de solidificação t_SL [s]", color=INK)
        ax.set_title(f"{name} — β̄ na janela exp.: {bw[win].min():.3f}–{bw[win].max():.3f}", color=INK, fontsize=11, loc="left")
        style(ax)
        ax.legend(fontsize=7.5, frameon=False, loc="upper left", labelcolor=INK)
        for lab, tt, S in (("MRB Γ¹ˢᵗ contínuo", t, S1), ("MRB MD sem β", t, Sm), ("MRB MD com β̄^(−1/3), teoria", tb, Sb),
                           ("MRB MD com β̄^(−1/2), ajuste", tf, Sf)):
            a, b = power_law_fit(tt, S)
            pred = a * exp.tSL_s ** b
            inside = ((pred <= exp.SDAS_um + exp.err_plus_um) & (pred >= exp.SDAS_um - exp.err_minus_um)).mean()
            rms = np.sqrt(np.mean(np.log(pred / exp.SDAS_um) ** 2)); bias = np.mean(np.log(pred / exp.SDAS_um))
            lines.append(f"| {name} | {lab} | {a:.2f}·t^{b:.3f} | {inside:.0%} | {rms:.3f} | {bias:+.3f} |")
    axes[0].set_ylabel("SDAS [μm]", color=INK)
    fig.suptitle("Redistribuição de soluto β: Γ da MD (Al puro, sem soluto) corrigido por β̄^(−1/3) da MRB (teoria)", color=INK, fontsize=11)
    fig.text(0.01, 0.005, "β = 1 alavanca (RB), β = 0 Scheil, 0 < β < 1 difusão finita (MRB). β̄ = β ponderado da MRB, "
             "autoconsistente com o SDAS corrigido. Pretas = RB: contínua = contínuo, tracejada = MD.", fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    fig.savefig(OUT / "fig_sdas_beta_md_continuo.png", dpi=150)
    plt.close(fig)
    return lines


def fig_mcmd(G_boro, G_nep_pure=8.582e-8, G_nep_alloy=8.247e-8):
    """SDAS with the MD Γ measured with infinite solute diffusion (md/e4_mcmd: NEP89, Al-3Cu, SGC MC, β = 1) against
    pure Al (NEP89 and Borovikov), the continuum (Ferreira 2024) and the experiments; RB black (solid = continuum,
    dashed = MD with β = 1)."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    lines = []
    for ax, (alloy, (name, ref)) in zip(axes, ALLOYS.items()):
        src = ContinuumSource(alloy=alloy, dSv_mode="closure", theta_mode="homogeneous")
        s = src.state_at_gradient(src.grad_T_exp)
        fm = FerreiraModel()
        g1, g2 = fm.gamma_tensor(s, 1), fm.gamma_tensor(s, 2)
        system = ALLOY_SDAS[alloy]
        exp = pd.read_csv(ROOT / "data" / "experimental" / f"sdas_{alloy}.csv", comment="#")
        t, S1 = sdas_from_gamma(g1, system, "MRB")
        _, S2 = sdas_from_gamma(g2, system, "MRB")
        _, R2 = sdas_from_gamma(g2, system, "RB")
        _, Sb = sdas_from_gamma(G_boro, system, "MRB")
        _, Sp = sdas_from_gamma(G_nep_pure, system, "MRB")
        _, Sa = sdas_from_gamma(G_nep_alloy, system, "MRB")
        _, Ra = sdas_from_gamma(G_nep_alloy, system, "RB")
        ax.plot(t, R2 * 1e6, color=INK, lw=1.8, label=f"RB (equilíbrio), contínuo Γ²ⁿᵈ = {g2:.2e} m·K")
        ax.plot(t, Ra * 1e6, color=INK, lw=1.8, ls=(0, (6, 3)), label=f"RB (equilíbrio), MD Al-3Cu β = 1, Γ = {G_nep_alloy:.2e} m·K")
        ax.plot(t, S2 * 1e6, color=C2, lw=2, label=f"contínuo Γ²ⁿᵈ, MRB")
        ax.plot(t, S1 * 1e6, color=C1, lw=2, label=f"contínuo Γ¹ˢᵗ = {g1:.2e} m·K, MRB")
        ax.plot(t, Sb * 1e6, color=C4, lw=1.6, ls=(0, (2, 2)), label=f"MD Al puro (Borovikov), Γ = {G_boro:.2e} m·K, MRB")
        ax.plot(t, Sp * 1e6, color=C3, lw=1.6, ls=(0, (2, 2)), label=f"MD Al puro (NEP89), Γ = {G_nep_pure:.2e} m·K, MRB")
        ax.plot(t, Sa * 1e6, color=C3, lw=2.6, label=f"MD Al-3Cu com MC, β = 1 (NEP89), Γ = {G_nep_alloy:.2e} m·K, MRB")
        ax.errorbar(exp.tSL_s, exp.SDAS_um, yerr=[exp.err_minus_um, exp.err_plus_um], fmt="s", color=INK, ms=6,
                    mfc=INK, mec="white", mew=1, elinewidth=1.2, capsize=0, label=f"experimental, {ref}", zorder=6)
        ax.set(xscale="log", yscale="log")
        ax.set_xlabel("tempo local de solidificação t_SL [s]", color=INK)
        ax.set_title(name, color=INK, fontsize=12, loc="left")
        style(ax)
        ax.legend(fontsize=7.5, frameon=False, loc="upper left", labelcolor=INK)
        for lab, S in (("MRB Γ¹ˢᵗ contínuo", S1), ("MRB MD Al puro Borovikov", Sb), ("MRB MD Al puro NEP89", Sp),
                       ("MRB MD Al-3Cu β = 1 NEP89", Sa), ("RB MD Al-3Cu β = 1 NEP89", Ra)):
            a, b = power_law_fit(t, S)
            pred = a * exp.tSL_s ** b
            inside = ((pred <= exp.SDAS_um + exp.err_plus_um) & (pred >= exp.SDAS_um - exp.err_minus_um)).mean()
            rms = np.sqrt(np.mean(np.log(pred / exp.SDAS_um) ** 2)); bias = np.mean(np.log(pred / exp.SDAS_um))
            lines.append(f"| {name} | {lab} | {a:.2f}·t^{b:.3f} | {inside:.0%} | {rms:.3f} | {bias:+.3f} |")
    axes[0].set_ylabel("SDAS [μm]", color=INK)
    fig.suptitle("SDAS: contínuo (Ferreira, 2024) × MD com difusão infinita do soluto (MD + Monte Carlo, β = 1)", color=INK, fontsize=11)
    fig.text(0.01, 0.005, "MD Al-3Cu: NEP89 + MC semi-grande canônico (Δμ = 1,013 eV, T_L = 914,3 K, k ≈ 0,41), seeding R = 2 nm "
             "(md/e4_mcmd). Pretas = RB: contínua = contínuo, tracejada = MD β = 1.", fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    fig.savefig(OUT / "fig_sdas_mcmd_md_continuo.png", dpi=150)
    plt.close(fig)
    return lines


def fig_md_dl(G_md, D0L_md=9.383e-8, QL_md=23855.0):
    """Al-3Cu-5Nb-0,1Fe: SDAS with everything the MD provides so far — Γ from seeding (Borovikov, ∇T = 0), D_L of Cu from
    NEP89 (md/e4_diffusion/nep89/dl_nep89.md) and the theoretical β̄^(−1/3) correction — against the continuum
    (Ferreira 2024, script data), RB and the experiments. D_S, k and the other solutes remain the script's."""
    import copy
    alloy, (name, ref) = "Al3Cu5Nb01Fe", ALLOYS["Al3Cu5Nb01Fe"]
    sy0 = ALLOY_SDAS[alloy]
    sy1 = copy.deepcopy(sy0)
    cu = [s for s in sy1.solutes if s.name == "Cu"][0]
    cu.D0L, cu.QL = D0L_md, QL_md
    src = ContinuumSource(alloy=alloy, dSv_mode="closure", theta_mode="homogeneous")
    st = src.state_at_gradient(src.grad_T_exp)
    g1, g2 = FerreiraModel().gamma_tensor(st, 1), FerreiraModel().gamma_tensor(st, 2)
    exp = pd.read_csv(ROOT / "data" / "experimental" / f"sdas_{alloy}.csv", comment="#")
    t, S1 = sdas_from_gamma(g1, sy0, "MRB")
    _, S2 = sdas_from_gamma(g2, sy0, "MRB")
    _, R2 = sdas_from_gamma(g2, sy0, "RB")
    _, Rm = sdas_from_gamma(G_md, sy1, "RB")
    _, Sm = sdas_from_gamma(G_md, sy1, "MRB")
    tb, Geff, Sb, bw = gamma_beta_corrected(G_md, sy1, 1 / 3)
    fig, ax = plt.subplots(figsize=(8, 5.4))
    ax.plot(t, R2 * 1e6, color=INK, lw=1.8, label=f"RB (equilíbrio), contínuo Γ²ⁿᵈ = {g2:.2e} m·K")
    ax.plot(t, Rm * 1e6, color=INK, lw=1.8, ls=(0, (6, 3)), label=f"RB (equilíbrio), MD: Γ = {G_md:.2e} m·K, D_L(Cu) da MD")
    ax.plot(t, S2 * 1e6, color=C2, lw=2, label="contínuo Γ²ⁿᵈ, MRB (Ferreira, 2024)")
    ax.plot(t, S1 * 1e6, color=C1, lw=2, label=f"contínuo Γ¹ˢᵗ = {g1:.2e} m·K, MRB (Ferreira, 2024)")
    ax.plot(t, Sm * 1e6, color=C4, lw=1.2, ls=(0, (2, 2)), label="MD: Γ e D_L(Cu) da MD, MRB, sem β̄^(−1/3)")
    ax.plot(tb, Sb * 1e6, color=C4, lw=2.8, label="MD: Γ e D_L(Cu) da MD, MRB com β̄^(−1/3) (teoria)")
    ax.errorbar(exp.tSL_s, exp.SDAS_um, yerr=[exp.err_minus_um, exp.err_plus_um], fmt="s", color=INK, ms=6,
                mfc=INK, mec="white", mew=1, elinewidth=1.2, capsize=0, label=f"experimental, {ref}", zorder=6)
    ax.set(xscale="log", yscale="log")
    ax.set_xlabel("tempo local de solidificação t_SL [s]", color=INK)
    ax.set_ylabel("SDAS [μm]", color=INK)
    ax.set_title(f"{name}: SDAS do contínuo × SDAS com entradas da dinâmica molecular", color=INK, fontsize=11, loc="left")
    style(ax)
    ax.legend(fontsize=8, frameon=False, loc="upper left", labelcolor=INK)
    fig.text(0.01, 0.005, f"MD: Γ por seeding (Borovikov, ∇T = 0); D_L(Cu) = {D0L_md:.2e}·exp(−{QL_md:.0f}/RT) m²/s (NEP89; script "
             "1,06e-7·exp(−24000/RT)).\nD_S, k e demais solutos: do script.", fontsize=7.5, color=INK2)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(OUT / "fig_sdas_md_entradas.png", dpi=150)
    plt.close(fig)
    lines = []
    for lab, tt, S in (("MRB Γ¹ˢᵗ contínuo", t, S1), ("MRB MD (Γ, D_L) sem β", t, Sm), ("MRB MD (Γ, D_L) com β̄^(−1/3)", tb, Sb),
                       ("RB MD (Γ, D_L)", t, Rm), ("RB Γ²ⁿᵈ contínuo", t, R2)):
        a, b = power_law_fit(tt, S)
        p = a * exp.tSL_s ** b
        inside = ((p <= exp.SDAS_um + exp.err_plus_um) & (p >= exp.SDAS_um - exp.err_minus_um)).mean()
        lines.append(f"| {name} | {lab} | {a:.2f}·t^{b:.3f} | {inside:.0%} | {np.sqrt(np.mean(np.log(p/exp.SDAS_um)**2)):.3f} | "
                     f"{np.mean(np.log(p/exp.SDAS_um)):+.3f} |")
    return lines


def fig_isolado(G_md, D0L_md=9.383e-8, QL_md=23855.0):
    """Only: experiments, the author's original RB (Γ²ⁿᵈ), MRB 1st and 2nd order (Ferreira 2024) and the MD curve
    equivalent to the 2nd order (Γ = ΔT·r_c/2 with the critical radius measured by seeding; MRB with β̄^(−1/3);
    D_L of Cu from NEP89 where the alloy has Cu)."""
    import copy
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    lines = []
    for ax, (alloy, (name, ref)) in zip(axes, ALLOYS.items()):
        sy0 = ALLOY_SDAS[alloy]
        sy1 = copy.deepcopy(sy0)
        for s_ in sy1.solutes:
            if s_.name == "Cu":
                s_.D0L, s_.QL = D0L_md, QL_md
        src = ContinuumSource(alloy=alloy, dSv_mode="closure", theta_mode="homogeneous")
        st = src.state_at_gradient(src.grad_T_exp)
        g1, g2 = FerreiraModel().gamma_tensor(st, 1), FerreiraModel().gamma_tensor(st, 2)
        exp = pd.read_csv(ROOT / "data" / "experimental" / f"sdas_{alloy}.csv", comment="#")
        t, S1 = sdas_from_gamma(g1, sy0, "MRB")
        _, S2 = sdas_from_gamma(g2, sy0, "MRB")
        _, R2 = sdas_from_gamma(g2, sy0, "RB")
        tb, Geff, Sb, _ = gamma_beta_corrected(G_md, sy1, 1 / 3)
        ax.plot(t, R2 * 1e6, color=INK, lw=1.8, label="RB, Rappaz e Boettinger (1999), Γ²ⁿᵈ")
        ax.plot(t, S1 * 1e6, color=C1, lw=2, label=f"MRB, Ferreira (2024), Γ¹ˢᵗ = {g1:.2e} m·K")
        ax.plot(t, S2 * 1e6, color=C2, lw=2, label=f"MRB, Ferreira (2024), Γ²ⁿᵈ = {g2:.2e} m·K")
        md_lab = "MD equivalente de 2ª ordem: Γ = ΔT·r_c/2 (seeding)" + (", D_L(Cu) da MD" if any(s_.name == "Cu" for s_ in sy1.solutes) else "")
        ax.plot(tb, Sb * 1e6, color=C4, lw=2.8, label=md_lab + ", β̄^(−1/3)")
        ax.errorbar(exp.tSL_s, exp.SDAS_um, yerr=[exp.err_minus_um, exp.err_plus_um], fmt="s", color=INK, ms=6,
                    mfc=INK, mec="white", mew=1, elinewidth=1.2, capsize=0, label=f"experimental, {ref}", zorder=6)
        ax.set(xscale="log", yscale="log")
        ax.set_xlabel("tempo local de solidificação t_SL [s]", color=INK)
        ax.set_title(name, color=INK, fontsize=12, loc="left")
        style(ax)
        ax.legend(fontsize=7.5, frameon=False, loc="upper left", labelcolor=INK)
        for lab, tt, S in (("RB Γ²ⁿᵈ", t, R2), ("MRB Γ¹ˢᵗ", t, S1), ("MRB Γ²ⁿᵈ", t, S2), ("MD equiv. 2ª ordem", tb, Sb)):
            a, b = power_law_fit(tt, S)
            p = a * exp.tSL_s ** b
            inside = ((p <= exp.SDAS_um + exp.err_plus_um) & (p >= exp.SDAS_um - exp.err_minus_um)).mean()
            lines.append(f"| {name} | {lab} | {a:.2f}·t^{b:.3f} | {inside:.0%} | {np.mean(np.log(p/exp.SDAS_um)):+.3f} |")
    axes[0].set_ylabel("SDAS [μm]", color=INK)
    fig.suptitle("SDAS: experimental × RB × MRB de 1ª e 2ª ordem (Ferreira, 2024) × dinâmica molecular", color=INK, fontsize=11)
    fig.text(0.01, 0.005, f"MD: Γ = {G_md:.2e} m·K (seeding, Al puro, Borovikov), corrigido por β̄^(−1/3) da MRB; "
             "D_L(Cu) do NEP89 na liga com Cu; demais dados de difusão e partição: dos scripts.", fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    fig.savefig(OUT / "fig_sdas_isolado.png", dpi=150)
    plt.close(fig)
    return lines


def fig_isolado_en(G_md, D0L_md=9.383e-8, QL_md=23855.0):
    """English version of fig_isolado (same curves and data)."""
    import copy
    names = {"Al3Cu5Nb01Fe": ("Al-3Cu-5Nb-0.1Fe", "Mendes et al. (2023)"),
             "Al08Si06Mg02Fe": ("Al-0.8Si-0.6Mg-0.2Fe", "Marques et al. (2025)")}
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    for ax, (alloy, (name, ref)) in zip(axes, names.items()):
        sy0 = ALLOY_SDAS[alloy]
        sy1 = copy.deepcopy(sy0)
        has_cu = False
        for s_ in sy1.solutes:
            if s_.name == "Cu":
                s_.D0L, s_.QL = D0L_md, QL_md
                has_cu = True
        src = ContinuumSource(alloy=alloy, dSv_mode="closure", theta_mode="homogeneous")
        st = src.state_at_gradient(src.grad_T_exp)
        g1, g2 = FerreiraModel().gamma_tensor(st, 1), FerreiraModel().gamma_tensor(st, 2)
        exp = pd.read_csv(ROOT / "data" / "experimental" / f"sdas_{alloy}.csv", comment="#")
        t, S1 = sdas_from_gamma(g1, sy0, "MRB")
        _, S2 = sdas_from_gamma(g2, sy0, "MRB")
        _, R2 = sdas_from_gamma(g2, sy0, "RB")
        tb, _, Sb, _ = gamma_beta_corrected(G_md, sy1, 1 / 3)
        ax.plot(t, R2 * 1e6, color=INK, lw=1.8, label="RB, Rappaz and Boettinger (1999), Γ²ⁿᵈ")
        ax.plot(t, S1 * 1e6, color=C1, lw=2, label=f"MRB, Ferreira (2024), Γ¹ˢᵗ = {g1:.2e} m·K")
        ax.plot(t, S2 * 1e6, color=C2, lw=2, label=f"MRB, Ferreira (2024), Γ²ⁿᵈ = {g2:.2e} m·K")
        ax.plot(tb, Sb * 1e6, color=C4, lw=2.8,
                label="Molecular dynamics, 2nd-order equivalent: Γ = ΔT·r_c/2 (seeding)"
                      + (", MD D_L(Cu)" if has_cu else "") + ", β̄^(−1/3)")
        ax.errorbar(exp.tSL_s, exp.SDAS_um, yerr=[exp.err_minus_um, exp.err_plus_um], fmt="s", color=INK, ms=6,
                    mfc=INK, mec="white", mew=1, elinewidth=1.2, capsize=0, label=f"Experimental, {ref}", zorder=6)
        ax.set(xscale="log", yscale="log")
        ax.set_xlabel("local solidification time t_SL [s]", color=INK)
        ax.set_title(name, color=INK, fontsize=12, loc="left")
        style(ax)
        ax.legend(fontsize=7.5, frameon=False, loc="upper left", labelcolor=INK)
    axes[0].set_ylabel("SDAS [μm]", color=INK)
    fig.suptitle("SDAS: experiments × RB × 1st- and 2nd-order MRB (Ferreira, 2024) × molecular dynamics", color=INK, fontsize=11)
    fig.text(0.01, 0.005, f"MD: Γ = {G_md:.2e} m·K (seeding, pure Al, Borovikov 2024 EAM), corrected by the MRB β̄^(−1/3); "
             "D_L(Cu) from NEP89 in the Cu-bearing alloy; other diffusion and partition data from the scripts.", fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    fig.savefig(OUT / "fig_sdas_isolated_en.png", dpi=200)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    md = md_seeding()
    G_md = md.Gamma.mean()
    # all MD estimates of Γ: seeding at ∇T = 0 (4 radii) and NEMD T*(G) (md/e4_nemd/nemd_summary.md)
    G_all = list(md.Gamma)
    f = ROOT / "md" / "e4_nemd" / "nemd_summary.md"
    if f.exists():
        G_all += [float(x) * 1e-7 for x in re.findall(r"^\| [\d.e+]+ \| \d+–\d+ \| [\d.]+ \| [\d.]+ \| [\d.]+ \| ([\d.]+) \|",
                                                       f.read_text(encoding="utf-8"), re.M)]
    for l in fig_sdas(G_md) + fig_gamma(md):
        print(l)
    print(f"\nΓ_MD: média {G_md:.3e}, faixa {min(G_all):.3e}–{max(G_all):.3e} m·K ({len(G_all)} estimativas)\n")
    print("| liga | curva | SDAS [μm] | dentro das barras | rms_log | bias_log |\n|---|---|---|---|---|---|")
    for l in fig_confronto(G_md, min(G_all), max(G_all)):
        print(l)
    print()
    for l in fig_gap(G_md, min(G_all), max(G_all)):
        print(l)
    print("\n| liga | curva | SDAS [μm] | dentro das barras | rms_log | bias_log |\n|---|---|---|---|---|---|")
    for l in fig_gap_sdas(G_md):
        print(l)
    print()
    for l in fig_beta(G_md):
        print(l)
    print()
    for l in fig_mcmd(G_md):
        print(l)
    print()
    for l in fig_md_dl(G_md):
        print(l)
    print()
    for l in fig_isolado(G_md):
        print(l)
    fig_isolado_en(G_md)


if __name__ == "__main__":
    main()
