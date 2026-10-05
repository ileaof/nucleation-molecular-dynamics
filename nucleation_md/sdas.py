# -*- coding: utf-8 -*-
"""
Secondary dendrite arm spacing from the thermal-field tensor Γ (MRS 2026, slide 18):

    SDAS = 5.5·(M·t_SL)^(1/3)
    M    = -Γ·Σ w_j β_j / Σ m_j(1 - k_ef,j)(c_f,j - c_0,j)/D_j · ln[ Σ m_j(1-k_ef,j)c_f,j/D_j / Σ m_j(1-k_ef,j)c_0,j/D_j ]

Literal port of `diff_length_scale`, `calc` (MRB, Ferreira 2024) and `calc_RB`
(Rappaz & Boettinger 1999) of the SDAS scripts:
    scriptsdasmodelAl08Si06Mg02Fe_zoom_horizontal_new_script_Ferreira_2024.py
    scriptsdasmodelAl3Cu5Nb01Fe_zoom_horizontal_B_MRS_Meeting_2026.py
Both scripts share the same functions; only the data differ (ALLOY_SDAS below).
"""
from dataclasses import dataclass, field

import numpy as np

R_GAS = 8.31451


def diffusion(T, D0, Q, R=R_GAS):
    return D0 * np.exp(-Q / (R * T))


@dataclass
class Solute:
    name: str
    c0: float      # wt%
    cf: float      # wt%
    k0: float
    mL: float      # K/wt%
    D0L: float
    QL: float
    D0S: float
    QS: float


@dataclass
class SDASSystem:
    TL: float
    Ts: float
    solutes: list
    VL: callable           # P [mm] -> V_L [mm/s]
    tSL: callable          # P [mm] -> t_SL [s]
    P: np.ndarray = field(default_factory=lambda: np.arange(1, 101, 1))
    R: float = R_GAS

    def DL(self):
        return [diffusion(self.TL, s.D0L, s.QL, self.R) for s in self.solutes]

    def DS(self):
        return [diffusion(self.Ts, s.D0S, s.QS, self.R) for s in self.solutes]


def diff_length_scale(system, VL):
    """`diff_length_scale`: mean diffusion lengths, effective k and solute weights."""
    V = VL / 1000.0
    DL = system.DL()
    deltam, kef = [], []
    for s, D in zip(system.solutes, DL):
        delta = D / V
        a = np.trapezoid(V, delta)                          # trapz(V, delta) as in the script
        dm = abs((1.0 / (V.max() - V.min())) * a)
        deltam.append(dm)
        kef.append(s.k0 / (s.k0 + (1.0 - s.k0) * np.exp(-V * dm / D)))
    c = np.array([s.c0 / 100 for s in system.solutes])
    w = c / c.sum()
    return deltam, kef, w


def sdas_mrb(Gamma, system, kef, w, tSL, tol=1e-10, legacy_beta3=True):
    """
    `calc`: fixed-point iteration on SDAS. Returns (SDAS [m], betas).
    legacy_beta3 - the script computes β₃ = 2γ₃Fo₂/(2Fo₃ + γ₃) (Fo of solute 2 in the
                   numerator); kept by default to reproduce the published curves.
    """
    DL, DS = system.DL(), system.DS()
    SDAS2, err = 1.0, 1.0
    while np.all(err > tol):
        SDAS = SDAS2
        Fo = [D * tSL / ((SDAS / 2) ** 2.0) for D in DS]
        gama = [(4.0 * F * s.k0) / (4.0 * F * s.k0 + 1.0) for F, s in zip(Fo, system.solutes)]
        beta = [(2.0 * g * F) / (2.0 * F + g) for g, F in zip(gama, Fo)]
        if legacy_beta3 and len(beta) >= 3:
            beta[2] = (2.0 * gama[2] * Fo[1]) / (2.0 * Fo[2] + gama[2])
        bw = sum(b * wj for b, wj in zip(beta, w))
        p1 = sum(s.mL * (1.0 - k) * (s.cf - s.c0) / D for s, k, D in zip(system.solutes, kef, DL))
        s2 = sum(s.mL * (1.0 - k) * s.cf / D for s, k, D in zip(system.solutes, kef, DL))
        t3 = sum(s.mL * (1.0 - k) * s.c0 / D for s, k, D in zip(system.solutes, kef, DL))
        M = bw * ((-Gamma / p1) * np.log(s2 / t3))
        SDAS2 = 5.5 * M ** (1 / 3) * tSL ** (1 / 3)
        err = np.abs(SDAS - SDAS2) / SDAS
    return SDAS, beta


def sdas_rb(Gamma, system, tSL):
    """`calc_RB`: Rappaz & Boettinger with equilibrium k0."""
    DL = system.DL()
    p1 = sum(s.mL * (1.0 - s.k0) * (s.cf - s.c0) / D for s, D in zip(system.solutes, DL))
    s2 = sum(s.mL * (1.0 - s.k0) * s.cf / D for s, D in zip(system.solutes, DL))
    t3 = sum(s.mL * (1.0 - s.k0) * s.c0 / D for s, D in zip(system.solutes, DL))
    return 5.5 * ((-Gamma / p1) * np.log(s2 / t3)) ** (1 / 3) * tSL ** (1 / 3)


def sdas_from_gamma(Gamma, system, model="MRB"):
    """SDAS(t_SL) [m] on the experimental kinetics of `system` for a given Γ [m K]."""
    VL, tSL = system.VL(system.P), system.tSL(system.P)
    if model == "RB":
        return tSL, sdas_rb(Gamma, system, tSL)
    _, kef, w = diff_length_scale(system, VL)
    return tSL, sdas_mrb(Gamma, system, kef, w, tSL)[0]


def power_law_fit(tSL, SDAS_m):
    """SDAS [μm] = a·t_SL^b, least squares in log-log as in the scripts."""
    b, ln_a = np.polyfit(np.log(tSL), np.log(SDAS_m * 1e6), 1)
    return np.exp(ln_a), b


# --------------------------------------------------------------- alloy data (copied)

ALLOY_SDAS = {
    # scriptsdasmodelAl08Si06Mg02Fe_zoom_horizontal_new_script_Ferreira_2024.py
    # (solute 1 is Mg, stored in the script under the *_Cu names)
    "Al08Si06Mg02Fe": SDASSystem(
        TL=925.32, Ts=525.0 + 273.15,
        solutes=[Solute("Mg", 0.6, 34.2, 0.36, -5.1, 9.9e-5, 71600.0, 6.23e-6, 115000.0),
                 Solute("Si", 0.8, 12.7, 0.11, -6.2, 1.34e-7, 30000.0, 2.48e-4, 137000.0),
                 Solute("Fe", 0.2, 1.81, 0.03, -4.2, 2.34e-7, 35000.0, 5.30e-3, 183400.0)],
        VL=lambda P: 1.2 * P ** (-0.27), tSL=lambda P: 3.58 * P ** 1.09,
        P=np.arange(1, 101, 1)),
    # scriptsdasmodelAl3Cu5Nb01Fe_zoom_horizontal_B_MRS_Meeting_2026.py
    "Al3Cu5Nb01Fe": SDASSystem(
        TL=936.58, Ts=817.75,
        solutes=[Solute("Cu", 3.0, 28.816, 0.1308, -3.64103, 1.06e-7, 24000.0, 6.5e-5, 136100.0),
                 Solute("Nb", 5.0, 11.02, 2.049, -2.96, 1.1e-5, 28400.0, 2.0e-4, 134500.0),
                 Solute("Fe", 0.2, 1.81, 0.03, -4.2, 2.34e-7, 35000.0, 5.30e-3, 183400.0)],
        VL=lambda P: 3.4 * P ** (-0.56), tSL=lambda P: 4.93116 * P ** 0.87606,
        P=np.arange(1, 251, 1)),
}
