# -*- coding: utf-8 -*-
"""
Nucleation models with a common interface.

    delta_G(r, state)            ΔG(r)                                  [J]
    r_critical(state, order)     critical radius                        [m]
    barrier(state)               ΔG(r_c) = ΔG*                          [J]
    gamma_tensor(state, order)   thermal-field tensor Γ = A·∇T          [K m]
    rate(state, D, lam, N_V)     nucleation rate I                      [m^-3 s^-1]

Notation in code (both are written ΔG_C in the paper):
    dG_conf  - configurational free energy in the denominator of r_c,  γ_SL·(∂f/∂r)/f
    dG_star  - critical free energy ΔG(r_c) used in the rate
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .geometry import cap_area

K_B = 1.380649e-23  # J/K


class NucleationModel(ABC):

    @abstractmethod
    def delta_G(self, r, state):
        ...

    @abstractmethod
    def r_critical(self, state, order=1):
        ...

    def barrier(self, state, order=1):
        return self.delta_G(self.r_critical(state, order), state)

    @abstractmethod
    def gamma_tensor(self, state, order=1):
        ...

    def grad_T(self, state, order=None):
        """∇T = Γ / A, A = solid–liquid area of the cap at r_c (order None = model default)."""
        order = getattr(self, "order", 1) if order is None else order
        rc = self.r_critical(state, order)
        return self.gamma_tensor(state, order) / cap_area(rc, state.theta(rc))

    @abstractmethod
    def rate(self, state, D, lam, N_V, **kwargs):
        ...


# --------------------------------------------------------------------------- CNT

class CNT(NucleationModel):
    """
    Classical nucleation theory: ΔS_V, ΔT, γ_SL, θ frozen at state.r_eval.

        ΔG(r) = [(1/3)πr³ΔS_VΔT + πr²γ_SL]·f(θ)
        r_c   = -2γ_SL / (ΔS_VΔT) = -2ΔG_S/ΔG_V
        ΔG*   = 4πγ³f / (3(ΔS_VΔT)²)
        Γ     = -γ_SL/ΔS_V
        I     = (D·A/λ⁴)·(N/V)·exp(-ΔG*/(k_B·T))
    """

    def _frozen(self, state):
        r = state.r_eval
        return state.DSv(r), state.DT(r), state.gamma(r), state.f(r)

    def _gamma_at(self, r, state):
        return self._frozen(state)[2]

    def delta_G(self, r, state):
        DSv, DT, _, f = self._frozen(state)
        g = self._gamma_at(r, state)
        return (np.pi * r ** 3 * DSv * DT / 3.0 + np.pi * r ** 2 * g) * f

    def r_critical(self, state, order=1):
        if order != 1:
            raise ValueError("CNT has a single critical radius (order=1)")
        DSv, DT, g, _ = self._frozen(state)
        return -2.0 * g / (DSv * DT)

    def gamma_tensor(self, state, order=1):
        DSv, _, g, _ = self._frozen(state)
        return -g / DSv

    def rate(self, state, D, lam, N_V, T=None, **kwargs):
        rc = self.r_critical(state)
        A = cap_area(rc, state.theta(state.r_eval))
        T = state.temperature() if T is None else T
        return D * A / lam ** 4 * N_V * np.exp(-self.barrier(state) / (K_B * T))


class Tolman(CNT):
    """
    CNT with Tolman curvature correction γ(r) = γ_∞/(1 + 2δ/r), γ_∞ = γ_SL(r_eval).
    δ (Tolman length) has no default: it must come from data.

    Usual use in nucleation: CNT with the effective γ_eff = γ(r_c), self-consistent with
    Gibbs–Thomson/Laplace, which gives
        r_c = r_c,CNT - 2δ,   ΔG(r) = CNT form with γ_eff  (so dΔG/dr = 0 at r_c).
    """

    def __init__(self, delta):
        self.delta = float(delta)

    def _rc(self, state):
        return super().r_critical(state) - 2.0 * self.delta

    def _gamma_at(self, r, state):
        return self._frozen(state)[2] / (1.0 + 2.0 * self.delta / self._rc(state))

    def r_critical(self, state, order=1):
        if order != 1:
            raise ValueError("Tolman has a single critical radius (order=1)")
        return self._rc(state)

    def gamma_tensor(self, state, order=1):
        return -self._gamma_at(None, state) / self._frozen(state)[0]


# ------------------------------------------------------------------ Ferreira (2024)

@dataclass
class Coefficients:
    """Terms of a·r² + 3b·r + 6γ_SL·f = 0 evaluated at radius r."""
    r: float
    a: float
    b: float
    c: float          # 6·γ_SL·f
    f: float
    df_dr: float
    dG_V: float       # ΔS_V·ΔT
    dG_S: float       # γ_SL
    ddG_S_dr: float   # ∂γ_SL/∂r
    dG_conf: float    # γ_SL·(∂f/∂r)/f
    DSv: float
    dDSv_dr: float
    DT: float
    dDT_dr: float

    @property
    def denominator(self):
        """ΔG_V + ∂ΔG_S/∂r + ΔG_C  (= b/f)."""
        return self.dG_V + self.ddG_S_dr + self.dG_conf


class FerreiraModel(NucleationModel):
    """
    Ferreira (2024), Physica B, doi:10.1016/j.physb.2024.416494.

        ΔG(r) = [(1/3)πr³ΔS_VΔT + πr²γ_SL]·f(θ),  ΔS_V, ΔT, γ_SL, θ functions of r
        dΔG/dr = (πr/3)·(a·r² + 3b·r + 6γ_SL·f)
        a = (∂ΔS_V/∂r·ΔT + ΔS_V·∂ΔT/∂r)·f + ΔS_V·ΔT·∂f/∂r
        b = (ΔS_V·ΔT + ∂γ_SL/∂r)·f + γ_SL·∂f/∂r
        r_c¹ = -2γ_SL·f/b = -2ΔG_S/(ΔG_V + ∂ΔG_S/∂r + ΔG_C)
        r_c² = -(3/2)·b/a
        Γ¹ˢᵗ = -ΔG_S/(ΔS_V + ∂ΔS_S/∂r + ΔS_C) = ΔT·r_c¹/2       (entropies S ≡ G/ΔT)
        Γ²ⁿᵈ = -(3/4)·[ΔG_V + ∂ΔG_S/∂r + ΔG_C]/(∂ΔS_V/∂r) = ΔT·r_c²/2
        Γ    = A·∇T   ->   ∇T = Γ/A  (only the thermal gradient displaces the equilibrium)
        I    = (D·A/λ⁴)·(N/V)·exp(ΔG*/ΔG*_Eq)

    reference      - State of the displaced-equilibrium reference (state 0); ΔG*_Eq = barrier(reference)
    dG_eq          - alternatively ΔG*_Eq given directly [J]
    order          - critical radius used by barrier() and rate()
    exponent_sign  - +1 reproduces exp(ΔG*/ΔG*_Eq) as written; kept as a parameter until the
                     sign convention is confirmed (docs/CONVENTIONS.md)
    """

    def __init__(self, reference=None, dG_eq=None, order=2, exponent_sign=+1):
        if order not in (1, 2):
            raise ValueError("order must be 1 or 2")
        self.order = order
        self.exponent_sign = exponent_sign
        self.reference = reference
        if dG_eq is not None:
            self.dG_eq = float(dG_eq)
        elif reference is not None:
            self.dG_eq = self.barrier(reference)
        else:
            self.dG_eq = None

    # -- coefficients ---------------------------------------------------------

    def coefficients(self, state, r=None):
        r = state.r_eval if r is None else r
        DSv, dDSv = state.DSv(r), state.DSv.d(r)
        DT, dDT = state.DT(r), state.DT.d(r)
        g, dg = state.gamma(r), state.gamma.d(r)
        f, df = state.f(r), state.df_dr(r)
        a = (dDSv * DT + DSv * dDT) * f + DSv * DT * df
        b = (DSv * DT + dg) * f + g * df
        return Coefficients(r=r, a=a, b=b, c=6.0 * g * f, f=f, df_dr=df,
                            dG_V=DSv * DT, dG_S=g, ddG_S_dr=dg, dG_conf=g * df / f,
                            DSv=DSv, dDSv_dr=dDSv, DT=DT, dDT_dr=dDT)

    # -- interface ------------------------------------------------------------

    def delta_G(self, r, state):
        return (np.pi * r ** 3 * state.DSv(r) * state.DT(r) / 3.0
                + np.pi * r ** 2 * state.gamma(r)) * state.f(r)

    def dG_dr(self, r, state):
        """Analytic dΔG/dr = (πr/3)·(a r² + 3b r + 6γf), coefficients taken at r."""
        k = self.coefficients(state, r)
        return np.pi * r / 3.0 * (k.a * r ** 2 + 3.0 * k.b * r + k.c)

    def r_critical(self, state, order=None, solve="explicit", bracket=None):
        """
        solve="explicit"        coefficients at state.r_eval (the formulas as written)
        solve="self_consistent" r such that r = r_c(r); needs bracket=(r_lo, r_hi)
        """
        order = self.order if order is None else order
        if solve == "explicit":
            return self._rc_formula(self.coefficients(state), order)
        if solve == "self_consistent":
            if bracket is None:
                raise ValueError("self_consistent solve needs bracket=(r_lo, r_hi)")
            return brentq(lambda r: r - self._rc_formula(self.coefficients(state, r), order),
                          *bracket, xtol=1e-18, rtol=1e-12)
        raise ValueError(f"unknown solve mode {solve!r}")

    @staticmethod
    def _rc_formula(k, order):
        if order == 1:
            return -2.0 * k.dG_S * k.f / k.b
        if order == 2:
            return -1.5 * k.b / k.a if k.a != 0.0 else np.inf
        raise ValueError("order must be 1 or 2")

    def exact_roots(self, state, r=None):
        """
        Roots of a·r² + 3b·r + 6γf = 0 with the coefficients frozen at r (default r_eval):
            r± = [-3b ± √(9b² - 24·a·γ·f)] / (2a)
        Returns (r_minus, r_plus, discriminant); roots are complex when the discriminant < 0.
        Reported for comparison only: r_c² stays defined as -3b/(2a).
        """
        k = self.coefficients(state, r)
        disc = 9.0 * k.b ** 2 - 4.0 * k.a * k.c
        sq = np.lib.scimath.sqrt(disc)
        return (-3.0 * k.b - sq) / (2.0 * k.a), (-3.0 * k.b + sq) / (2.0 * k.a), disc

    def barrier(self, state, order=None):
        return self.delta_G(self.r_critical(state, order), state)

    def gamma_tensor(self, state, order=None):
        """
        Thermal-field tensor (MRS 2026, slide 10), entropies defined as S ≡ G/ΔT:
            Γ¹ˢᵗ = -ΔG_S / (ΔS_V + ∂ΔS_S/∂r + ΔS_C)            = -ΔT·γf/b      = ΔT·r_c¹/2
            Γ²ⁿᵈ = -(3/4)[ΔG_V + ∂ΔG_S/∂r + ΔG_C] / (∂ΔS_V/∂r)  = -(3/4)·ΔT·b/a = ΔT·r_c²/2
        where ∂ΔS_V/∂r is the total derivative (1/ΔT)·∂(ΔS_VΔT f)/∂r / f = a/(f·ΔT).
        """
        order = self.order if order is None else order
        k = self.coefficients(state)
        if order == 1:
            return -k.dG_S / (k.denominator / k.DT)
        if order == 2:
            if k.a == 0.0:                  # parabola degenerates: Γ²ⁿᵈ undefined
                return np.inf
            return -0.75 * k.denominator / (k.a / (k.f * k.DT))
        raise ValueError("order must be 1 or 2")

    def rate(self, state, D, lam, N_V, order=None, **kwargs):
        if self.dG_eq is None:
            raise ValueError("FerreiraModel needs a reference state or dG_eq for the rate")
        rc = self.r_critical(state, order)
        A = cap_area(rc, state.theta(rc))
        dG_star = self.delta_G(rc, state)
        return D * A / lam ** 4 * N_V * np.exp(self.exponent_sign * dG_star / self.dG_eq)
