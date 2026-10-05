# -*- coding: utf-8 -*-
"""
Surface free energy ΔG_S(r) = γ_SL(r) providers. Each returns a Quantity.

metric_2sphere      - Ferreira, metric / Ricci-scalar form on the 2-sphere (R = 2/r₀²),
                      transcribed literally:
    ΔG_S = ΔG_S0·3/(1 + r²/r₀² + r²sin²θ/(r₀²sin²θ₀)) - (1/(4πr₀²))·∫_{A0}^{A1} Σ dA

gurtin_murdoch_legacy - the form the continuum script actually evaluates (L329, L355, L373):
    γ = γ₀·(r₀/r)² - Σ_GM(r),  Σ_GM = -2σ₀·ln(y₀/y),  y = (3λ+2μ)·r + 2(λ₀+μ₀)

Open points (docs/CONVENTIONS.md): which form is canonical; whether θ in the metric
term is the wetting angle θ(r); A0, A1 limits of the Σ integral.
"""
import numpy as np
from scipy.integrate import quad

from .quantity import Quantity


def metric_2sphere(dG_S0, r0, theta0=None, theta=None, Sigma=None, area=None):
    """
    dG_S0   - ΔG_S0 = γ₀ of the reference state [J m^-2]
    r0      - reference radius r₀ [m]
    theta0  - θ₀ [rad]; with theta=None the angular term is taken at θ = θ₀ (ratio = 1)
    theta   - Quantity θ(r) or None
    Sigma   - callable Σ(r) [N m^-1] (surface-stress trace) or None (no integral term)
    area    - callable A(r); default sphere 4πr². The integral runs from A0 = A(r₀) to A1 = A(r)
    """
    if theta is not None and (theta0 is None or np.isclose(np.sin(theta0), 0.0)):
        raise ValueError("theta(r) given: theta0 must be set and sin(theta0) != 0")
    area = (lambda r: 4.0 * np.pi * r ** 2) if area is None else area

    def angular(r):
        if theta is None:
            return 1.0
        return np.sin(theta(r)) ** 2 / np.sin(theta0) ** 2

    def stress_term(r):
        if Sigma is None:
            return 0.0
        # ∫_{A(r0)}^{A(r)} Σ dA = ∫_{r0}^{r} Σ(r')·dA/dr' dr'
        dA = lambda x: (area(x * (1 + 1e-7)) - area(x * (1 - 1e-7))) / (2e-7 * x)
        val, _ = quad(lambda x: Sigma(x) * dA(x), r0, r, limit=200)
        return val / (4.0 * np.pi * r0 ** 2)

    def g(r):
        r = np.asarray(r, dtype=float)
        if r.ndim:
            return np.array([g(x) for x in r])
        x2 = (r / r0) ** 2
        return dG_S0 * 3.0 / (1.0 + x2 + x2 * angular(r)) - stress_term(r)

    return Quantity(g, None, "gamma_metric_2sphere", rel_step=1e-4)


# ------------------------------------------------------- continuum (legacy) form

def gm_surface_stress(r, r0, sigma0, lam, mu, lam0, mu0):
    """`sigma_func` (L329): returns (tau, Σ_GM)."""
    y = (3 * lam + 2 * mu) * r + 2 * (lam0 + mu0)
    y0 = (3 * lam + 2 * mu) * r0 + 2 * (lam0 + mu0)
    surf_stress = -2 * sigma0 * np.log(y0 / y)
    return -sigma0 + surf_stress, surf_stress


def gurtin_murdoch_legacy(gamma0, r0, sigma0, lam, mu, lam0, mu0):
    """
    γ(r) = γ₀(r₀/r)² - Σ_GM(r)   (`gamma_func_r`, L355)
    dγ/dr from `dgammadr_func_r` (L373), as coded:
        -2γ₀r₀²/r³ - 2σ₀/(r(1+α)),  α = 2(σ₀ + 2(λ₀+μ₀))/(r(3λ+2μ))
    """
    def g(r):
        return gamma0 / (r / r0) ** 2 - gm_surface_stress(r, r0, sigma0, lam, mu, lam0, mu0)[1]

    def dg(r):
        alpha = 2 * (sigma0 + 2 * (lam0 + mu0)) / (r * (3 * lam + 2 * mu))
        return -2 * gamma0 * r0 ** 2 / r ** 3 - 2 * sigma0 / (r * (1 + alpha))

    return Quantity(g, dg, "gamma_gurtin_murdoch_legacy")
