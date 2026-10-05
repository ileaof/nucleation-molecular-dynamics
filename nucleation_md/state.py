# -*- coding: utf-8 -*-
"""
Thermodynamic state seen by a nucleation model.

Sign convention (see docs/CONVENTIONS.md):
    ΔS_V < 0   [J m^-3 K^-1]  solidification, ΔS_V = -ΔH_m·ρ_s/T_m in the limit
    ΔT   > 0   [K]            undercooling, T = T_m - ΔT
    γ_SL > 0   [J m^-2]       ΔG_S
    θ in (0, π] [rad]         wetting angle, θ = π homogeneous
"""
from dataclasses import dataclass
from typing import Optional

import numpy as np

from .geometry import f_theta, df_dtheta, df_legacy
from .quantity import Quantity

DF_MODES = ("chain", "continuum_legacy")


@dataclass
class State:
    """
    DSv, DT, gamma, theta - Quantity fields of r
    r_eval   - radius [m] where the coefficients are evaluated by the explicit
               critical-radius formulas and where CNT freezes the fields
    T_m      - melting / liquidus temperature [K]; needed only for kT-based rates
    df_mode  - "chain":            ∂f/∂r = f'(θ)·∂θ/∂r  [1/m]  (follows from ΔG(r))
               "continuum_legacy": ∂f/∂r = dfthetadr(θ) of the continuum script
    """
    DSv: Quantity
    DT: Quantity
    gamma: Quantity
    theta: Quantity
    r_eval: float
    T_m: Optional[float] = None
    df_mode: str = "chain"
    label: str = ""

    def __post_init__(self):
        if self.df_mode not in DF_MODES:
            raise ValueError(f"df_mode must be one of {DF_MODES}, got {self.df_mode!r}")
        if not self.r_eval > 0.0:
            raise ValueError("r_eval must be > 0")

    def f(self, r):
        return f_theta(self.theta(r))

    def df_dr(self, r):
        th = self.theta(r)
        if self.df_mode == "chain":
            return df_dtheta(th) * self.theta.d(r)
        return df_legacy(th)

    def temperature(self, r=None):
        if self.T_m is None:
            raise ValueError("State.T_m is required for this quantity")
        return self.T_m - self.DT(self.r_eval if r is None else r)


def constant_state(DSv, DT, gamma, theta=np.pi, r_eval=1.0, T_m=None, label=""):
    """State with r-independent fields (the CNT assumptions)."""
    return State(Quantity.constant(DSv, "DSv"), Quantity.constant(DT, "DT"),
                 Quantity.constant(gamma, "gamma"), Quantity.constant(theta, "theta"),
                 r_eval=r_eval, T_m=T_m, label=label)
