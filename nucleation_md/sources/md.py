# -*- coding: utf-8 -*-
"""
MD source: builds a State from quantities measured on seeds of several radii.

Each field may be a scalar (r-independent, e.g. ΔT imposed by the thermostat) or a
sequence aligned with `radii` (one value per seed; ∂/∂r from the spline through them).
The LAMMPS post-processing of stages E2–E5 produces these arrays.
"""
import numpy as np

from ..quantity import Quantity
from ..state import State


def _field(radii, value, name):
    if np.ndim(value) == 0:
        return Quantity.constant(value, name)
    value = np.asarray(value, dtype=float)
    if radii is None or len(radii) != len(value):
        raise ValueError(f"{name}: needs one value per seed radius")
    if len(value) < 4:
        raise ValueError(f"{name}: at least 4 seed radii are needed for a spline derivative")
    return Quantity.tabulated(radii, value, name)


def state_from_seeds(r_eval, DSv, DT, gamma, theta=np.pi, radii=None, T_m=None,
                     df_mode="chain", label="MD"):
    """
    r_eval - radius at which the explicit formulas are evaluated [m]
    DSv    - ΔS_V [J m^-3 K^-1] (< 0)
    DT     - ΔT [K] (> 0)
    gamma  - γ_SL [J m^-2]
    theta  - θ [rad]
    radii  - seed radii [m] when any field is given per seed
    """
    return State(_field(radii, DSv, "DSv"), _field(radii, DT, "DT"),
                 _field(radii, gamma, "gamma"), _field(radii, theta, "theta"),
                 r_eval=r_eval, T_m=T_m, df_mode=df_mode, label=label)
