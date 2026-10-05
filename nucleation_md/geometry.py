# -*- coding: utf-8 -*-
"""
Spherical-cap geometry.

    f(θ)     = 2 - 3cosθ + cos³θ          (f(π) = 4: homogeneous sphere)
    df/dθ    = 3 sin³θ
    V_cap    = (1/3)·π·r³·f(θ)
    A_SL     = 2·π·r²·(1 - cosθ)          (solid–liquid area; 4πr² at θ = π)
"""
import numpy as np


def f_theta(theta):
    c = np.cos(theta)
    return 2.0 - 3.0 * c + c ** 3


def df_dtheta(theta):
    return 3.0 * np.sin(theta) ** 3


def df_legacy(theta):
    """
    `dfthetadr(theta)` of the continuum script (L391), reproduced verbatim:
        -3·(2 - 3cosθ + cos³θ) - (1 - cosθ)·(2 - cosθ - cos²θ)
    Function of θ only, dimensionless (-16 at θ = π). Kept for regression against
    the continuum tables; it is not df/dθ·dθ/dr.
    """
    c = np.cos(theta)
    return -3.0 * (2.0 - 3.0 * c + c ** 3) - (1.0 - c) * (2.0 - c - c ** 2)


def cap_volume(r, theta):
    return np.pi * r ** 3 * f_theta(theta) / 3.0


def cap_area(r, theta):
    return 2.0 * np.pi * r ** 2 * (1.0 - np.cos(theta))
