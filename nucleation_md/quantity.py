# -*- coding: utf-8 -*-
"""
Scalar fields q(r) with their radial derivative dq/dr.

Every r-dependent input of the model (ΔS_V, ΔT, γ_SL, θ) is a Quantity, so the
same model code runs on analytic continuum expressions, on continuum tables or
on MD seed data. All values are SI.
"""
import numpy as np
from scipy.interpolate import CubicSpline


def _scalar(x):
    """Return a numpy scalar for 0-d input, the array otherwise."""
    return np.asarray(x, dtype=float)[()]


class Quantity:
    """
    q(r) and dq/dr.

    func      - callable r -> q(r)
    deriv     - callable r -> dq/dr; if None a central difference is used,
                dq/dr = (q(r+h) - q(r-h)) / (2h),  h = rel_step * r
    name      - label used in reports
    """

    def __init__(self, func, deriv=None, name="", rel_step=1e-5):
        self._f = func
        self._df = deriv
        self.name = name
        self.rel_step = rel_step

    def __call__(self, r):
        return _scalar(self._f(r))

    def d(self, r):
        if self._df is not None:
            return _scalar(self._df(r))
        r = np.asarray(r, dtype=float)
        h = self.rel_step * np.abs(r)
        return _scalar((self._f(r + h) - self._f(r - h)) / (2.0 * h))

    @property
    def has_analytic_derivative(self):
        return self._df is not None

    @classmethod
    def constant(cls, value, name=""):
        v = float(value)
        return cls(lambda r: np.full_like(np.asarray(r, dtype=float), v),
                   lambda r: np.zeros_like(np.asarray(r, dtype=float)),
                   name)

    @classmethod
    def tabulated(cls, r, y, name=""):
        """
        Cubic spline through (r_k, y_k), e.g. MD seeds of several radii or rows of the
        continuum table. dq/dr is the spline derivative. Outside the data range the
        spline extrapolates; callers should stay inside [min r, max r].
        """
        r = np.asarray(r, dtype=float)
        y = np.asarray(y, dtype=float)
        order = np.argsort(r)
        r, y = r[order], y[order]
        if np.any(np.diff(r) <= 0.0):
            raise ValueError(f"{name}: radii must be distinct")
        spline = CubicSpline(r, y)
        dspline = spline.derivative()
        q = cls(spline, dspline, name)
        q.r_range = (r[0], r[-1])
        return q
