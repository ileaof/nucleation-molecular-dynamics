# -*- coding: utf-8 -*-
"""
Drop-in replacement for `scipy.misc.derivative`, which was deprecated and
removed in SciPy 1.12. Same signature and finite-difference weight tables
as the original (n = 1..4, order = 3..9), so call sites are unchanged.

    derivative(func, x0, dx=1.0, n=1, args=(), order=3)

For the default n=1, order=3 this evaluates the central 3-point stencil
    f'(x0) = (f(x0+dx) - f(x0-dx)) / (2*dx)
exactly as the old scipy.misc.derivative did. Uses plain Python so it has no
dependency on the removed `numpy.asfarray` and works under NumPy 2.x.
"""


def derivative(func, x0, dx=1.0, n=1, args=(), order=3):
    x0 = float(x0)
    if n == 0:
        return func(*((x0,) + args))
    if order < n + 1:
        raise ValueError("'order' must be at least 'n+1'.")
    if order % 2 == 0:
        raise ValueError("'order' must be odd.")
    if n == 1:
        if order == 3:
            weights = [-0.5, 0.0, 0.5]
        elif order == 5:
            weights = [1 / 12, -8 / 12, 0.0, 8 / 12, -1 / 12]
        elif order == 7:
            weights = [-1 / 60, 9 / 60, -45 / 60, 0.0, 45 / 60, -9 / 60, 1 / 60]
        elif order == 9:
            weights = [3 / 840, -32 / 840, 168 / 840, -672 / 840, 0.0,
                       672 / 840, -168 / 840, 32 / 840, -3 / 840]
        else:
            raise NotImplementedError("order > 9 not supported")
    elif n == 2:
        if order == 3:
            weights = [1.0, -2.0, 1.0]
        elif order == 5:
            weights = [-1 / 12, 16 / 12, -30 / 12, 16 / 12, -1 / 12]
        elif order == 7:
            weights = [2 / 180, -27 / 180, 270 / 180, -490 / 180,
                       270 / 180, -27 / 180, 2 / 180]
        elif order == 9:
            weights = [-9 / 5040, 128 / 5040, -1008 / 5040, 8064 / 5040,
                       -14350 / 5040, 8064 / 5040, -1008 / 5040, 128 / 5040,
                       -9 / 5040]
        else:
            raise NotImplementedError("order > 9 not supported")
    elif n == 3:
        if order == 5:
            weights = [-0.5, 1.0, 0.0, -1.0, 0.5]
        elif order == 7:
            weights = [1 / 8, -9 / 8, 45 / 8, 0.0, -45 / 8, 9 / 8, -1 / 8]
        elif order == 9:
            weights = [-9 / 80, 108 / 80, -1008 / 80, 0.0, 5040 / 80,
                       -1008 / 80, 108 / 80, -9 / 80]
        else:
            raise RuntimeError("order < 5 or > 9 not supported for 3rd derivative")
    elif n == 4:
        if order == 5:
            weights = [1.0, -4.0, 6.0, -4.0, 1.0]
        elif order == 7:
            weights = [-3 / 24, 32 / 24, -168 / 24, 672 / 24, -1032 / 24,
                       672 / 24, -168 / 24, 32 / 24, -3 / 24]
        elif order == 9:
            weights = [9 / 720, -128 / 720, 1344 / 720, -6720 / 720,
                       20160 / 720, -6720 / 720, 1344 / 720, -128 / 720, 9 / 720]
        else:
            raise RuntimeError("order < 5 or > 9 not supported for 4th derivative")
    else:
        raise RuntimeError("only derivatives of order 1 to 4 supported")
    val = 0.0
    ho = order // 2
    for k in range(order):
        x = x0 + (k - ho) * dx
        val += weights[k] * func(*((x,) + args))
    return val / (dx ** n)