# -*- coding: utf-8 -*-
"""
γ₀, ε₁, ε₂ from the three CFM stiffnesses (cubic expansion, Hoyt, Asta, Karma 2001), with a Monte Carlo over
their statistical errors — as in md/e2_melting/E2_results.md.

    python combine_gamma0.py g100 e100 g110 e110 g111 e111
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "e2_melting"))
from analyze_cfm import combine, GEOMETRY, _AB   # noqa: E402

v = list(map(float, sys.argv[1:7]))
g = {"100": v[0], "110": v[2], "111": v[4]}
e = {"100": v[1], "110": v[3], "111": v[5]}
g0, e1, e2 = combine(g)
rng = np.random.default_rng(1)
mc = np.array([combine({k: rng.normal(g[k], e[k]) for k in g}) for _ in range(20000)])
s = mc.std(0)
print(f"γ₀ = {g0:.4f} ± {s[0]:.4f} J/m², ε₁ = {e1:.4f} ± {s[1]:.4f}, ε₂ = {e2:.4f} ± {s[2]:.4f}")
for o, (n, t) in GEOMETRY.items():
    print(f"  γ({o}) = {g0*(1+e1*_AB(n)[0]+e2*_AB(n)[1]):.4f} J/m²")
