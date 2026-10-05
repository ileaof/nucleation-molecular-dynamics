# -*- coding: utf-8 -*-
"""
Capillary fluctuation method (Hoyt, Asta & Karma, PRL 86, 5530, 2001).

For a quasi-2D ribbon of length L (x) and thickness W (y), the interface height h(x)
has Fourier amplitudes A(k) with
        <|A(k)|²> = k_B T / (W L γ̃ k²),      γ̃ = γ + γ''  (interface stiffness)
The dump from in.cfm holds the CNA label (1 = fcc); solid columns along z give h(x)
for the two interfaces of the periodic solid slab.

Stiffnesses of (100)[100], (110)[001], (111)[1-10] are combined with the cubic expansion
(Hoyt et al. 2001)
        γ(n̂) = γ₀ [1 + ε₁(Σn_i⁴ - 3/5) + ε₂(3Σn_i⁴ + 66 n_x²n_y²n_z² - 17/7)]
to obtain γ₀ (the orientation average).

    python analyze_cfm.py cfm_<pot>_<orient>.dump T [nbins]
"""
import sys

import numpy as np

KB = 1.380649e-23


def frames(path):
    with open(path) as fh:
        while True:
            line = fh.readline()
            if not line:
                return
            if line.startswith("ITEM: TIMESTEP"):
                step = int(fh.readline())
                fh.readline(); n = int(fh.readline())
                fh.readline()
                box = [list(map(float, fh.readline().split()[:2])) for _ in range(3)]
                fh.readline()
                data = np.loadtxt([fh.readline() for _ in range(n)])
                yield step, np.array(box), data


def heights(box, data, nbins):
    """h(x) of the two interfaces of the solid slab, from the fcc fraction profile per x-column."""
    (xlo, xhi), _, (zlo, zhi) = box
    L, Lz = xhi - xlo, zhi - zlo
    x = (data[:, 1] - xlo) % L
    z = (data[:, 3] - zlo) % Lz
    sol = data[:, 4] == 1
    ib = np.minimum((x / L * nbins).astype(int), nbins - 1)
    zc = np.median(z[sol])                       # slab centre (shift to avoid wrapping)
    zs = (z - zc + Lz / 2) % Lz
    h1, h2 = np.empty(nbins), np.empty(nbins)
    for b in range(nbins):
        m = ib == b
        zz, ss = zs[m], sol[m]
        bins = np.linspace(0, Lz, 61)
        tot, _ = np.histogram(zz, bins)
        fs, _ = np.histogram(zz[ss], bins)
        if tot.sum() == 0 or fs.sum() == 0:
            return L, None, None                 # empty column or no solid: frame unusable
        phi = np.where(tot > 0, fs / np.maximum(tot, 1), 0.0)
        mid = 0.5 * (bins[1:] + bins[:-1])
        # Gibbs-like positions: solid thickness = integral of phi
        thick = np.trapezoid(phi, mid)
        centre = np.sum(phi * mid) / max(phi.sum(), 1e-12)
        h1[b], h2[b] = centre - thick / 2, centre + thick / 2
    return L, h1, h2


def main(path, T, nbins=40):
    amps, L, W, skipped = [], None, None, 0
    for k, (step, box, data) in enumerate(frames(path)):
        if k < 5:                                # skip first frames (relaxation)
            continue
        if not np.any(data[:, 4] == 1):
            skipped += 1                         # no solid at all (slab melted)
            continue
        L, h1, h2 = heights(box, data, nbins)
        W = box[1, 1] - box[1, 0]
        if h1 is None:
            skipped += 1
            continue
        for h in (h1, h2):
            A = np.fft.rfft(h - h.mean()) / nbins
            amps.append(np.abs(A[1:]) ** 2)
    if skipped:
        print(f"  WARNING: {skipped} frames skipped (no solid or empty columns)")
    if not amps:
        raise RuntimeError(f"{path}: no usable frames")
    amps = np.array(amps)
    kvec = 2 * np.pi * np.arange(1, amps.shape[1] + 1) / L          # 1/Å
    mean = amps.mean(axis=0)                                          # Å²
    n_fit = max(4, len(kvec) // 4)                                    # long-wavelength modes only
    y = 1.0 / (mean[:n_fit] * kvec[:n_fit] ** 2)                      # = W L γ̃ / kT  [1/Å²·Å⁻²...]
    stiff = KB * T * y.mean() / (W * L) * 1e20                        # J/m² (Å² -> m²)
    print(f"{path}: {len(amps)} interface samples, L = {L:.1f} Å, W = {W:.1f} Å")
    print(f"  stiffness γ̃ = {stiff:.4f} J/m²  (mean of {n_fit} lowest modes; spread {KB*T*y.std()/(W*L)*1e20:.4f})")
    return stiff




# ------------------------------------------------------------ anisotropy fit
# in.cfm geometry: (normal n, fluctuation direction t)
GEOMETRY = {"100": ((0, 0, 1), (1, 0, 0)),
            "110": ((1, 1, 0), (0, 0, 1)),
            "111": ((1, 1, 1), (1, -1, 0))}


def _AB(n):
    n = np.asarray(n, float) / np.linalg.norm(n)
    s4 = np.sum(n ** 4)
    return s4 - 3 / 5, 3 * s4 + 66 * np.prod(n ** 2) - 17 / 7


def stiffness_coefficients(n, t, h=1e-3):
    """γ̃/γ₀ = 1 + ε₁·a + ε₂·b, with a = A + A'', b = B + B'' (rotation of n towards t)."""
    n = np.asarray(n, float) / np.linalg.norm(n)
    t = np.asarray(t, float) / np.linalg.norm(t)
    rot = lambda th: np.cos(th) * n + np.sin(th) * t
    vals = [np.array(_AB(rot(th))) for th in (-h, 0.0, h)]
    second = (vals[0] - 2 * vals[1] + vals[2]) / h ** 2
    a, b = vals[1] + second
    return a, b


def combine(stiff):
    """stiff = {"100": γ̃, "110": γ̃, "111": γ̃}  ->  γ₀, ε₁, ε₂ (cubic expansion of Hoyt et al. 2001)."""
    M, y = [], []
    for key, g in stiff.items():
        a, b = stiffness_coefficients(*GEOMETRY[key])
        M.append([1.0, a, b]); y.append(g)
    g0, g0e1, g0e2 = np.linalg.solve(np.array(M), np.array(y))
    return g0, g0e1 / g0, g0e2 / g0


def main_all(pot, T, nbins=40):
    stiff = {o: main(f"cfm_{pot}_{o}.dump", T, nbins) for o in GEOMETRY}
    g0, e1, e2 = combine(stiff)
    print(f"\nγ₀ = {g0:.4f} J/m²,  ε₁ = {e1:.4f},  ε₂ = {e2:.4f}")
    for o, (n, t) in GEOMETRY.items():
        print(f"  γ({o}) = {g0*(1+e1*_AB(n)[0]+e2*_AB(n)[1]):.4f} J/m²")
    return g0, e1, e2


if __name__ == "__main__":
    # python analyze_cfm.py <pot> <T>        (all three orientations)
    # python analyze_cfm.py <dump> <T>       (one file)
    if sys.argv[1].endswith(".dump"):
        main(sys.argv[1], float(sys.argv[2]))
    else:
        main_all(sys.argv[1], float(sys.argv[2]))
