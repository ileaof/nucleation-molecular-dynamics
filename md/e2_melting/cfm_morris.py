# -*- coding: utf-8 -*-
"""
CFM interface heights from the Morris (2002) orientational order parameter, replacing the
CNA label (too noisy at T_m: white-noise floor in <|A(k)|²>).

    φ_i  = (1/12) Σ_j min_k |r_ij - v_k|²     over the 12 nearest neighbours j,
           v_k = ideal fcc bond vectors (a/2)<110> rotated to the ribbon frame
    φ̄_i = mean of φ over atoms within r_cg (spatial coarse-graining)
    h(x) = z where the column profile of φ̄ crosses (φ_s + φ_l)/2, for both interfaces

    python cfm_morris.py <orient> [T] [nbins] [pot]     (reads cfm_<pot>_<orient>.dump in the cwd)
"""
import sys

import numpy as np
from scipy.spatial import cKDTree

from analyze_cfm import frames, combine, GEOMETRY, _AB, KB

A0 = {"borovikov": 4.016, "nep89": 3.9899}   # 0 K fcc lattice of each potential
ORIENT = {"100": [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
          "110": [(0, 0, 1), (1, -1, 0), (1, 1, 0)],
          "111": [(1, -1, 0), (1, 1, -2), (1, 1, 1)]}


def ideal_bonds(orient, a, T_expand=1.0):
    R = np.array([np.array(v, float) / np.linalg.norm(v) for v in ORIENT[orient]])   # rows: box axes in crystal frame
    nn = []
    for i in (-1, 1):
        for j in (-1, 1):
            nn += [(i, j, 0), (i, 0, j), (0, i, j)]
    nn = np.array(nn, float) * a / 2 * T_expand
    return nn @ R.T                                                                   # crystal -> box frame


def order_parameter(pos, box, bonds, r_cg=4.0):
    L = box[:, 1] - box[:, 0]
    p = (pos - box[:, 0]) % L
    tree = cKDTree(p, boxsize=L)
    _, idx = tree.query(p, k=13)
    rij = p[idx[:, 1:]] - p[:, None, :]
    rij -= L * np.round(rij / L)
    d2 = ((rij[:, :, None, :] - bonds[None, None, :, :]) ** 2).sum(-1).min(-1)       # (N, 12)
    phi = d2.mean(1)
    pairs = tree.query_ball_point(p, r_cg)
    phibar = np.array([phi[nb].mean() for nb in pairs])
    return p, phibar


def heights(p, phibar, L, nbins, dz=1.0):
    nz = int(L[2] / dz)
    ix = np.minimum((p[:, 0] / L[0] * nbins).astype(int), nbins - 1)
    iz = np.minimum((p[:, 2] / L[2] * nz).astype(int), nz - 1)
    grid_s = np.zeros((nbins, nz)); grid_n = np.zeros((nbins, nz))
    np.add.at(grid_s, (ix, iz), phibar); np.add.at(grid_n, (ix, iz), 1)
    prof = np.where(grid_n > 0, grid_s / np.maximum(grid_n, 1), np.nan)
    # fill empty z-bins by interpolation along z
    for b in range(nbins):
        good = ~np.isnan(prof[b])
        prof[b] = np.interp(np.arange(nz), np.where(good)[0], prof[b, good], period=nz)
    lo, hi = np.percentile(prof, 10), np.percentile(prof, 90)        # solid ~ low φ, liquid ~ high φ
    mid = 0.5 * (lo + hi)
    solid = (prof < mid).mean(0) > 0.5                                 # z-bins solid in most columns
    ang = 2 * np.pi * np.where(solid)[0] / nz                          # circular mean (slab may wrap)
    zc = int(round((np.arctan2(np.sin(ang).mean(), np.cos(ang).mean()) % (2 * np.pi)) / (2 * np.pi) * nz)) % nz
    shift = nz // 2 - zc
    prof = np.roll(prof, shift, axis=1)
    h1, h2 = np.empty(nbins), np.empty(nbins)
    z = (np.arange(nz) + 0.5) * dz
    for b in range(nbins):
        f = prof[b] - mid                                              # < 0 inside the solid
        c = nz // 2
        if f[c] >= 0:
            return None, None
        i = c
        while i > 0 and f[i] < 0:
            i -= 1
        j = c
        while j < nz - 1 and f[j] < 0:
            j += 1
        if i == 0 or j == nz - 1:
            return None, None
        h1[b] = z[i] + (z[i + 1] - z[i]) * (0 - f[i]) / (f[i + 1] - f[i])
        h2[b] = z[j - 1] + (z[j] - z[j - 1]) * (0 - f[j - 1]) / (f[j] - f[j - 1])
    return h1, h2


def stiffness(orient, T=936.0, nbins=40, pot="borovikov", skip=5, nmodes=6):
    bonds = ideal_bonds(orient, A0[pot])
    amps, used, bad = [], 0, 0
    for k, (step, box, data) in enumerate(frames(f"cfm_{pot}_{orient}.dump")):
        if k < skip:
            continue
        L = box[:, 1] - box[:, 0]
        p, phibar = order_parameter(data[:, 1:4], box, bonds)
        h1, h2 = heights(p, phibar, L, nbins)
        if h1 is None:
            bad += 1
            continue
        used += 1
        for h in (h1, h2):
            A = np.fft.rfft(h - h.mean()) / nbins
            amps.append(np.abs(A[1:]) ** 2)
    amps = np.array(amps)
    kv = 2 * np.pi * np.arange(1, amps.shape[1] + 1) / L[0]
    m, se = amps.mean(0), amps.std(0) / np.sqrt(len(amps))
    W = L[1]
    g_k = KB * T / (W * L[0] * m * kv ** 2) * 1e20
    np.savez(f"spectrum_{pot}_{orient}.npz", k=kv, m=m, se=se, L=L[0], W=W, T=T, nframes=len(amps))
    print(f"({orient}) frames used {used}, rejected {bad}; L = {L[0]:.1f} Å, W = {W:.1f} Å")
    for i in range(8):
        print(f"   mode {i+1}: <|A|²> = {m[i]:.4f} ± {se[i]:.4f} Å²   γ̃_k = {g_k[i]:.4f} J/m²")
    # capillary window: modes kmin..kmax (long modes not relaxed in 400 ps; short modes reach the
    # height-resolution floor). Log-log slope of <|A|²> vs k must be close to -2 there.
    w = slice(3, 10)
    slope = np.polyfit(np.log(kv[w]), np.log(m[w]), 1)[0]
    g = np.mean(g_k[w]); g_se = np.std(g_k[w]) / np.sqrt(len(g_k[w]))
    print(f"   window modes 4-10: log-log slope = {slope:.2f} (ideal -2);  γ̃ = {g:.4f} ± {g_se:.4f} J/m²")
    return g, g_se


if __name__ == "__main__":
    o = sys.argv[1]
    T = float(sys.argv[2]) if len(sys.argv) > 2 else 936.0
    nb = int(sys.argv[3]) if len(sys.argv) > 3 else 40
    pot = sys.argv[4] if len(sys.argv) > 4 else "borovikov"
    stiffness(o, T, nb, pot)
