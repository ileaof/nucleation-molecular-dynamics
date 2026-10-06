# -*- coding: utf-8 -*-
"""
E4 NEMD analysis of one run directory (log.lammps, Tz.txt, T3d.txt, solid.dump), per 10 ps window:

  seed      largest connected cluster of solid atoms (neighbours < 3.6 Å); N_cr = N / f_s(T0);
            r_eq = (3 N_cr a³/4 / 4π)^(1/3); centre z_c; extents to the cold (−z) and hot (+z) sides along the axis;
            gyration-tensor eigenvalues → asphericity b = λ_z − (λ_x + λ_y)/2 normalised by Rg²
  gradient  LOCAL ∇T at the interface: weighted linear fit T = T_s + g·(x − x_c) over the 3D bins of a liquid shell
            r_eq + 2 Å < d < r_eq + 14 Å around the seed centre (author: the ∇T of Γ = A·∇T is the one at the interface)
  ΔT        T_m − T_s (shell mean at the interface), also cold-side and hot-side halves
  Γ         ΔT·r_eq/2  versus  A·|∇T| with A = 4π r_eq² (θ = π, homogeneous)
plus the far-field gradient from the 1D profile (between the slabs, outside the seed).

    python analyze_nemd.py <run_dir> [T_m] [T0]
"""
import re
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

CAL_T = np.array([760.0, 820.0, 880.0, 920.0])        # ../e4_seeding/in.calib (q6·q6 > 0.5, >= 7 bonds)
CAL_F = np.array([0.9285, 0.8948, 0.8728, 0.8472])


def read_ave_chunk(path, ncols):
    """{timestep: array(nchunks, ncols)} from a fix ave/chunk file."""
    out, cur, rows, need = {}, None, [], 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) == 3 and need == 0:
                if cur is not None:
                    out[cur] = np.array(rows, float)
                cur, need, rows = int(p[0]), int(p[1]), []
                continue
            rows.append(p[:ncols]); need -= 1
    if cur is not None and rows:
        out[cur] = np.array(rows, float)
    return out


def read_dump(path):
    with open(path) as fh:
        while True:
            line = fh.readline()
            if not line:
                return
            if line.startswith("ITEM: TIMESTEP"):
                step = int(fh.readline())
                fh.readline(); n = int(fh.readline())
                fh.readline()
                box = np.array([list(map(float, fh.readline().split()[:2])) for _ in range(3)])
                fh.readline()
                d = np.array([fh.readline().split() for _ in range(n)], float) if n else np.zeros((0, 5))
                yield step, box, d


def largest_cluster(pos, L, rc=3.6):
    if len(pos) == 0:
        return pos
    t = cKDTree(pos % L, boxsize=L)
    pairs = t.query_pairs(rc, output_type="ndarray")
    from scipy.sparse import coo_matrix
    m = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(len(pos), len(pos)))
    _, lab = connected_components(m, directed=False)
    big = np.bincount(lab).argmax()
    return pos[lab == big]


def unwrap_cluster(p, L):
    ref = p[0]
    d = p - ref
    d -= L * np.round(d / L)
    q = ref + d
    c = q.mean(0)
    d = q - c
    d -= L * np.round(d / L)
    return c + d


def main(run, TM=936.0, T0=None):
    run = Path(run)
    log = (run / "log.lammps").read_text(errors="ignore")
    a = float(re.search(r"lattice at T0: a = ([\d.]+)", log).group(1))
    G = float(re.search(r"G = ([\d.eE+-]+) K/A", log).group(1))
    m = re.search(r"slabs: cold ([\d.]+) K at z = [\d.]+, hot ([\d.]+) K", log)
    T0 = T0 or (float(m.group(1)) + float(m.group(2))) / 2           # slabs at T0 ∓ G·L/4
    fs = float(np.interp(T0, CAL_T, CAL_F))
    va = a**3 / 4

    T3 = read_ave_chunk(run / "T3d.txt", 7)       # chunk x y z count temp dens
    Tz = read_ave_chunk(run / "Tz.txt", 5)        # chunk z count temp dens
    frames = list(read_dump(run / "solid.dump"))
    out = []
    for step in sorted(T3):
        fr = [f for f in frames if step - 5000 < f[0] <= step]
        if not fr:
            continue
        geo = []
        for st, box, d in fr:
            L = box[:, 1] - box[:, 0]
            p = largest_cluster(d[:, 1:4], L)
            if len(p) < 20:
                continue
            q = unwrap_cluster(p, L)
            c = q.mean(0)
            gy = np.cov((q - c).T)
            lam = np.sort(np.linalg.eigvalsh(gy))
            ez = np.linalg.eigh(gy)
            Rg2 = np.trace(gy)
            dz = q[:, 2] - c[2]
            geo.append(dict(N=len(p) / fs, c=c, zlo=np.percentile(-dz, 98), zhi=np.percentile(dz, 98),
                            rxy=np.percentile(np.hypot(q[:, 0] - c[0], q[:, 1] - c[1]), 98),
                            b=(gy[2, 2] - (gy[0, 0] + gy[1, 1]) / 2) / Rg2, L=L))
        if not geo:
            out.append(dict(step=step, N=0)); continue
        N = np.mean([g["N"] for g in geo]); c = np.mean([g["c"] for g in geo], 0); L = geo[0]["L"]
        r = (3 * N * va / (4 * np.pi)) ** (1 / 3)
        b3 = T3[step]
        x = b3[:, 1:4]; cnt = b3[:, 4]; T = b3[:, 5]
        d = x - c; d -= L * np.round(d / L)
        dist = np.linalg.norm(d, axis=1)
        sh = (dist > r + 2) & (dist < r + 14) & (cnt > 0)
        A = np.column_stack([np.ones(sh.sum()), d[sh]])
        w = np.sqrt(cnt[sh])
        coef, *_ = np.linalg.lstsq(A * w[:, None], T[sh] * w, rcond=None)
        Ts, g = coef[0], coef[1:]                                    # K, K/Å
        Tcold = np.average(T[sh & (d[:, 2] < 0)], weights=cnt[sh & (d[:, 2] < 0)])
        Thot = np.average(T[sh & (d[:, 2] > 0)], weights=cnt[sh & (d[:, 2] > 0)])
        # far field: 1D profile between the slabs, |z − z_c| in (r + 20, L/4 − 15)
        z1 = Tz[step]; zz = z1[:, 1] - c[2]; ok = (np.abs(zz) > r + 20) & (np.abs(zz) < L[2] / 4 - 15) & (z1[:, 2] > 0)
        gfar = np.polyfit(zz[ok], z1[ok, 3], 1)[0] if ok.sum() > 3 else np.nan
        dT = TM - Ts
        Ar = 4 * np.pi * (r * 1e-10) ** 2
        out.append(dict(step=step, N=N, r=r / 10, zc=c[2], zlo=np.mean([g_["zlo"] for g_ in geo]),
                        zhi=np.mean([g_["zhi"] for g_ in geo]), rxy=np.mean([g_["rxy"] for g_ in geo]),
                        b=np.mean([g_["b"] for g_ in geo]), Ts=Ts, Tcold=Tcold, Thot=Thot, gloc=g * 1e10,
                        gfar=gfar * 1e10, dT=dT, Gam=dT * r * 1e-10 / 2, AgradT=Ar * np.linalg.norm(g) * 1e10))
    print(f"# {run.name}: G imposto = {G} K/Å = {G*1e10:.2e} K/m, T0 = {T0}, T_m = {TM}, f_s = {fs:.3f}\n")
    print("| t [ps] | N_cr | r_eq [nm] | z_c [Å] | ext. fria/quente/lateral [Å] | b | T_int [K] | T fria/quente [K] "
          "| ∇T local (x,y,z) [10⁹ K/m] | ∇T far [10⁹ K/m] | ΔT [K] | ΔT·r/2 [m·K] | A·|∇T| [m·K] |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    t0 = out[0]["step"] - 5000 if out else 0
    for o in out:
        if not o.get("N"):
            print(f"| {(o['step']-t0)*0.002:.0f} | 0 | — |"); continue
        gl = o["gloc"] / 1e9
        print(f"| {(o['step']-t0)*0.002:.0f} | {o['N']:.0f} | {o['r']:.3f} | {o['zc']:.1f} | {o['zlo']:.1f}/{o['zhi']:.1f}/{o['rxy']:.1f} "
              f"| {o['b']:+.3f} | {o['Ts']:.1f} | {o['Tcold']:.1f}/{o['Thot']:.1f} | {gl[0]:+.2f}, {gl[1]:+.2f}, {gl[2]:+.2f} "
              f"| {o['gfar']/1e9:+.2f} | {o['dT']:.1f} | {o['Gam']:.3e} | {o['AgradT']:.3e} |")

    ok = [o for o in out if o.get("N")]
    if len(ok) >= 3:
        t = np.array([(o["step"] - t0) * 0.002 for o in ok])
        g = np.array([o["gloc"] for o in ok]) / 1e9
        gm, ge = g.mean(0), g.std(0, ddof=1) / np.sqrt(len(g))
        Ar = np.array([4 * np.pi * (o["r"] * 1e-9) ** 2 for o in ok])
        Agz = Ar * g[:, 2] * 1e9
        summ = dict(run=run.name, G=G * 1e10, T0=T0, dNdt=np.polyfit(t, [o["N"] for o in ok], 1)[0],
                    dzdt=np.polyfit(t, [o["zc"] for o in ok], 1)[0], gm=gm, ge=ge,
                    gam=np.mean([o["Gam"] for o in ok]), Agz=Agz.mean(), Agz_e=Agz.std(ddof=1) / np.sqrt(len(Agz)),
                    r=np.mean([o["r"] for o in ok]), dT=np.mean([o["dT"] for o in ok]),
                    dTside=np.mean([o["Thot"] - o["Tcold"] for o in ok]))
        print(f"\nRESUMO | {run.name} | G = {G*1e10:.2e} K/m | dN/dt = {summ['dNdt']:+.1f} át/ps | dz_c/dt = {summ['dzdt']:+.3f} Å/ps"
              f" | <∇T local> = ({gm[0]:+.2f}±{ge[0]:.2f}, {gm[1]:+.2f}±{ge[1]:.2f}, {gm[2]:+.2f}±{ge[2]:.2f})·10⁹ K/m"
              f" | <ΔT·r/2> = {summ['gam']:.3e} | <A·(∇T)_z> = {summ['Agz']:.2e} ± {summ['Agz_e']:.1e} m·K"
              f" | T_quente − T_fria = {summ['dTside']:+.1f} K")
        return summ


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], float(a[1]) if len(a) > 1 else 936.0, float(a[2]) if len(a) > 2 else None)
