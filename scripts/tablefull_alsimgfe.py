# -*- coding: utf-8 -*-
"""
Full nucleation table (manuscript Tables 1-2 format) for Al-0.8Si-0.6Mg-0.2Fe.

The author's AlSiMgFe nucleation script writes only tabela_paper.csv (10 columns). The
AlCuNbFe script (same loop) also writes tablefull. This tool makes a TEMPORARY COPY of the
AlSiMgFe script, inserts the tablefull blocks of the AlCuNbFe script verbatim (references to
its lines in the comments), runs it, and stores the outputs in data/continuum/Al08Si06Mg02Fe/.
The original script is never modified. Physics is unchanged: only extra arrays and output.

    python scripts/tablefull_alsimgfe.py
"""
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
# in-project copy of C:\MacbookPro\ufpa\PPGEM\SDAS Al_Si_Mg_Fe\Sim Luane Final (copied 2026-10-04, byte-identical)
SRC_DIR = ROOT / "Continuum_Mechanics_Solution" / "Sim_Luane_Final"
SCRIPT = "Nucleation_model_FCC_A1_new_paper_2024_dev_full_GLH_AlSiMgFe_signal_beta_2nd_Sim_Luane.py"
OUT = ROOT / "data" / "continuum" / "Al08Si06Mg02Fe"

HEADER = ("i;dTdr;DT;r_hom_1st;r_het_1st;theta_1st;dfthetadr_1st;r_hom_2nd;r_het_2nd;theta_2nd;"
          "dfthetadr_2nd;sigma;surface_stress;gam_hom;dfgamdr_ana;GT_het;GT_het_2nd;DSv_hom;DSv_het;"
          "DSs_hom;DSs_het;DSc;GB;GS;GC;dDSv_homdr")

PATCHES = [
    # (anchor in the AlSiMgFe script, text inserted AFTER the anchor line)
    ("rc_v = zeros(n)\n",
     "# --- inserted (AlCuNbFe script, array declarations)\n"
     "dfthetadr_2nd_v = zeros(n); gb_v = zeros(n); gs_v = zeros(n); gc_v = zeros(n)\n"
     "DSv_hom_v = zeros(n); dDSv_homdr_v = zeros(n); DSv_het_v = zeros(n)\n"
     "DSs_hom_v = zeros(n); DSs_het_v = zeros(n); DSc_v = zeros(n)\n"),
    ("file = open('tabela_paper.csv','w')\n",
     "# --- inserted (AlCuNbFe script L414-L432 and L853/L865)\n"
     "def GB(Delta_Sv_hom, DT_):\n    return Delta_Sv_hom * DT_\n"
     "def GS(gammasl):\n    return gammasl\n"
     "def GC(Delta_Sv_hom, DT_, rc_het_2nd, gammasl, theta_het_2nd):\n"
     "    return gammasl / ftheta(theta_het_2nd) * dfthetadr(theta_het_2nd)\n"
     "filefull = open('tablefull_AlSiMgFe.csv','w')\n"
     f"filefull.write('{HEADER}\\n')\n"),
    ("       DeltaSV_Total_hom_v[i] = dDeltaSvhomDTdr_v[i] * dreq / DT_v[i]\n",
     "       # --- inserted (AlCuNbFe script L900, L937-L941)\n"
     "       DSv_hom_v[0] = Delta_Sv_hom\n"
     "       DSv_hom_v[i] = DSv_hom_v[i-1] + dDeltaSvhomDTdr_v[i] * dreq / DT_v[i]\n"
     "       dDSv_homdr_v[i] = (DSv_hom_v[i] - DSv_hom_v[i-1]) / (rc_v[i]-rc_v[i-1])\n"
     "       dfthetadr_2nd_v[i] = dfthetadr(theta_2nd_v[i])\n"),
    ("    Delta_Sv_het = Delta_Sv_hom * ftheta(theta_v[i])\n",
     "    DSv_het_v[i] = DSv_hom_v[i] * ftheta(theta_2nd_v[i])/4.0   # inserted (AlCuNbFe L962)\n"),
    ("    DeltaSV_Total_het_v[i] = dDeltaSvhomDTfthetadr_v[i] * dreq/ ( ftheta(theta_v[i])*DT_v[i])\n",
     "    # --- inserted (AlCuNbFe script L973-L979)\n"
     "    DSs_hom_v[i] = dfgamdr_ana_v[i] / DT_v[i]\n"
     "    DSs_het_v[i] = DSs_hom_v[i] * ftheta(theta_2nd_v[i])/4.0\n"
     "    DSc_v[i] = gam_hom_v[i] / DT_v[i] * 1.0/ftheta(theta_2nd_v[i]) * dfthetadr(theta_2nd_v[i])\n"
     "    gb_v[i] = GB(DSv_hom_v[i], DT_v[i])\n"
     "    gs_v[i] = GS( gam_hom_v[i])\n"
     "    gc_v[i] = GC(DSv_hom_v[i], DT_v[i],r_het_2nd_v[i], gam_hom_v[i],  theta_2nd_v[i])\n"),
    ("{GT_het_v[i]:g};{GT_het_2nd_v[i]:g}\\n')\n",
     "     filefull.write(f'{i};{dTdr_v[i]:g};{DT_v[i]:g};{r_hom_v[i]:.5e};{r_het_v[i]:.5e};{theta_v[i]:g};"
     "{dfthetadr(theta_v[i]):g};{r_hom_2nd_v[i]:.5e};{r_het_2nd_v[i]:.5e};{theta_2nd_v[i]:g};"
     "{dfthetadr(theta_2nd_v[i]):g};{sigma_v[i]:g};{surface_stress_v[i]:g};{gam_hom_v[i]:g};"
     "{dfgamdr_ana_v[i]:g};{GT_het_v[i]:.5e};{GT_het_2nd_v[i]:.5e};{DSv_hom_v[i]:g};{DSv_het_v[i]:g};"
     "{DSs_hom_v[i]:g};{DSs_het_v[i]:g};{DSc_v[i]:g};{gb_v[i]:g};{gs_v[i]:g};{gc_v[i]:g};"
     "{dDSv_homdr_v[i]:.5e}\\n')   # inserted (AlCuNbFe L1022, without Vm/Density/mu_Al)\n"),
    ("file.close()\n#quit()\n", "filefull.close()   # inserted\n"),
]


def patched_source():
    s = (SRC_DIR / SCRIPT).read_text(encoding="utf-8")
    for anchor, add in PATCHES:
        n = s.count(anchor)
        if n != 1:
            raise RuntimeError(f"anchor found {n} times: {anchor!r}")
        s = s.replace(anchor, anchor + add)
    return s


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        for f in SRC_DIR.glob("*.py"):
            shutil.copy2(f, tmp / f.name)
        (tmp / "patched.py").write_text(patched_source(), encoding="utf-8")
        env = dict(os.environ, MPLBACKEND="Agg", PYTHONIOENCODING="utf-8")
        r = subprocess.run([sys.executable, "patched.py"], cwd=tmp, env=env, capture_output=True,
                           text=True, encoding="utf-8", errors="replace")
        if r.returncode:
            raise RuntimeError(r.stderr[-3000:])
        shutil.copy2(tmp / "tablefull_AlSiMgFe.csv", OUT / "tablefull_AlSiMgFe.csv")
        shutil.copy2(tmp / "tabela_paper.csv", OUT / "tabela_paper_from_patched.csv")
        (OUT / "patched_script_for_reference.py").write_text(patched_source(), encoding="utf-8")

    # check: columns shared with the unpatched output must be identical
    full = pd.read_csv(OUT / "tablefull_AlSiMgFe.csv", sep=";")
    ref = pd.read_csv(OUT / "tabela_paper.csv", sep=";", skiprows=1, header=None,
                      names=["i", "dTdr", "DT", "r_hom_2nd", "r_het_2nd", "theta_2nd", "sigma",
                             "surface_stress", "gam_hom", "GT_het", "GT_het_2nd"])
    rows = full["i"].isin(ref["i"])
    for c in ["dTdr", "DT", "r_hom_2nd", "r_het_2nd", "theta_2nd", "sigma", "surface_stress",
              "gam_hom", "GT_het", "GT_het_2nd"]:
        a = full.loc[rows, c].to_numpy(float)
        b = ref.set_index("i").loc[full.loc[rows, "i"], c].to_numpy(float)
        print(f"{c:15s} max rel diff vs original output: {np.max(np.abs(a - b) / np.abs(b)):.2e}")
    print(f"{len(full)} rows ->", OUT / "tablefull_AlSiMgFe.csv")


if __name__ == "__main__":
    main()
