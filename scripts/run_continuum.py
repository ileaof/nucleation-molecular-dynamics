# -*- coding: utf-8 -*-
"""
Run the author's continuum scripts in a temporary copy (the originals and their folder
are never written) and store their outputs inside the project:

    data/continuum/<alloy>/<table>.csv           nucleation table written by the script
    data/experimental/sdas_<alloy>.csv           SDAS points embedded in the SDAS script

    python scripts/run_continuum.py Al08Si06Mg02Fe
"""
import os
import runpy
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "Continuum_Mechanics_Solution" / "Scripts"

RUNS = {
    "Al08Si06Mg02Fe": dict(
        nucleation="Nucleation_model_FCC_A1_new_paper_2024_dev_full_GLH_AlSiMgFe_signal_beta_2nd_Sim_Luane.py",
        table="tabela_paper.csv",
        sdas="scriptsdasmodelAl08Si06Mg02Fe_zoom_horizontal_new_script_Ferreira_2024.py",
        source="Marques et al. (2025)",
    ),
    "Al3Cu5Nb01Fe": dict(
        nucleation="Nucleation_model_FCC_A1_new_paper_B_MRS_Meeting_2026_dev_full_GLH_AlCuNbFe_signal_beta_2nd_Sim_2026_paper.py",
        table="tablefull_paper_B_MRS_2026.csv",
        sdas=None,                       # does not execute (L342); points are parsed from the text
        sdas_text="scriptsdasmodelAl3Cu5Nb01Fe_zoom_horizontal_B_MRS_Meeting_2026.py",
        source="Mendes et al. (2023)",
    ),
}


def _copy_scripts(tmp):
    for f in SCRIPTS.glob("*.py"):
        shutil.copy2(f, tmp / f.name)


def run_nucleation(alloy):
    cfg = RUNS[alloy]
    out = ROOT / "data" / "continuum" / alloy
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        _copy_scripts(tmp)
        env = dict(os.environ, MPLBACKEND="Agg", PYTHONIOENCODING="utf-8")
        res = subprocess.run([sys.executable, cfg["nucleation"]], cwd=tmp, env=env,
                             capture_output=True, text=True, encoding="utf-8", errors="replace")
        if res.returncode != 0:
            raise RuntimeError(res.stderr[-2000:])
        shutil.copy2(tmp / cfg["table"], out / cfg["table"])
        (out / "nucleation.log").write_text(res.stdout, encoding="utf-8")
    return out / cfg["table"]


def sdas_namespace(alloy):
    """Execute the SDAS script (temporary copy, Agg backend) and return its globals."""
    cfg = RUNS[alloy]
    if cfg["sdas"] is None:
        raise ValueError(f"no runnable SDAS script registered for {alloy}")
    os.environ["MPLBACKEND"] = "Agg"
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        _copy_scripts(tmp)
        sys.path.insert(0, str(tmp))
        cwd = os.getcwd()
        try:
            os.chdir(tmp)
            ns = runpy.run_path(str(tmp / cfg["sdas"]), run_name="sdas_script")
        finally:
            os.chdir(cwd)
            sys.path.remove(str(tmp))
            import matplotlib.pyplot as plt
            plt.close("all")
    return ns


def _parse_points(path):
    """xsc/ysc/yscp/yscm assignments that are active (before the first triple-quoted block)."""
    import re
    text = path.read_text(encoding="utf-8").split('"""')[0]
    arr = {k: np.zeros(20) for k in ("xsc", "ysc", "yscp", "yscm")}
    for name, k, val in re.findall(r"^(xsc|ysc|yscp|yscm)\[(\d+)\]\s*=\s*([-+0-9.eE]+)", text, re.M):
        arr[name][int(k)] = float(val)
    return arr


def extract_sdas_points(alloy):
    cfg = RUNS[alloy]
    ns = sdas_namespace(alloy) if cfg["sdas"] else _parse_points(SCRIPTS / cfg["sdas_text"])
    n = int(np.count_nonzero(ns["xsc"]))
    out = ROOT / "data" / "experimental"
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"sdas_{alloy}.csv"
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"# SDAS experimental points from {cfg['sdas'] or cfg['sdas_text']} ({cfg['source']})\n")
        fh.write("tSL_s,SDAS_um,err_plus_um,err_minus_um\n")
        for k in range(n):
            fh.write(f"{ns['xsc'][k]},{ns['ysc'][k]},{ns['yscp'][k]},{ns['yscm'][k]}\n")
    return path


if __name__ == "__main__":
    for alloy in sys.argv[1:] or ["Al08Si06Mg02Fe"]:
        if alloy != "Al3Cu5Nb01Fe":      # its tablefull is read from the author's folder
            print(run_nucleation(alloy))
        print(extract_sdas_points(alloy))
