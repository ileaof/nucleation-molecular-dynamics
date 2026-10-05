"""Every tablefull column used by nucleation_md is recomputed from the script's equations
(scripts/verify_continuum.py) and must match to the printing precision (6 digits)."""
import runpy
from pathlib import Path

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "verify_continuum.py"
CHECKS = runpy.run_path(str(SCRIPT), run_name="verify")["checks"]


@pytest.mark.parametrize("name", list(CHECKS))
def test_column_recomputed_from_continuum_equations(name):
    calc, tab = CHECKS[name]
    tab = np.asarray(tab, dtype=float)
    assert np.max(np.abs(calc - tab) / np.abs(tab)) < 1e-5
