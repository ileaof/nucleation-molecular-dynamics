"""
Sign convention (docs/CONVENTIONS.md) checked on every continuum row, plus the
mandatory trend test of I against experimental ΔT and ∇T.
"""
from pathlib import Path

import numpy as np
import pytest

from nucleation_md import FerreiraModel
from nucleation_md.sources.continuum import ContinuumSource

EXP_DIR = Path(__file__).resolve().parents[1] / "data" / "experimental"


@pytest.fixture(scope="module")
def src():
    return ContinuumSource()


@pytest.mark.parametrize("mode", ["closure", "integrated", "constant"])
def test_input_signs_on_continuum_rows(mode):
    """Inputs obey the convention in every ΔS_V mode of the source."""
    src = ContinuumSource(dSv_mode=mode)
    fm = FerreiraModel()
    for i in src.rows:
        k = fm.coefficients(src.state(i))
        assert k.DSv < 0.0, i
        assert k.DT > 0.0, i
        assert k.dG_V < 0.0, i
        assert k.dG_S > 0.0, i


def test_exponent_sign_is_explicit(src):
    ref = src.reference_state()
    s = src.state(38)
    plus = FerreiraModel(reference=ref, exponent_sign=+1).rate(s, 1e-9, 2.5e-10, 1e28)
    minus = FerreiraModel(reference=ref, exponent_sign=-1).rate(s, 1e-9, 2.5e-10, 1e28)
    fm = FerreiraModel(reference=ref)
    assert plus / minus == pytest.approx(np.exp(2 * fm.barrier(s) / fm.dG_eq), rel=1e-12)


@pytest.mark.skipif(not EXP_DIR.exists() or not any(EXP_DIR.glob("rate_trend_*.csv")),
                    reason="experimental data (Mendes 2023, Marques 2025) not in data/experimental/ "
                           "and the experimental observable for the I-trend is not yet defined")
def test_rate_trend_matches_experiment():
    pytest.fail("define the experimental observable and its expected trend (see docs/CONVENTIONS.md)")
