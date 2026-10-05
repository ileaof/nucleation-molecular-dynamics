"""SDAS model port (nucleation_md.sdas) and the second alloy (Al08Si06Mg02Fe)."""
import contextlib
import io
import sys
from pathlib import Path

import numpy as np
import pytest

from nucleation_md import FerreiraModel
from nucleation_md.sdas import ALLOY_SDAS, diff_length_scale, sdas_mrb, sdas_rb, sdas_from_gamma, power_law_fit
from nucleation_md.sources.continuum import ContinuumSource

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))


@pytest.fixture(scope="module")
def script_ns():
    from run_continuum import sdas_namespace
    with contextlib.redirect_stdout(io.StringIO()):
        return sdas_namespace("Al08Si06Mg02Fe")


def test_port_reproduces_sdas_script(script_ns):
    system = ALLOY_SDAS["Al08Si06Mg02Fe"]
    tSL = system.tSL(system.P)
    assert np.allclose(tSL, script_ns["tSL"], rtol=1e-14)
    _, kef, w = diff_length_scale(system, system.VL(system.P))
    G2, G1 = script_ns["Gibbs_Thomson_het"], script_ns["Gibbs_Thomson_het_1st"]
    assert np.allclose(sdas_mrb(G2, system, kef, w, tSL)[0], script_ns["SDAS_Eq"], rtol=1e-12)
    assert np.allclose(sdas_mrb(G1, system, kef, w, tSL)[0], script_ns["SDAS_het_1st_Eq"], rtol=1e-12)
    assert np.allclose(sdas_rb(G2, system, tSL), script_ns["SDAS_RB"], rtol=1e-12)


@pytest.mark.parametrize("alloy, Gamma, a, b", [
    ("Al08Si06Mg02Fe", 1.57452e-06, 6.71, 0.387),   # slide 30
    ("Al08Si06Mg02Fe", 7.63702e-07, 5.81, 0.387),
    ("Al3Cu5Nb01Fe", 1.06852e-06, 9.38, 0.379),     # slide 27
    ("Al3Cu5Nb01Fe", 3.887e-07, 7.67, 0.377),
])
def test_power_law_fits_of_the_slides(alloy, Gamma, a, b):
    fa, fb = power_law_fit(*sdas_from_gamma(Gamma, ALLOY_SDAS[alloy]))
    assert fa == pytest.approx(a, abs=0.006)
    assert fb == pytest.approx(b, abs=0.0006)


@pytest.fixture(scope="module")
def src2():
    return ContinuumSource(alloy="Al08Si06Mg02Fe", dSv_mode="closure", theta_mode="homogeneous")


def test_second_alloy_tensor_columns(src2):
    """rtol 3e-5: the table is printed with 6 digits."""
    fm = FerreiraModel()
    for i in src2.rows:
        s, row = src2.state(i), src2.row(i)
        assert fm.r_critical(s, 2) == pytest.approx(s.r_eval, rel=1e-10), i
        assert fm.gamma_tensor(s, 1) == pytest.approx(row["GT_het"], rel=3e-5), i
        assert fm.gamma_tensor(s, 2) == pytest.approx(row["GT_het_2nd"], rel=3e-5), i
        assert fm.grad_T(s, 2) == pytest.approx(row["dTdr"], rel=3e-5), i


def test_second_alloy_experimental_gradient_slide29(src2):
    """∇T = 7902.38 K/m (Marques et al. 2025): Γ¹ˢᵗ = 7.63702e-7, Γ²ⁿᵈ = 1.57452e-6 m·K (script L881-)."""
    fm = FerreiraModel()
    s = src2.state_at_gradient(7902.38)
    assert s.r_eval == pytest.approx(3.9821e-06, rel=1e-4)
    assert fm.gamma_tensor(s, 1) == pytest.approx(7.63702e-07, rel=1e-4)
    assert fm.gamma_tensor(s, 2) == pytest.approx(1.57452e-06, rel=1e-4)
