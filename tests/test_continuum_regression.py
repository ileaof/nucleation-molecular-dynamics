"""
Regression against the continuum tables (tablefull_paper_B_MRS_2026.csv, printed with %g,
i.e. ~6 significant digits -> rtol 1e-5 / 1e-4).

The `characterization` tests document how the continuum script builds its solution;
see docs/CONTINUUM_VERIFICATION.md.
"""
import numpy as np
import pytest

from nucleation_md import FerreiraModel
from nucleation_md.geometry import f_theta, df_legacy
from nucleation_md.sources.continuum import ContinuumSource, delta_Sv_constant
from nucleation_md.surface import metric_2sphere


@pytest.fixture(scope="module")
def src():
    """ΔS_V as tabulated (DSv_hom, the value the script uses for GB)."""
    return ContinuumSource(dSv_mode="integrated")


def test_table_rows(src):
    assert src.rows == list(range(1, 72))
    assert np.isclose(src.experimental["i"].iloc[0], 38.5)


def test_gamma_matches_table(src):
    t = src.table
    r = t["r_hom_2nd"].to_numpy()
    assert np.allclose(src.gamma(r), t["gam_hom"], rtol=1e-5)
    assert np.allclose(src.gamma.d(r), t["dfgamdr_ana"], rtol=1e-5)


def test_gamma_analytic_derivative_close_to_numeric(src):
    """`dgammadr_func_r` (L373) is not the exact derivative of `gamma_func_r`: its α carries σ₀;
    the difference is ~3e-6 relative over the table."""
    r = src.table["r_hom_2nd"].to_numpy()
    h = 1e-6 * r
    num = (src.gamma(r + h) - src.gamma(r - h)) / (2 * h)
    assert np.allclose(src.gamma.d(r), num, rtol=1e-5)


@pytest.mark.parametrize("i", [1, 10, 38, 71])
def test_decomposition_matches_GB_GS_GC(src, i):
    """GB = ΔS_V·ΔT, GS = γ, GC = γ/f·dfthetadr(θ) (legacy ∂f/∂r)."""
    row = src.row(i)
    k = FerreiraModel().coefficients(src.state(i, df_mode="continuum_legacy"))
    assert k.dG_V == pytest.approx(row["GB"], rel=1e-4)
    assert k.dG_S == pytest.approx(row["GS"], rel=1e-5)
    assert k.dG_conf == pytest.approx(row["GC"], rel=1e-4)


def test_spline_fields_reproduce_rows(src):
    t = src.table
    r = t["r_hom_2nd"].to_numpy()
    assert np.allclose(src.DT(r), t["DT"], rtol=1e-12)
    assert np.allclose(src.DSv(r), t["DSv_hom"], rtol=1e-12)
    assert np.allclose(src.theta(r), t["theta_2nd"], rtol=1e-12)
    assert np.allclose(df_legacy(t["theta_2nd"]), t["dfthetadr_2nd"], rtol=1e-4)


@pytest.mark.parametrize("i", [1, 2, 20, 38, 70, 71])
def test_closure_makes_second_order_radius_equal_r(i):
    """Continuum design (L902, L905): with the closure ∂ΔS_V/∂r and a homogeneous nucleus,
    FerreiraModel's -3b/(2a) returns the loop radius r."""
    src = ContinuumSource(dSv_mode="closure", theta_mode="homogeneous")
    s = src.state(i)
    assert FerreiraModel(order=2).r_critical(s) == pytest.approx(s.r_eval, rel=1e-10)


@pytest.mark.characterization
def test_characterization_DT_is_cnt_inverse(src):
    """Tabulated ΔT is Gibbs–Thomson, ΔT = 2Γ/r with Γ = -γ(r)/ΔS_V (value returned at L890)."""
    t = src.table
    r = t["r_hom_2nd"].to_numpy()
    assert np.allclose(t["DT"], -2 * src.gamma(r) / (r * delta_Sv_constant()), rtol=1e-5)


@pytest.mark.characterization
def test_characterization_second_order_radius_equals_loop_radius(src):
    t = src.table
    assert np.allclose(t["r_hom_2nd"], t["r_hom_1st"], rtol=1e-6)


def test_metric_form_equals_gamma0_at_reference():
    g = metric_2sphere(0.296, 6.3e-6)
    assert g(6.3e-6) == pytest.approx(0.296, rel=1e-12)
    # 3/(1 + 2x²) with x = r/r0 when θ = θ0 and no stress term
    assert g(3.15e-6) == pytest.approx(0.296 * 3 / 1.5, rel=1e-12)


def test_metric_form_stress_integral_constant_sigma():
    """Σ constant: (1/(4πr0²))·Σ·(4πr² - 4πr0²) = Σ·(r²/r0² - 1)."""
    r0, S = 6.3e-6, 0.5
    g = metric_2sphere(0.296, r0, Sigma=lambda r: S)
    r = 5e-6
    expect = 0.296 * 3 / (1 + 2 * (r / r0) ** 2) - S * ((r / r0) ** 2 - 1)
    assert g(r) == pytest.approx(expect, rel=1e-9)


# ----------------------------------------------- thermal-field tensor Γ = A·∇T

@pytest.fixture(scope="module")
def src_hom():
    """Closure ∂ΔS_V/∂r and homogeneous nucleus: the continuum's own construction."""
    return ContinuumSource(dSv_mode="closure", theta_mode="homogeneous")


def test_gamma_tensors_match_GT_columns(src_hom):
    """rtol 3e-5: ΔT and GT are both printed with 6 digits (max deviation found 1.0e-5)."""
    fm = FerreiraModel()
    for i in src_hom.rows:
        s, row = src_hom.state(i), src_hom.row(i)
        assert fm.gamma_tensor(s, 1) == pytest.approx(row["GT_het"], rel=3e-5), i
        assert fm.gamma_tensor(s, 2) == pytest.approx(row["GT_het_2nd"], rel=3e-5), i


def test_gradient_column_is_gamma_over_area(src_hom):
    fm = FerreiraModel(order=2)
    for i in src_hom.rows:
        assert fm.grad_T(src_hom.state(i)) == pytest.approx(src_hom.row(i)["dTdr"], rel=3e-5), i


def test_experimental_gradient_point_slide25(src_hom):
    """∇T = 3294 K/m (Mendes et al. 2023): Γ¹ˢᵗ = 3.887e-7, Γ²ⁿᵈ = 1.06852e-6 m·K (slide 25)."""
    fm = FerreiraModel()
    s = src_hom.state_at_gradient(3294.0)
    exp = src_hom.experimental.iloc[0]
    assert s.r_eval == pytest.approx(exp["r_hom_2nd"], rel=1e-4)
    assert fm.gamma_tensor(s, 1) == pytest.approx(3.887e-7, rel=2e-3)
    assert fm.gamma_tensor(s, 2) == pytest.approx(1.06852e-6, rel=2e-3)
    assert fm.grad_T(s, 2) == pytest.approx(3294.0, rel=2e-3)
