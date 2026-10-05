"""
Model tests. Numerical inputs are the continuum script values (CONTINUUM_PARAMS) or
smooth synthetic fields; they test algebra, not physics.
"""
import numpy as np
import pytest

from nucleation_md import CNT, Tolman, FerreiraModel, Quantity, State, constant_state
from nucleation_md.sources.continuum import CONTINUUM_PARAMS, delta_Sv_constant

P = CONTINUUM_PARAMS
DSV = delta_Sv_constant()
GAMMA0 = P["gamma_0"]
DT0 = 0.1


def synthetic_state(df_mode="chain", r_eval=5e-6):
    """All four fields depend on r (smooth, analytic) to exercise every term of a and b."""
    r0 = P["req"]
    DSv = Quantity(lambda r: DSV * (1 + 0.3 * (r / r0 - 1) ** 2))
    DT = Quantity(lambda r: DT0 * (r0 / r) ** 0.5)
    gamma = Quantity(lambda r: GAMMA0 * (r0 / r) ** 2 - 0.05 * np.log(r0 / r))
    theta = Quantity(lambda r: 2.0 + 0.4 * np.tanh(r / r0))
    return State(DSv, DT, gamma, theta, r_eval=r_eval, df_mode=df_mode)


# ------------------------------------------------------------------ CNT limit

@pytest.mark.parametrize("theta", [np.pi, 2.0, 1.2])
def test_ferreira_reproduces_cnt_without_surface_potential_and_configurational_terms(theta):
    """Mandatory: ∂ΔG_S/∂r = 0 and ΔG_C = 0 (constant γ, constant θ) -> CNT exactly."""
    s = constant_state(DSV, DT0, GAMMA0, theta=theta, r_eval=3e-6)
    fm, cnt = FerreiraModel(order=1), CNT()
    k = fm.coefficients(s)
    assert k.ddG_S_dr == 0.0 and k.dG_conf == 0.0
    assert fm.r_critical(s, order=1) == pytest.approx(cnt.r_critical(s), rel=1e-14)
    assert fm.barrier(s, order=1) == pytest.approx(cnt.barrier(s), rel=1e-14)
    assert fm.gamma_tensor(s, order=1) == pytest.approx(cnt.gamma_tensor(s), rel=1e-14)
    for r in (1e-7, 1e-6, 1e-5):
        assert fm.delta_G(r, s) == pytest.approx(cnt.delta_G(r, s), rel=1e-14)
    # with every field constant a = 0: the parabola degenerates and r_c² is undefined
    assert k.a == 0.0
    assert fm.r_critical(s, order=2) == np.inf


def test_cnt_closed_forms():
    s = constant_state(DSV, DT0, GAMMA0)
    g = DSV * DT0
    assert CNT().r_critical(s) == pytest.approx(-2 * GAMMA0 / g)
    assert CNT().barrier(s) == pytest.approx(16 * np.pi * GAMMA0 ** 3 / (3 * g ** 2))
    assert CNT().gamma_tensor(s) == pytest.approx(-GAMMA0 / DSV)


# ------------------------------------------------- the parabola follows from ΔG(r)

@pytest.mark.parametrize("r", [2e-6, 5e-6, 8e-6])
def test_dG_dr_equals_parabola(r):
    """dΔG/dr = (πr/3)(a r² + 3b r + 6γf) with a, b exactly as defined."""
    s, fm = synthetic_state(), FerreiraModel()
    h = 1e-6 * r
    num = (fm.delta_G(r + h, s) - fm.delta_G(r - h, s)) / (2 * h)
    assert fm.dG_dr(r, s) == pytest.approx(num, rel=1e-7)


def test_denominator_decomposition():
    s, fm = synthetic_state(), FerreiraModel()
    k = fm.coefficients(s)
    assert k.b / k.f == pytest.approx(k.dG_V + k.ddG_S_dr + k.dG_conf, rel=1e-14)
    assert fm.r_critical(s, order=1) == pytest.approx(-2 * k.dG_S / k.denominator, rel=1e-14)


def test_second_order_radius_is_vertex_definition():
    """r_c² stays -3b/(2a) (user's definition), not a root of the quadratic."""
    s, fm = synthetic_state(), FerreiraModel()
    k = fm.coefficients(s)
    assert fm.r_critical(s, order=2) == -1.5 * k.b / k.a


def test_exact_roots_solve_quadratic():
    s, fm = synthetic_state(), FerreiraModel()
    k = fm.coefficients(s)
    r_minus, r_plus, disc = fm.exact_roots(s)
    for root in (r_minus, r_plus):
        res = k.a * root ** 2 + 3 * k.b * root + k.c
        assert abs(res) <= 1e-9 * (abs(k.a * root ** 2) + abs(3 * k.b * root) + abs(k.c))
    # vertex = mean of the two roots (real or complex conjugate)
    assert np.real((r_minus + r_plus) / 2) == pytest.approx(fm.r_critical(s, order=2), rel=1e-12)


def test_exact_root_tends_to_first_order_when_a_vanishes():
    s = constant_state(DSV, DT0, GAMMA0, r_eval=3e-6)
    s.DSv = Quantity(lambda r: DSV * (1 + 1e-9 * r / P["req"]))     # a -> 0
    fm = FerreiraModel()
    r_minus, r_plus, _ = fm.exact_roots(s)
    finite = r_plus if abs(r_plus) < abs(r_minus) else r_minus
    assert np.real(finite) == pytest.approx(fm.r_critical(s, order=1), rel=1e-6)


def test_self_consistent_radius():
    s, fm = synthetic_state(), FerreiraModel(order=1)
    rc = fm.r_critical(s, order=1, solve="self_consistent", bracket=(1e-8, 5e-5))
    k = fm.coefficients(s, rc)
    assert rc == pytest.approx(-2 * k.dG_S * k.f / k.b, rel=1e-10)


# --------------------------------------------------------------------------- Γ

def test_first_order_tensor_is_gibbs_thomson_when_DT_uniform():
    """With ∂ΔT/∂r = 0:  Γ¹ˢᵗ = ΔT·r_c¹/2."""
    base = synthetic_state()
    s = State(base.DSv, Quantity.constant(DT0), base.gamma, base.theta, r_eval=base.r_eval)
    fm = FerreiraModel()
    assert fm.gamma_tensor(s, order=1) == pytest.approx(DT0 * fm.r_critical(s, order=1) / 2, rel=1e-12)


@pytest.mark.parametrize("order", [1, 2])
def test_tensor_is_gibbs_thomson_of_each_order(order):
    """Slide 10 (MRS 2026): Γ = ΔT·r_c/2 at both orders, with ∂ΔT/∂r ≠ 0 too."""
    s, fm = synthetic_state(), FerreiraModel()
    k = fm.coefficients(s)
    assert k.dDT_dr != 0.0
    assert fm.gamma_tensor(s, order) == pytest.approx(k.DT * fm.r_critical(s, order) / 2, rel=1e-12)


def test_tensor_area_times_gradient():
    """Γ = A·∇T  <=>  grad_T() = Γ/A."""
    from nucleation_md.geometry import cap_area
    s, fm = synthetic_state(), FerreiraModel()
    rc = fm.r_critical(s, 2)
    assert fm.grad_T(s, 2) * cap_area(rc, s.theta(rc)) == pytest.approx(fm.gamma_tensor(s, 2), rel=1e-14)


# ---------------------------------------------------------------------- Tolman

def test_tolman():
    s = constant_state(DSV, DT0, GAMMA0)
    assert Tolman(0.0).r_critical(s) == pytest.approx(CNT().r_critical(s))
    t = Tolman(1e-7)
    rc = t.r_critical(s)
    assert rc == pytest.approx(CNT().r_critical(s) - 2e-7)
    h = 1e-6 * rc
    assert abs(t.delta_G(rc + h, s) - t.delta_G(rc - h, s)) / (2 * h) < 1e-9 * t.barrier(s) / rc
    g = s.DSv(0) * s.DT(0)
    assert t.barrier(s) == pytest.approx(16 * np.pi * t._gamma_at(None, s) ** 3 / (3 * g ** 2), rel=1e-12)


# ------------------------------------------------------------------------ rate

def test_rate_requires_reference():
    with pytest.raises(ValueError):
        FerreiraModel().rate(synthetic_state(), D=1e-9, lam=2.5e-10, N_V=1e28)


def test_rate_at_reference_is_prefactor_times_e():
    """exp(ΔG*/ΔG*_Eq) = e at the reference state, as written in the model."""
    s = synthetic_state()
    fm = FerreiraModel(reference=s)
    from nucleation_md.geometry import cap_area
    rc = fm.r_critical(s)
    pref = 1e-9 * cap_area(rc, s.theta(rc)) / 2.5e-10 ** 4 * 1e28
    assert fm.rate(s, D=1e-9, lam=2.5e-10, N_V=1e28) == pytest.approx(pref * np.e, rel=1e-12)
