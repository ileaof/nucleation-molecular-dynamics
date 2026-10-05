import numpy as np
import pytest

from nucleation_md import Quantity
from nucleation_md.geometry import f_theta, df_dtheta, df_legacy, cap_area, cap_volume


def test_constant_quantity():
    q = Quantity.constant(2.5)
    assert q(1e-9) == 2.5
    assert q.d(1e-9) == 0.0
    assert np.allclose(q(np.array([1.0, 2.0])), 2.5)


def test_numeric_derivative():
    q = Quantity(lambda r: np.sin(r / 1e-9))
    r = 0.7e-9
    assert q.d(r) == pytest.approx(np.cos(r / 1e-9) / 1e-9, rel=1e-8)


def test_tabulated_reproduces_cubic():
    r = np.linspace(1e-9, 5e-9, 9)
    y = 3.0 + 2e9 * r - 4e18 * r ** 2 + 1e26 * r ** 3
    q = Quantity.tabulated(r[::-1], y[::-1])       # unsorted input on purpose
    x = 2.3e-9
    assert q(x) == pytest.approx(3.0 + 2e9 * x - 4e18 * x ** 2 + 1e26 * x ** 3, rel=1e-10)
    assert q.d(x) == pytest.approx(2e9 - 8e18 * x + 3e26 * x ** 2, rel=1e-8)


def test_tabulated_rejects_duplicates():
    with pytest.raises(ValueError):
        Quantity.tabulated([1.0, 1.0, 2.0, 3.0], [0, 1, 2, 3])


@pytest.mark.parametrize("theta", [0.3, 1.0, np.pi / 2, 2.5, 3.06])
def test_df_dtheta(theta):
    h = 1e-6
    num = (f_theta(theta + h) - f_theta(theta - h)) / (2 * h)
    assert df_dtheta(theta) == pytest.approx(num, rel=1e-6, abs=1e-9)


def test_sphere_limits():
    r = 2e-9
    assert f_theta(np.pi) == pytest.approx(4.0)
    assert cap_area(r, np.pi) == pytest.approx(4 * np.pi * r ** 2)
    assert cap_volume(r, np.pi) == pytest.approx(4 / 3 * np.pi * r ** 3)
    assert df_legacy(np.pi) == pytest.approx(-16.0)   # value seen in the continuum table
