import pytest

from nucleation_md.validation import potential_gate, PotentialNotValidated

REF = {"T_m": 936.58, "dH_m": 335300.0, "gamma_0": 0.296}
TOL = {"T_m": 0.02, "dH_m": 0.10, "gamma_0": 0.20}


def test_gate_passes():
    md = {"T_m": 930.0, "dH_m": 320000.0, "gamma_0": 0.27}
    assert potential_gate(md, REF, TOL).passed


def test_gate_blocks():
    md = {"T_m": 870.0, "dH_m": 320000.0, "gamma_0": 0.27}
    with pytest.raises(PotentialNotValidated):
        potential_gate(md, REF, TOL)


def test_gate_requires_explicit_tolerances():
    with pytest.raises(KeyError):
        potential_gate(REF, REF, {"T_m": 0.01})
