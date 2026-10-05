# -*- coding: utf-8 -*-
"""
Gate: the interatomic potential must reproduce the reference T_m, ΔH_m and γ₀
before any model comparison (stage E2). Tolerances have no defaults: they are
part of the validation protocol and must be stated explicitly.
"""
from dataclasses import dataclass


class PotentialNotValidated(RuntimeError):
    pass


@dataclass
class GateResult:
    passed: bool
    rows: list   # (key, measured, reference, rel_err, rtol, ok)

    def table(self):
        lines = [f"{'quantity':<10}{'MD':>14}{'reference':>14}{'rel.err':>10}{'rtol':>8}  ok"]
        for k, m, r, e, t, ok in self.rows:
            lines.append(f"{k:<10}{m:>14.6g}{r:>14.6g}{e:>10.3%}{t:>8.1%}  {'yes' if ok else 'NO'}")
        return "\n".join(lines)


REQUIRED = ("T_m", "dH_m", "gamma_0")


def potential_gate(measured, reference, rtol, raise_on_fail=True):
    """
    measured, reference - dicts with at least T_m [K], dH_m [J/kg or J/m3, same unit], gamma_0 [J/m2]
    rtol                - dict of relative tolerances, one per key
    """
    rows = []
    for key in REQUIRED:
        for d, name in ((measured, "measured"), (reference, "reference"), (rtol, "rtol")):
            if key not in d:
                raise KeyError(f"{name} lacks {key!r}")
    for key in reference:
        m, r, t = float(measured[key]), float(reference[key]), float(rtol[key])
        err = abs(m - r) / abs(r)
        rows.append((key, m, r, err, t, err <= t))
    result = GateResult(all(row[-1] for row in rows), rows)
    if raise_on_fail and not result.passed:
        raise PotentialNotValidated("potential fails the reference gate:\n" + result.table())
    return result
