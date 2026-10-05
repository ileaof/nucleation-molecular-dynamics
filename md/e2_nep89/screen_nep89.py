# Static screening of NEP89 for the Al-Nb system: fcc Al, bcc Nb, Al3Nb (D0_22, L1_2) at 0 K.
import numpy as np
from ase.build import bulk
from ase import Atoms
from ase.filters import ExpCellFilter
from ase.optimize import BFGS
from calorine.calculators import CPUNEP

POT = "../../potentials/NEP89_20250409/nep89_20250409.txt"

def relax(atoms, mask=None):
    atoms.calc = CPUNEP(POT)
    BFGS(ExpCellFilter(atoms, mask=mask), logfile=None).run(fmax=1e-4, steps=500)
    return atoms.get_potential_energy() / len(atoms), atoms

eAl, al = relax(bulk("Al", "fcc", a=4.05, cubic=True))
eNb, nb = relax(bulk("Nb", "bcc", a=3.30, cubic=True))
a, c = 3.841, 8.609
frac = [(0,0,0,"Nb"),(.5,.5,.5,"Nb"),(0,0,.5,"Al"),(.5,.5,0,"Al"),(0,.5,.25,"Al"),(.5,0,.25,"Al"),(0,.5,.75,"Al"),(.5,0,.75,"Al")]
d022 = Atoms([f[3] for f in frac], scaled_positions=[f[:3] for f in frac], cell=[a, a, c], pbc=True)
l12 = Atoms(["Nb","Al","Al","Al"], scaled_positions=[(0,0,0),(.5,.5,0),(.5,0,.5),(0,.5,.5)], cell=[4.0]*3, pbc=True)
d022_fixed = d022.copy(); d022_fixed.calc = CPUNEP(POT); e_fixed = d022_fixed.get_potential_energy() / 8
eD, d = relax(d022)
eL, l = relax(l12)
ref = 0.75 * eAl + 0.25 * eNb
print(f"fcc Al: a = {al.cell[0,0]:.4f} Å, E = {eAl:.4f} eV/atom   (exp a = 4.05)")
print(f"bcc Nb: a = {nb.cell[0,0]:.4f} Å, E = {eNb:.4f} eV/atom   (exp a = 3.30)")
print(f"Al3Nb D0_22 relaxed: a = {d.cell[0,0]:.4f}, c = {d.cell[2,2]:.4f} Å, ΔH_f = {eD-ref:+.4f} eV/atom  (exp a=3.841, c=8.609; ΔH_f ≈ -0.4)")
print(f"Al3Nb D0_22 at exp lattice: ΔH_f = {e_fixed-ref:+.4f} eV/atom")
print(f"Al3Nb L1_2 relaxed: a = {l.cell[0,0]:.4f} Å, ΔH_f = {eL-ref:+.4f} eV/atom")
