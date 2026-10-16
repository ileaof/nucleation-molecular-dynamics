# Self-diffusion of liquid Al (Borovikov 2024) from the mean-squared displacement, NPT at T, P = 0.
#   lmp -sf gpu -pk gpu 1 -in in.msd -var T 813
variable pot   index borovikov
variable T     index 813
variable n     index 12          # 6912 atoms
variable seed  index 11
variable nprod index 50000       # 100 ps
units        metal
atom_style   atomic
boundary     p p p
lattice      fcc 4.016
region       box block 0 ${n} 0 ${n} 0 ${n}
create_box   1 box
create_atoms 1 box
mass         1 26.982
include      ../e2_melting/pot_${pot}.lmp
timestep     0.002
thermo       2500
velocity     all create 1500 ${seed} mom yes rot yes dist gaussian
fix          m all npt temp 1500 1500 0.2 iso 0.0 0.0 2.0
run          5000
unfix        m
fix          c all npt temp ${T} ${T} 0.2 iso 0.0 0.0 2.0
run          10000
reset_timestep 0
compute      msd all msd com yes
variable     t equal step*dt
fix          out all ave/time 50 1 50 v_t c_msd[4] file msd_${T}.txt
run          ${nprod}
