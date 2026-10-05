# Capillary fluctuation method (Hoyt, Asta, Karma 2001): quasi-2D solid-liquid ribbon at T_m.
# Interface normal z, fluctuation direction x, thin y.  orient: 100 | 110 | 111 (interface plane)
#   mpirun -np 4 lmp -in in.cfm -var pot borovikov -var Tm 930 -var orient 100
variable pot    index borovikov
variable Tm     index 930
variable orient index 100
variable ny     index 3   # y cells (thickness); E2 used 3 | 2 | 3
variable dump   index cfm_${pot}_${orient}.dump
variable seed   index 3301
variable neq    index 25000
variable nprod  index 200000
variable nstage index 10000
variable ndump  index 2000
variable a0 equal 4.05
if "${pot} == borovikov" then "variable a0 equal 4.016"
if "${pot} == mendelev" then "variable a0 equal 4.045"
units        metal
atom_style   atomic
boundary     p p p
# (x, y, z) = (fluctuation, thin, interface normal)
if "${orient} == 100" then "lattice fcc ${a0} orient x 1 0 0 orient y 0 1 0 orient z 0 0 1" "region box block 0 60 0 ${ny} 0 30"
if "${orient} == 110" then "lattice fcc ${a0} orient x 0 0 1 orient y 1 -1 0 orient z 1 1 0" "region box block 0 60 0 ${ny} 0 21"
if "${orient} == 111" then "lattice fcc ${a0} orient x 1 -1 0 orient y 1 1 -2 orient z 1 1 1" "region box block 0 42 0 ${ny} 0 17"
create_box   1 box
create_atoms 1 box
mass         1 26.982
include      ../e2_melting/pot_${pot}.lmp
timestep     0.002
thermo       1000
# solid at Tm, lateral lattice fixed afterwards
velocity     all create ${Tm} ${seed} mom yes rot yes dist gaussian
fix          eq all npt temp ${Tm} ${Tm} 0.2 aniso 0 0 2.0
run          ${nstage}
unfix        eq
variable     zmid equal (zlo+zhi)/2
region       up block INF INF INF INF ${zmid} INF units box
group        liq region up
group        sol subtract all liq
velocity     sol set 0 0 0
fix          fr sol setforce 0 0 0
fix          m liq nvt temp 1500 1500 0.2
run          ${nstage}
unfix        m
fix          c liq nvt temp ${Tm} ${Tm} 0.2
run          ${nstage}
unfix        c
unfix        fr
velocity     sol create ${Tm} $(v_seed+1) mom yes rot yes dist gaussian
fix          rel all npt temp ${Tm} ${Tm} 0.2 z 0 0 2.0
run          ${neq}
unfix        rel
# production: NPH along z (interface position drifts slowly if T != T_m)
compute      cna all cna/atom $(0.854*v_a0)
variable     issol atom c_cna==1
compute      ns all reduce sum v_issol
variable     fsol equal c_ns/atoms
fix          prod all nph z 0 0 2.0
thermo_style custom step temp press pzz pe lz v_fsol
dump         d all custom ${ndump} ${dump} id x y z c_cna
dump_modify  d sort id
run          ${nprod}
