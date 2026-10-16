# Difusão no líquido Al-3 wt% Cu, NEP89 (GPUMD, 16 384 át., 212 Cu, NPT P = 0, 200 ps)

| T [K] | D_Cu MD [m²/s] | D_Cu script [m²/s] | MD/script | D_Al MD [m²/s] |
|---|---|---|---|---|
| 900.00 | 3.985e-09 ± 1.6e-11 | 4.290e-09 | 0.93 | 4.611e-09 ± 1.2e-12 |
| 936.58 | 4.270e-09 ± 6.7e-13 | 4.862e-09 | 0.88 | 5.188e-09 ± 4.4e-12 |
| 1000.00 | 5.235e-09 ± 8.4e-12 | 5.912e-09 | 0.89 | 6.118e-09 ± 1.2e-11 |
| 1100.00 | 7.005e-09 ± 9.9e-11 | 7.685e-09 | 0.91 | 7.731e-09 ± 9.1e-12 |
| 1200.00 | 8.601e-09 ± 2.4e-11 | 9.564e-09 | 0.90 | 9.396e-09 ± 1.7e-12 |

Arrhenius MD: Cu D0 = 9.383e-08 m²/s, Q = 23855 J/mol; Al D0 = 7.868e-08 m²/s, Q = 21212 J/mol
Script (Cu): D0 = 1.060e-07 m²/s, Q = 24000 J/mol

## SDAS Al-3Cu-5Nb-0,1Fe: D_L do Cu do script × da MD (T_L = 936.58 K: script 4.862e-09, MD 4.385e-09 m²/s)

| curva | D_L do Cu | SDAS [μm] | dentro das barras | bias_log |
|---|---|---|---|---|
| MRB Γ¹ˢᵗ contínuo | script | 7.67·t^0.377 | 100% | -0.000 |
| MRB Γ¹ˢᵗ contínuo | MD | 7.54·t^0.376 | 100% | -0.019 |
| MRB Γ²ⁿᵈ contínuo | script | 9.38·t^0.379 | 94% | +0.213 |
| MRB Γ²ⁿᵈ contínuo | MD | 9.22·t^0.379 | 94% | +0.194 |
| RB Γ²ⁿᵈ contínuo | script | 27.59·t^0.333 | 0% | +1.074 |
| RB Γ²ⁿᵈ contínuo | MD | 26.76·t^0.333 | 0% | +1.044 |
