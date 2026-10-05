# CALPHAD — T_L, T_S e entalpia de fusão das ligas (COST507-modified.tdb, pycalphad)

Base: `C:\Users\ileao\OneDrive\Documentos\OpenCalphad\OC6\COST507-modified.tdb` (só leitura). Equilíbrio com todas as fases da base para os elementos de cada liga.

Definições de ΔH (a escolha é do autor):
- **ΔH_a** = H_LIQ(x₀) − H_FCC(x₀) a T_L(FCC_A1): líquido único × CFC única da composição nominal (análogo ao ΔH_m de metal puro; é o que a MD de Al puro mede).
- **ΔH_b** = H_eq(T_L) − H_eq(T_S): calor total liberado na solidificação de equilíbrio (latente + sensível).
- **ΔH_c** = ΔH_b − ∫ c_p,sólido dT entre T_S e T_L: só a parte latente (c_p do conjunto sólido de equilíbrio logo abaixo de T_S).

| liga | T_L [K] | 1º sólido | T_L(FCC_A1) [K] | T_S [K] | ΔH_a (x₀ / x_L) [kJ/kg] | ΔH_b [kJ/kg] | ΔH_c [kJ/kg] | T_L / ΔH do script |
|---|---|---|---|---|---|---|---|---|
| Al puro | 933.46 | FCC_A1 | 933.46 | 933.46 | 397.0 (397.0) | 397.1 | 397.1 | — |
| Al-0,8Si-0,6Mg-0,2Fe | 925.15 | FCC_A1 | 925.14 | 873.60 | 385.2 (385.2) | 452.3 | 390.3 | 925.32 / 335.3 |
| Al-3Cu-5Nb-0,1Fe | 1598.66 | AL3M_D022 | 925.36 | 856.50 | 386.3 (381.1) | 424.2 | 346.2 | 936.58 / 335.3 |

Composições e fases:
- Al puro: M = 26.9815 g/mol; x = AL 1.00000; líquido a T_L(FCC_A1): ; sólidos logo abaixo de T_S: FCC_A1 1.0000
- Al-0,8Si-0,6Mg-0,2Fe: M = 27.0001 g/mol; x = AL 0.98468, SI 0.00769, MG 0.00667, FE 0.00097; líquido a T_L(FCC_A1): FE 0.2000 wt.%, MG 0.6000 wt.%, SI 0.8000 wt.%; sólidos logo abaixo de T_S: ALFESI_BETA 0.0059, FCC_A1 0.9941
- Al-3Cu-5Nb-0,1Fe: M = 28.4993 g/mol; x = AL 0.97070, CU 0.01345, NB 0.01534, FE 0.00051; líquido a T_L(FCC_A1): CU 3.3094 wt.%, FE 0.1103 wt.%, NB 0.0041 wt.%; sólidos logo abaixo de T_S: FCC_A1 0.9401, AL13FE4 0.0018, AL3M_D022 0.0581
