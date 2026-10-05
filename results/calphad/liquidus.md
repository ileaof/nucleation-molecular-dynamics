# CALPHAD — liquidus para a MD

## Binário Al-Nb — ALFENB-B2-2SL.TDB

Fases da base para Al-Nb: AL3NB, AL8FE5, ALNB2, ALNB3, B2_2SL, FCC_A1, HCP_A3, LIQUID, L_C14, MU_4SL

| composição | T_L estável [K] | 1º sólido | T_L(FCC_A1) metaestável [K] | x_Nb(FCC) | k | ΔH_m [J/mol] | ΔH_m [J/kg] |
|---|---|---|---|---|---|---|---|
| Al puro | 933.47 | FCC_A1 | 933.47 | — | — | — | — |
| Al-5.00wt%Nb (x=0.01506) | 1539.2 | AL3NB | 1059.68 | 0.07274 | 4.869 | 15028 | 537374 |
| Al-5.16wt%Nb (x=0.01556) | 1544.2 | AL3NB | 1061.70 | 0.07354 | 4.765 | 15041 | 537192 |
| x_Nb = 0.03468 (script L723) | 1675.1 | AL3NB | 1121.85 | 0.09670 | 2.813 | 14989 | 512498 |

## Liga nominal Al-3Cu-5Nb-0,1Fe — COST507-modified.tdb (só LIQUID + FCC_A1)

A COST507 não tem intermetálicos de Nb; o resultado é o liquidus metaestável da FCC_A1.

| composição | T_L(FCC_A1) [K] |
|---|---|
| Al-3Cu-5Nb-0,1Fe | 1155.80 |
| Al-5Nb (binário, COST507) | 1157.73 |
| Al-3Cu (binário, COST507) | 924.80 |

Referência DTA (script de nucleação, L218): T_L = 936,58 K.
## Caminho de solidificação estável, Al-5wt%Nb (ALFENB-B2-2SL.TDB)

| T [K] | fases (fração; x_Nb) |
|---|---|
| 1539,2 | liquidus: primeiro sólido = Al₃Nb |
| 1400 | LIQUID 0,961 (x_Nb 5,4e-3); Al₃Nb 0,039 |
| 1200 | LIQUID 0,943 (x_Nb 8,1e-4); Al₃Nb 0,057 |
| 1000 | LIQUID 0,940 (x_Nb 5,6e-5); Al₃Nb 0,060 |
| 935,0 | LIQUID 0,940 (x_Nb 1,8e-5); Al₃Nb 0,060 |
| **934,52** | **FCC_A1 aparece** (peritético L + Al₃Nb → FCC_A1); FCC com x_Nb ≈ 1,5e-3 |
| 900 | FCC_A1 0,944 (x_Nb 1,0e-3); Al₃Nb 0,056 |

O α-Al (FCC_A1) se forma a partir de um líquido com x_Nb ≈ 1,8e-5, praticamente Al puro em Nb,
a 934,52 K, 2 K abaixo do T_L do DTA (936,58 K). O Nb está quase todo no Al₃Nb (~6 % em fração molar).
Diagramas: AlNb_full.png, AlNb_Al_rich.png.
