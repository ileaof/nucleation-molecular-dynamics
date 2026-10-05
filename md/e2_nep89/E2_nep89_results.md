# E2 — validação do NEP89 (Al puro), mesmo protocolo da E2 do Borovikov

Potencial: NEP89 (`potentials/NEP89_20250409/nep89_20250409.txt`; Liang et al., arXiv 2504.21286), universal, 89 elementos.
Código: GPUMD (git `eb7fbbf`, RTX 4050). Protocolos portados de `md/e2_melting` (`in.coexist`, `in.enthalpy`, `in.cfm`);
fração CFC pelo mesmo `cna/atom` do LAMMPS (cutoff 0,854·a₀, a₀ = 3,9899 Å = rede do NEP89 a 0 K) via `rerun`.
Δt = 2 fs (deriva NVE −0,0008 eV/átomo/ns; 1 fs: 0,0000 — `drift_nep89.txt`).

## T_m — coexistência isotérmica (NPT só em z, 8192 átomos, interface (001), 100 ps) — `analyze_tm_nep.py`

| T [K] | df_fcc/dt [1/ns] |
|---|---|
| 800 | +0,569 |
| 900 | +0,854 |
| 910 | +0,104 |
| 920 | +0,095 |
| 925 | −0,315 |
| 930 | −0,710 |
| 935 | −0,883 |
| 940 | −1,376 |
| 960 | −1,872 |
| 980, 1000, 1100 | sólido some por inteiro (saturação; excluídas) |

**T_m = 919 ± 3 K**: troca de sinal entre 920 e 925 K (interpolação 921,2 K); ajuste linear em ±50 K (900–960 K, 8 corridas)
916,5 ± 1,8 K. Fora da janela a cinética não é linear (800 K cresce mais devagar que 900 K), por isso excluída do ajuste.
Confirmação independente: nas fitas da flutuação capilar (NPH, que relaxa para T_m) a (111) ficou em 917–924 K.

## ΔH_m e densidades — NPT, 2048 átomos, 100 ps — `analyze_h_nep.py`

| | NEP89 (920 K) | Borovikov (935 K) | contínuo | Al puro exp. |
|---|---|---|---|---|
| ΔH_m | **377,5 kJ/kg** (10,19 kJ/mol) | 325,5 | 335,3 | 397 |
| ρ_s | 2637,2 kg/m³ | 2551,5 | 2737,04 | ≈ 2550 |
| ρ_l | 2473,2 kg/m³ | 2397,9 | 2527,05 | ≈ 2375 |
| ΔV/V_s | 6,63 % | 6,41 % | — | ≈ 6,5–7 % |
| ΔS_V = −ΔH·ρ_s/T_m | −1,083×10⁶ J m⁻³ K⁻¹ | −8,878×10⁵ | −9,799×10⁵ | — |

Erro estatístico de H (médias em blocos): ≤ 0,3 meV/átomo (≤ 0,3 % de ΔH). O líquido não cristalizou em nenhuma T
(CNA do quadro final = 0); a ΔV/V e H mostram que o sólido não fundiu (o CNA a ~T_m dá 0,39–0,64 de CFC no sólido —
a subcontagem já conhecida).

| Liga (T_s) | ρ_s NEP89 | ρ_s script | ρ_l NEP89 | ρ_l script |
|---|---|---|---|---|
| Al-0,8Si-0,6Mg-0,2Fe (798,15 K) | 2670,1 (+4,5 %) | 2555,72 | 2513,8 (+5,7 %) | 2378,33 |
| Al-3Cu-5Nb-0,1Fe (817,75 K) | 2665,2 (−2,6 %) | 2737,04 | 2506,9 (−0,8 %) | 2527,05 |

(Borovikov: 2587,4 / 2443,5 a 798,15 K e 2582,5 / 2436,9 a 817,75 K.) O NEP89 dá Al mais denso porque a rede a 0 K
é 3,99 Å (exp. 4,05).

## γ₀ — flutuação capilar (NPH, 400 ps), parâmetro de Morris (`../e2_melting/cfm_morris.py ... nep89`)

Fitas a T_m (T média da produção 924–925 K), janela capilar nos modos 4–10, 40 colunas.
O GPUMD exige espessura periódica > 2,5·(r_c + 1 Å) = 17,5 Å (`src/force/nep.cu`), por isso (100) e (110) são mais
espessas que na E2 (`in.build_cfm`); comprimento L e altura iguais aos da E2.

| Orientação | W [Å] (E2 Borovikov) | inclinação log-log (ideal −2) | γ̃ NEP89 [J/m²] | γ̃ Borovikov [J/m²] |
|---|---|---|---|---|
| (100)[100] | 20,4 (12,4) | −2,18 | 0,097 ± 0,003 | 0,094 ± 0,004 |
| (110)[001] | 23,1 (11,7) | −2,55 | 0,104 ± 0,007 | 0,096 ± 0,008 |
| (111)[1-10] | 20,0 (20,2) | −1,91 | 0,138 ± 0,002 | 0,116 ± 0,005 |

Expansão cúbica (Hoyt, Asta, Karma 2001), Monte Carlo sobre os erros estatísticos (`combine_gamma0.py`,
que reproduz o 0,106 ± 0,004 da E2 com os dados do Borovikov):

**γ₀ = 0,120 ± 0,003 J/m²**, ε₁ = 0,055 ± 0,006, ε₂ = −0,001 ± 0,002
→ γ(100) = 0,123, γ(110) = 0,120, γ(111) = 0,118 J/m².

Incertezas sistemáticas não incluídas (como na E2): 400 ps, inclinação da (110) fora de −2. Ordem provável ±15–20 %.
Referências: Al exp. ≈ 0,13–0,17 J/m²; contínuo 0,183 (Al-Si-Mg-Fe) e 0,296 (Al-Cu-Nb-Fe, r₀ = 6,3 μm).

## Resumo × contínuo e × Borovikov

| | NEP89 | Borovikov | contínuo Al-Si-Mg-Fe | contínuo Al-Cu-Nb-Fe |
|---|---|---|---|---|
| T_m [K] | 919 ± 3 | 936 ± 3 | 925,32 | 936,58 |
| ΔH_m [kJ/kg] | 377,5 | 325,5 | 335,3 | 335,3 |
| ρ_s / ρ_l [kg/m³] | 2670,1 / 2513,8 (798,15 K); 2665,2 / 2506,9 (817,75 K) | 2587,4 / 2443,5; 2582,5 / 2436,9 | 2555,72 / 2378,33 | 2737,04 / 2527,05 |
| ΔS_V [J m⁻³K⁻¹] | −1,083×10⁶ | −8,878×10⁵ | −9,261×10⁵ | −9,799×10⁵ |
| γ₀ [J/m²] | 0,120 ± 0,003 | 0,106 ± 0,004 | 0,183 | 0,296 |

Arquivos: `coexist_<T>/`, `H_<fase>_<T>/`, `cfm_<orient>/` (trajetórias `prod.xyz`; o dump LAMMPS para o
`cfm_morris.py` se regenera com `python nep_tools.py todump cfm_<o>/prod.xyz cfm_<o>/cfm_nep89_<o>.dump`).
Decisão sobre o potencial: do autor.
