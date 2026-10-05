# E2 — validação do potencial (Al puro)

Potencial: Borovikov, Mendelev, Zarkevich, Smith, Lawson, Int. J. Plast. 178 (2024) 104004, EAM/FS Ni-Al-Nb.

## T_m — coexistência isotérmica (NPT, 8192 átomos, 100 ps, fração FCC por CNA)

| T [K] | df_fcc/dt [1/ns] |
|---|---|
| 900 | +0,958 |
| 930 | +0,265 |
| 960 | −0,761 |

T_m = 936 ± 3 K (ajuste linear 900–960 K: 935,4 K; interpolação 930–960 K: 937,7 K; a corrida de 990 K fundiu por inteiro e satura). Comparação: Mendelev 2008 Al, 932,3 K. Referências: DTA 936,58 K; CALPHAD (peritético, OC6P) 934,51 K; Al puro 933,47 K.

## ΔH_m e densidades — NPT a 935 K, 2048 átomos, 100 ps

| | MD | referência do contínuo | Al puro exp. |
|---|---|---|---|
| ΔH_m | 325,5 kJ/kg (8,782 kJ/mol) | 335,3 kJ/kg | 397 kJ/kg |
| ρ_s | 2551,5 kg/m³ | 2737,04 | ≈ 2550 |
| ρ_l | 2397,9 kg/m³ | 2527,05 | ≈ 2375 |
| ΔV/V_s | 6,41 % | — | ≈ 6,5–7 % |
| ΔS_V = −ΔH·ρ_s/T_m | −8,878×10⁵ J m⁻³ K⁻¹ | −9,799×10⁵ | — |

## Massas específicas nas temperaturas eutéticas dos scripts de SDAS (NPT, 2048 átomos)

| Liga (T_s) | ρ_s MD | ρ_s script | ρ_l MD | ρ_l script |
|---|---|---|---|---|
| Al-0,8Si-0,6Mg-0,2Fe (798,15 K) | 2587,4 | 2555,72 (+1,2 %) | 2443,5 | 2378,33 (+2,7 %) |
| Al-3Cu-5Nb-0,1Fe (817,75 K) | 2582,5 | 2737,04 (−5,7 %) | 2436,9 | 2527,05 (−3,6 %) |

MD = Al puro; os valores do Al-Cu-Nb-Fe incluem Cu e Al₃Nb (mais densos). O líquido super-resfriado não cristalizou (ΔV < 0,1 %).

## γ₀ — método de flutuação capilar

Fitas a T_m (NPH, 400 ps; T média 924–936 K), alturas da interface pelo parâmetro de ordem de Morris
(cfm_morris.py; o CNA a T_m subconta o sólido e gera piso de ruído). Ajuste na janela capilar, modos 4–10.

| Orientação | W [Å] | inclinação log-log (ideal −2) | γ̃ [J/m²] |
|---|---|---|---|
| (100)[100] | 12,4 | −2,15 | 0,094 ± 0,004 |
| (110)[001] | 11,7 | −2,68 | 0,096 ± 0,008 |
| (111)[1-10] | 20,2 | −1,64 | 0,116 ± 0,005 |

Expansão cúbica (Hoyt, Asta, Karma 2001), incerteza por Monte Carlo sobre os erros estatísticos:

**γ₀ = 0,106 ± 0,004 J/m²**, ε₁ = 0,033 ± 0,010, ε₂ = −0,001 ± 0,002
→ γ(100) = 0,107, γ(110) = 0,105, γ(111) = 0,105 J/m².

Incertezas sistemáticas não incluídas: fitas finas (1,2–2 nm), 400 ps (modos longos não relaxados),
inclinações fora de −2 em (110) e (111). Ordem de grandeza provável: ±20 %.
Referências: Al experimental ≈ 0,13–0,17 J/m²; contínuo γ₀ = 0,183 (Al-Si-Mg-Fe) e 0,296 (Al-Cu-Nb-Fe, r₀ = 6,3 μm).
1ª tentativa da (111) inválida (fita com 1 célula em [11-2] não é periódica; refeita com 3).

Pendente: tolerâncias do critério (decisão do autor).
