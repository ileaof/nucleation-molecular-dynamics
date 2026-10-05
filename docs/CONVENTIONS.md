# Convenções do pacote `nucleation_md`

Unidades SI em todo o código. As equações seguem Ferreira (2024), na forma enviada pelo autor;
onde ainda há ambiguidade, a escolha é um **parâmetro** com o padrão literal indicado abaixo.

## Sinais

| Grandeza | Sinal | Origem |
|---|---|---|
| ΔS_V [J m⁻³ K⁻¹] | < 0 | solidificação; −ΔH·ρ_s/T_L no contínuo (L222) |
| ΔT [K] | > 0 | sub-resfriamento, T = T_m − ΔT |
| ΔG_V = ΔS_V·ΔT | < 0 | força motriz |
| ΔG_S = γ_SL [J m⁻²] | > 0 | |
| ΔG* = ΔG(r_c) | > 0 nos dados do contínuo | resultado, não imposto |
| ΔG*_Eq | > 0 nos dados atuais | ΔG(r_c) no estado de referência |

Testes: `tests/test_sign_convention.py`.

## Nomes (os dois aparecem como ΔG_C no artigo)

| Código | Significado |
|---|---|
| `dG_conf` | energia livre configuracional no denominador de r_c: γ_SL·(∂f/∂r)/f |
| `dG_star` / `barrier()` | energia livre crítica ΔG(r_c), usada na taxa |
| `dG_eq` | ΔG*_Eq = `barrier(reference)` |

## Equações implementadas (`nucleation_md/models.py`)

- ΔG(r) = [(1/3)πr³ΔS_VΔT + πr²γ_SL]·f(θ), com f(θ) = 2 − 3cosθ + cos³θ (f(π) = 4)
- dΔG/dr = (πr/3)(a r² + 3b r + 6γf), verificado numericamente (`test_dG_dr_equals_parabola`)
- a, b, r_c¹ = −2γf/b, r_c² = −3b/(2a) exatamente como definidos; as raízes exatas são só reportadas
- b/f = ΔG_V + ∂ΔG_S/∂r + ΔG_C
- Γ¹ˢᵗ = −ΔG_S/(ΔS_V + ∂ΔS_S/∂r + ΔS_C) = ΔT·r_c¹/2, com entropias S ≡ G/ΔT (slide 10, MRS 2026)
- Γ²ⁿᵈ = −(3/4)(ΔG_V + ∂ΔG_S/∂r + ΔG_C)/(∂ΔS_V/∂r) = ΔT·r_c²/2, onde ∂ΔS_V/∂r é a derivada total
  (1/ΔT)·∂(ΔS_VΔT·f)/∂r / f = a/(f·ΔT) (slide 10)
- **Γ = A·∇T**: só o gradiente térmico desloca o equilíbrio; ∇T = Γ/A, A = 2πr_c²(1 − cosθ).
  No experimento transiente, ∇T = G = Ṫ/V_L.
- Γ¹ˢᵗ, Γ²ⁿᵈ e ∇T reproduzem as colunas GT_het, GT_het_2nd e dTdr do contínuo (rtol 3e-5)
- I = (D·A/λ⁴)·(N/V)·exp(s·ΔG*/ΔG*_Eq), `exponent_sign` s = +1 (literal)

## Escolhas configuráveis, pendentes de confirmação

| Opção | Padrão (literal) | Alternativa | Onde |
|---|---|---|---|
| ∂f/∂r | `df_mode="chain"`: f′(θ)·∂θ/∂r = 3sin³θ·∂θ/∂r | `"continuum_legacy"`: `dfthetadr(θ)` da L391 | `State` |
| sinal do expoente da taxa | `exponent_sign=+1` | −1 | `FerreiraModel` |
| r_c usado em barrier/rate | `order=2` | 1 | `FerreiraModel` |
| r_c explícito ou autoconsistente | `solve="explicit"` (coeficientes em `r_eval`) | `"self_consistent"`, r = r_c(r) | `FerreiraModel.r_critical` |
| forma de ΔG_S | ambas disponíveis | | `surface.metric_2sphere`, `surface.gurtin_murdoch_legacy` |
| θ no termo métrico | θ(r) do núcleo (ângulo de molhamento) | | `metric_2sphere(theta=...)` |
| limites de ∫Σ dA | A0 = A(r₀), A1 = A(r), A = 4πr² | `area=` | `metric_2sphere` |
| ΔS_V do contínuo | `dSv_mode="closure"` (L222 + L902), a construção do contínuo | `"integrated"` (coluna DSv_hom), `"constant"` | `ContinuumSource` |
| estado pelo gradiente | `ContinuumSource.state_at_gradient(∇T)` (interpolação da L1036) | `state(i)` | `ContinuumSource` |

## Teste de tendência com dados experimentais

`tests/test_rate_trend_matches_experiment` fica pulado até que
(1) os dados (Mendes et al. 2023; Marques et al. 2025) estejam em `data/experimental/` e
(2) esteja definido qual observável experimental representa a tendência de I (ex.: densidade de grãos, 1/SDAS³, …).

## Critério para o potencial interatômico

`nucleation_md.validation.potential_gate(measured, reference, rtol)` bloqueia qualquer comparação se
T_m, ΔH_m ou γ₀ da MD saírem da tolerância. As tolerâncias não têm valor padrão.
