# Verificação do contínuo

Script: `Continuum_Mechanics_Solution/Scripts/Nucleation_model_FCC_A1_..._2nd_Sim_2026_paper.py` (linhas = Lxxx).
Nada foi alterado no script.

## 1. Recomputação independente da tabela

`scripts/verify_continuum.py` refaz o laço do script (L870–L1022) a partir das suas equações e
compara com `tablefull_paper_B_MRS_2026.csv`. O teste `tests/test_continuum_recomputation.py`
mantém essa verificação permanente.

| Coluna | Equação | Desvio relativo máx. |
|---|---|---|
| gam_hom | γ = γ₀(r₀/r)² − Σ_GM (L879) | 4,8e-6 |
| dfgamdr_ana | ∂γ/∂r (L883) | 4,5e-6 |
| DT | ΔT = 2Γ/r, Γ = −γ/ΔS_V (L890) | 3,4e-6 |
| DSv_hom | recorrência da L937 | 3,8e-6 |
| dDSv_homdr | L939 | 3,8e-6 |
| dTdr | ∇T = ΔT/(8πr) = Γ/A (L923) | 3,4e-6 |
| GT_het_2nd | Γ²ⁿᵈ = A·∇T = ΔT·r/2 (L930) | 9,4e-6 |
| GT_het | Γ¹ˢᵗ = −ΔT·γf/b (θ = π) | 8,8e-6 |
| GB, GS, GC | ΔG_V, ΔG_S, ΔG_C (L977–L979) | ≤ 4,8e-6 |
| DSs_hom, DSc | ΔS_S, ΔS_C (L973, L975) | ≤ 4,8e-6 |

Todos os desvios estão na precisão de impressão do CSV (`%g`, 6 algarismos).

## 2. O `FerreiraModel` reproduz o contínuo

Com `ContinuumSource(dSv_mode="closure", theta_mode="homogeneous")`:

- r_c² = −3b/(2a) = r em todas as linhas (rtol 1e-10);
- Γ¹ˢᵗ e Γ²ⁿᵈ (slide 10, MRS 2026) = colunas `GT_het` e `GT_het_2nd`, e Γ²ⁿᵈ/A = coluna `dTdr` (rtol 3e-5);
- no ponto experimental ∇T = 3294 K/m (Mendes et al. 2023): Γ¹ˢᵗ = 3,887e-7 e Γ²ⁿᵈ = 1,0685e-6 m·K, os valores do slide 25.

## 3. Como o contínuo constrói a solução (fatos, para orientar a MD)

- O gradiente é a variável que desloca o equilíbrio. Cada linha corresponde a um ∇T = Γ/A;
  r_c² e ΔT decorrem de Γ = A·∇T = ΔT·r_c²/2.
- O ΔT tabelado é Gibbs–Thomson com Γ = −γ(r)/ΔS_V. No `brentq(f_DT, 1e-4, 1)` (L890),
  esse valor é onde `f_DT` muda de sinal.
- ∂ΔS_V/∂r (L902) é obtido da condição r_c² = r. Na MD, ∂ΔS_V/∂r será medido entre seeds,
  o que dá um teste independente dessa condição.
- ΔG* (L914, L925) usa o ΔS_V da L222. A coluna `DSv_hom` (recorrência da L937) entra em GB e nas figuras 7 e 10.
  No pacote, o modo `dSv_mode="integrated"` usa essa coluna; a derivada da sua spline não é o ∂ΔS_V/∂r da L902,
  por isso o padrão é `closure`.
- `dgammadr_func_r` (L373) difere da derivada numérica de `gamma_func_r` em ~3e-6 relativo.
- A forma métrica/Ricci de ΔG_S está nos slides 13–16; no script, o γ(r) avaliado é o da L355 (Gurtin–Murdoch).

## 4. Salto entre escalas

O valor de Γ (Gibbs–Thomson) vale nas duas escalas. O salto ocorre em ∇T, ΔT e r_c, ligados por
Γ = A·∇T = ΔT·r_c/2. Com o Γ²ⁿᵈ de cada experimento:

| Liga (∇T exp.) | Γ²ⁿᵈ [m·K] | r_c = 5 nm | r_c = 10 nm | r_c = 20 nm |
|---|---|---|---|---|
| Al-3Cu-5Nb-0,1Fe (3294 K/m) | 1,0685e-6 | ΔT 427 K, ∇T 3,4e9 | ΔT 214 K, ∇T 8,5e8 | ΔT 107 K, ∇T 2,1e8 |
| Al-0,8Si-0,6Mg-0,2Fe (7902 K/m) | 1,57e-6 | ΔT 628 K, ∇T 5,0e9 | ΔT 314 K, ∇T 1,2e9 | ΔT 157 K, ∇T 3,1e8 |

ΔT em K, ∇T em K/m. A janela atomística compatível com o mesmo Γ é r_c ~ 10–20 nm, ΔT ~ 100–300 K,
∇T ~ 10⁸–10⁹ K/m: a faixa natural da NEMD (`fix ehex`). Um seed de 10 nm tem ~2,5×10⁵ átomos de Al.

## 5. Segunda liga e cadeia até o SDAS

- **Al-0,8Si-0,6Mg-0,2Fe** (Marques et al. 2025): `scripts/run_continuum.py` roda
  `Nucleation_model_..._AlSiMgFe_..._Luane.py` numa cópia temporária e guarda `data/continuum/Al08Si06Mg02Fe/tabela_paper.csv`.
  O cabeçalho desse CSV tem 10 nomes e cada linha tem 11 valores (`GT_het` não está no cabeçalho, L762 × L875);
  a leitura segue o formato das linhas.
  Em todas as 77 linhas: r_c² = r, Γ¹ˢᵗ = `GT_het`, Γ²ⁿᵈ = `GT_het_2nd`, Γ²ⁿᵈ/A = `dTdr`.
  Em ∇T = 7902,38 K/m: Γ¹ˢᵗ = 7,637e-7 e Γ²ⁿᵈ = 1,5745e-6 m·K (slide 29).
- **SDAS** (`nucleation_md/sdas.py`): port literal de `diff_length_scale`, `calc` (MRB) e `calc_RB`.
  Reproduz o script do Al-Si-Mg-Fe com diferença 0 e os ajustes dos slides 27 e 30.
  Dois detalhes reproduzidos como estão no script, com efeito desprezível ou restrito:
  o laço termina quando algum ponto de t_SL converge (diferença < 1e-8 frente à convergência de todos),
  e β₃ usa Fo₂ no numerador (`legacy_beta3=True`).
- O script de SDAS do Al-3Cu-5Nb-0,1Fe não executa (L342); seus dados foram copiados para `ALLOY_SDAS`
  e os pontos experimentais extraídos do texto para `data/experimental/sdas_Al3Cu5Nb01Fe.csv`.
- Relatório ponta a ponta: `scripts/report_sdas.py` → `results/sdas/`.
