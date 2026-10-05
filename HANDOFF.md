# HANDOFF — Nucleation_MD (passagem para a máquina com GPU)

Data: 2026-10-04. Máquina de origem: notebook i7-1360P (12 núcleos / 16 threads, sem GPU NVIDIA),
Windows 11 + WSL2. Destino: máquina com **RTX 4050 Laptop (6 GB, Ada, sm_89)**.

Leitores: o autor (Prof. Ivaldo Leão Ferreira, UFPA) e a próxima sessão do Claude.
**A memória do Claude da máquina antiga não acompanha a pasta** — este arquivo a substitui.
Ler junto: `README.md`, `docs/CONVENTIONS.md`, `docs/CONTINUUM_VERIFICATION.md`, `md/e2_melting/E2_results.md`.

---

## 1. Objetivo do projeto

Construir a versão em dinâmica molecular (LAMMPS + PLUMED + Python) do modelo de nucleação de
Ferreira (2024), *Physica B*, doi:10.1016/j.physb.2024.416494, calculando atomisticamente as
grandezas do modelo e validando a solução do contínuo.

- ΔG(r) = [(1/3)πr³ΔS_VΔT + πr²γ_SL]·f(θ), f(θ) = 2 − 3cosθ + cos³θ; ΔS_V, ΔT, γ_SL, θ dependem de r.
- dΔG/dr = (πr/3)(a r² + 3b r + 6γf); r_c¹ = −2γf/b; r_c² = −3b/(2a).
- **Tensor de campo térmico Γ = A·∇T** (Gibbs–Thomson deslocado). Γ¹ˢᵗ = ΔT·r_c¹/2, Γ²ⁿᵈ = ΔT·r_c²/2.
- **Só o gradiente térmico desloca o equilíbrio** (palavra do autor). No experimento transiente ∇T = Ṫ/V_L.
- O SDAS experimental das duas ligas concorda com o modelo de 2ª ordem (MRB, Ferreira 2024).
- Entre escalas, **Γ vale nas duas** (palavra do autor); saltam ∇T, ΔT e r_c (r_c ~ nm na MD).
- Validação central da MD (confirmada pelo autor): seeds sob ∇T imposto por NEMD, r_c por seeding,
  verificar ΔT·r_c/2 = A·∇T = Γ.

### Preferências do autor (obrigatórias)
- **Não trocar as equações dele por formas da literatura**; onde houver ambiguidade, perguntar.
- **Perguntar antes de assumir qualquer valor físico** que ele não tenha fornecido.
- **As equações do contínuo são precisas** — não apresentar a construção do código dele como "erro";
  descrever fatos e perguntar.
- Trabalhar em etapas pequenas e testáveis, com testes unitários. Escreve em português.

---

## 2. Estado por etapa

| Etapa | Estado | Onde |
|---|---|---|
| E1 — pacote Python (CNT, Tolman, FerreiraModel, fontes, SDAS) | **concluída** (77 testes) | `nucleation_md/`, `tests/` |
| Verificação do contínuo (2 ligas) | **concluída** | `docs/CONTINUUM_VERIFICATION.md`, `scripts/verify_continuum.py` |
| Cadeia ∇T → Γ → SDAS | **concluída** (reproduz slides 25–30) | `nucleation_md/sdas.py`, `scripts/report_sdas.py` |
| CALPHAD (Al-Nb) | **concluída** | `scripts/calphad_liquidus.py`, `calphad_oc/`, `results/calphad/` |
| Instalação LAMMPS+PLUMED (CPU, WSL) | **concluída na máquina antiga** | `install/` |
| E2 — validação do potencial (Al) | **concluída para Borovikov 2024**; γ₀ abaixo da referência → teste de outros potenciais em aberto | `md/e2_melting/E2_results.md` |
| Tabelas contínuo × MD (Al-Si-Mg-Fe) | **rascunho pronto**; coluna MD parcial | `results/tables/`, `scripts/make_tables_docx.py` |
| E4+ — seeding, NEMD (Γ = A·∇T), PLUMED, FFS | **não iniciadas** (exigem GPU) | — |

---

## 3. Resultados principais

### 3.1 Contínuo (reproduzido pelo pacote, desvio ≤ 1e-5 = precisão do CSV)
- **Al-3Cu-5Nb-0,1Fe** (Mendes et al. 2023): ∇T_exp = 3294 K/m → Γ¹ˢᵗ = 3,887e-7, Γ²ⁿᵈ = 1,0685e-6 m·K,
  r_c² = 5,08 μm, ΔT = 0,4206 K. SDAS: MRB Γ¹ˢᵗ 7,67·t^0,377; Γ²ⁿᵈ 9,38·t^0,379; RB 27,59·t^0,333 μm.
- **Al-0,8Si-0,6Mg-0,2Fe** (Marques et al. 2025): ∇T_exp = 7902,38 K/m → Γ¹ˢᵗ = 7,637e-7, Γ²ⁿᵈ = 1,5745e-6 m·K,
  r_c² = 3,98 μm, ΔT = 0,7908 K. SDAS: 5,81·t^0,387; 6,71·t^0,387; RB 27,03·t^0,333.
- Métricas contra pontos experimentais (`results/sdas/report_sdas.md`): ambas as ordens MRB dentro da
  dispersão; Γ¹ˢᵗ mais próxima das médias, Γ²ⁿᵈ ~+20 %. Relatado ao autor sem alterar nada.

### 3.2 CALPHAD (OC6P e pycalphad, base ALFENB-B2-4SL/2SL, idênticas no canto rico em Al)
- Al-5 wt.% Nb: liquidus estável 1539 K (Al₃Nb primário); **peritético L + Al₃Nb → FCC_A1 a 934,51 K**;
  líquido nesse ponto com **0,0062 wt.% Nb**; FCC nascente com 0,516 wt.% Nb; Al₃Nb 53,4 wt.% Nb.
- Conclusão: o α-Al nasce de líquido praticamente sem Nb; o Nb está nas partículas de Al₃Nb (substrato).
- **Decisão do autor — opção (a):** líquido de Al, núcleo Al-0,5 wt.% Nb, substrato Al₃Nb (θ heterogêneo).
- **Ligas completas (2026-10-04, `scripts/calphad_alloys.py` → `results/calphad/alloys_TL_dH.md`, COST507-modified, sem a fase ALTI):**

  | liga | T_L(FCC_A1) | T_S (equil.) | ΔH_a = H_LIQ − H_FCC a T_L | ΔH_b = H(T_L) − H(T_S) | ΔH_c (só latente) | script |
  |---|---|---|---|---|---|---|
  | Al puro (controle) | 933,46 K | — | 397,0 kJ/kg | — | — | — |
  | Al-0,8Si-0,6Mg-0,2Fe | 925,14 K | 873,60 K | 385,2 | 452,3 | 390,3 | 925,32 K / 335,3 |
  | Al-3Cu-5Nb-0,1Fe | 925,36 K (Al₃Nb a 1598,7 K) | 856,50 K | 386,3 (x₀) / 381,1 (líquido a T_L) | 424,2 | 346,2 | 936,58 K / 335,3 |

  Fatos: o T_L do Al-Si-Mg-Fe bate com o script (Δ = 0,17 K). Nesta base os solutos baixam ΔH_a só ~3 % em relação ao Al puro.
  No Al-Cu-Nb-Fe o líquido que forma α-Al tem 3,31 Cu e 0,004 Nb (wt.%); a COST507 não dá o peritético de 934,5 K da ALFENB (que não tem Cu).
  **Autor (2026-10-05): referência = definição comparável à MD (ΔH_a); os 335,3 kJ/kg vêm de CALPHAD.** O manuscrito
  (`Full_Manuscript_Phase_Nucleation_edited_R2.docx`) cita Thermo-Calc + base **TTAl7** (comercial, não disponível aqui); a
  Tabela 1 da Luane registra "Latent heat (FCC_A1) 335300 J/kg". Outros scripts: 293 900 (Al-3Cu-0,5Mg), 260 300, 234 700.
  A COST507 não reproduz: ΔH_a(T) = 384–394 kJ/kg (ligas) e 397–406 (Al puro) entre 933 e 700 K; H(T_L) − H(T) cruza 335 kJ/kg
  só em f_s ≈ 0,79 (Al-Si-Mg-Fe) e 0,86 (Al-Cu-Nb-Fe), sem significado evidente. Falta saber qual grandeza do Thermo-Calc deu 335,3.

### 3.3 Potenciais testados
| Potencial | Veredito |
|---|---|
| **Borovikov, Mendelev et al. 2024, EAM/FS Ni-Al-Nb** (`potentials/NiAlNb_Borovikov2024/`) | **Escolhido.** Al₃Nb D0₂₂ estável (ΔH_f −0,314 eV/át; D0₂₂ < L1₂ por 2 meV), mas tetragonalidade não reproduzida (c/2a 1,02 × 1,12 exp) → substrato **congelado na rede experimental** (a = 3,841, c = 8,609 Å; ΔH_f −0,299 aí) |
| Fereidonnejad 2022, 2NN MEAM Al-Nb | Rejeitado: Al₃Nb com ΔH_f > 0 (+0,019 D0₂₂); o próprio NIST confirma |
| Jelinek 2012, MEAM AlSiMgCuFe | Rejeitado: T_m(Al) < 900 K; perde átomos na fusão a 1700 K |
| Farkas & Jones 1996, EAM Nb-Ti-Al | Rejeitado sem rodar: Al a₀ = 3,87 Å, energias de B2/Al₃Nb absurdas (NIST) |
| Mendelev 2008 Al (`potentials/Al_Mendelev2008/`) | Só Al; usado como comparação: T_m = 932,3 K |
| NEP Ti-Al-Nb (PRB 2024) | Indisponível (autores não publicaram) |
| **NEP89** (`potentials/NEP89_20250409/`, universal, 89 elementos) | Triagem estática: Al a = 3,99 Å; Al₃Nb ΔH_f −0,27 mas **L1₂ < D0₂₂** (errado). **E2 na GPU (2026-10-04, `md/e2_nep89/E2_nep89_results.md`): T_m = 919 ± 3 K, ΔH_m = 377,5 kJ/kg, ρ_s/ρ_l (920 K) = 2637/2473, γ₀ = 0,120 ± 0,003 J/m²** (Borovikov: 936, 325,5, 2552/2398, 0,106). Cobre Al, Nb, Cu, Fe, Si, Mg → permitiria ligas completas |

### 3.4 E2 — validação (Al, Borovikov 2024) — `md/e2_melting/E2_results.md`
| Grandeza | MD | Referência |
|---|---|---|
| T_m (coexistência isotérmica, 8192 át.) | **936 ± 3 K** | DTA 936,58; CALPHAD 934,51; Al 933,47 |
| ΔH_m (NPT a 935 K) | **325,5 kJ/kg** | contínuo 335,3 (−2,9 %); Al exp 397 |
| ρ_s / ρ_l a T_m | 2551,5 / 2397,9 kg/m³ | Al ≈ 2550 / 2375 |
| ρ_s / ρ_l a 798,15 K | 2587,4 / 2443,5 | Al-Si-Mg-Fe: 2555,72 / 2378,33 (+1,2 / +2,7 %) |
| ρ_s / ρ_l a 817,75 K | 2582,5 / 2436,9 | Al-Cu-Nb-Fe: 2737,04 / 2527,05 (liga com Cu e Al₃Nb) |
| ΔS_V | −8,878e5 J m⁻³K⁻¹ | contínuo −9,80e5 (densidades diferentes) |
| **γ₀** (flutuação capilar, 3 orientações) | **0,106 ± 0,004 J/m²** (só estatística; sistemático ~±20 %); ε₁ = 0,033, ε₂ ≈ 0 | Al exp ≈ 0,13–0,17; contínuo 0,183 (Al-Si-Mg-Fe), 0,296 (Al-Cu-Nb-Fe, r₀ = 6,3 μm) |

Tolerâncias aprovadas pelo autor: T_m ~1 %, ΔH_m ~10 %, ρ ~3 %. **γ₀ reprovado** → autor escolheu (c): testar outro potencial.

---

## 4. Decisão pendente (próximo passo imediato)

**2026-10-04: o autor escolheu a opção 1 (NEP89 na GPU) — feita; resultados em `md/e2_nep89/E2_nep89_results.md`.
Aguarda o veredito do autor sobre o potencial (NEP89 × Borovikov) antes da E4.**
**Feito (2026-10-05): opção 2 — γ₀ do Borovikov com fitas espessas** (`md/e2_cfm_thick/E2_thick_results.md`): W = 29–40 Å, 1 ns.
Reanálise comum das 9 fitas com piso de ruído (`fit_floor.py`, várias janelas): Borovikov fina 0,105–0,116, espessa 0,112–0,120,
NEP89 0,115–0,125 J/m². Espessura não muda γ₀ além de ~±10 %; NEP89 × Borovikov ~5 % (< sistemático). γ₀ não discrimina os potenciais.

Proposta enviada ao autor (aguarda escolha):
1. **Testar NEP89 de verdade na GPU** (T_m, ΔH_m, γ₀ com o mesmo protocolo) — GPUMD ou LAMMPS+NEP_CPU.
2. **Refazer γ₀ do Borovikov** com fitas mais espessas (W ≥ 2–3 nm) e 1–2 ns, para ver se 0,106 se sustenta.
3. Ajustar potencial próprio (fine-tuning do NEP89 com DFT de Al líquido, interface S/L e Al₃Nb) — projeto maior.
Recomendação dada: 1 e 2 em paralelo na GPU.

Outras perguntas em aberto ao autor:
- Tabela de propriedades da Luane (`Tabela 1 - Propriedades simulação Luane Final.docx`) × scripts:
  Fe 0,65 (tabela) × 0,2 wt.% (SDAS); título "Al-Cu-Si-Fe"; σ₀/γ₀ −0,914/0,154 (tabela) × 1,09/0,183
  (script de nucleação, que gera os Γ da tabela); V_L(P) 1,82P^−0,30 × 1,2P^−0,27; t_SL 0,93P^1,04 × 3,58P^1,09.
  Nas tabelas geradas usei 0,2 Fe e 1,09/0,183; V_L e t_SL ficaram de fora até confirmação.
- Origem dos 3,47 at.% Nb da FCC_A1 na L723 do script Al-Cu-Nb-Fe (CALPHAD dá ~0,15 at.% no peritético).
- k₀(Nb) = 2,049 no script de SDAS × ~2,8–4,9 no ramo metaestável do CALPHAD (só informativo).
- Tendência experimental de I (taxa): observável ainda não definido (teste `test_rate_trend_matches_experiment` pulado).
- Sinal do expoente da taxa: mantido literal (+ΔG*/ΔG*_Eq); o autor disse que a taxa é secundária.

Ideias discutidas para depois (aceitas como etapas futuras):
- Propriedades termofísicas vs T (ρ, c_p, η, D, γ_LV, γ_SV, C_ij, parâmetros de Gurtin–Murdoch σ₀, λ₀, μ₀).
  Condutividade térmica: MD só dá fônons; elétrons via Wiedemann–Franz (dados), Ziman (S(q) da MD, líquido)
  ou Kubo–Greenwood (DFT).
- Tensão × deformação vs T (cuidado: taxas 10⁷–10¹⁰ s⁻¹; usar policristal/defeitos).
- **Hipótese do autor: a entropia aumenta com o gradiente.** Teste: matriz G × V com o mesmo Ṫ = G·V
  (solidificação direcional NEMD), medindo entropia retida (`compute entropy/atom`, Frenkel–Ladd) e produzida
  (calor dos termostatos). Pergunta pendente: qual entropia é o critério (retida, produzida ou ambas).

---

## 5. Mapa de arquivos

| Caminho | Conteúdo |
|---|---|
| `nucleation_md/models.py` | `NucleationModel`, `CNT`, `Tolman`, `FerreiraModel` (coefficients, r_critical, exact_roots, gamma_tensor, grad_T, rate) |
| `nucleation_md/sources/continuum.py` | `ContinuumSource(alloy=…)` para Al3Cu5Nb01Fe (tablefull) e Al08Si06Mg02Fe (tabela_paper); `state_at_gradient` |
| `nucleation_md/sdas.py` | port literal do MRB/RB dos scripts de SDAS + dados das duas ligas |
| `nucleation_md/surface.py` | ΔG_S métrica/Ricci (literal) e Gurtin–Murdoch (do código) |
| `nucleation_md/validation.py` | critério do potencial (tolerâncias sem default) |
| `scripts/run_continuum.py` | roda scripts do autor em cópia temporária; extrai pontos de SDAS |
| `scripts/tablefull_alsimgfe.py` | cópia do script Al-Si-Mg-Fe que grava a tablefull (colunas comuns idênticas) |
| `scripts/make_tables_docx.py` | gera `results/tables/Tables_AlSiMgFe_continuum_vs_MD.docx` a partir de `results/md/properties_md.json` |
| `scripts/plot_scales.py` | figuras Γ×∇T, ΔT×r_c, SDAS (lê `results/md/<liga>/gamma_md.csv` quando existir) |
| `scripts/calphad_liquidus.py`, `calphad_oc/alnb_wt.OCM` | CALPHAD (bases em `C:\Users\ileao\OneDrive\Documentos\OpenCalphad\OC6`, fora do projeto) |
| `md/e2_melting/` | entradas e análises da E2: `in.coexist`, `in.enthalpy`, `in.cfm`, `analyze_tm.py`, **`cfm_morris.py`** (análise de γ₀ válida), `analyze_cfm.py` (versão CNA, ruidosa — usar só a parte `combine`) |
| `md/e2_substrate/` | testes do Al₃Nb (D0₂₂/L1₂) por potencial |
| `md/e2_nep89/screen_nep89.py` | triagem estática do NEP89 (calorine) |
| `md/gpu_benchmark/` | teste de velocidade na RTX 4050 (Windows OpenCL ou CUDA) |
| `install/` | builds LAMMPS 10Sep2025 + PLUMED 2.9.4 (CPU e CUDA) |
| `data/` | tabelas do contínuo (2ª liga) e pontos experimentais de SDAS |
| `results/` | E1, SDAS, figuras, CALPHAD, tabelas, `md/properties_md.json` |

Arquivos grandes: `md/e2_melting/*.dump` (~257 MB). **Manter** `cfm_borovikov_100/110/111.dump` (reanálise de γ₀);
os `coexist_*` e dumps antigos de alnb/jelinek podem ser apagados.

Scripts do autor (não modificar), **todos dentro do projeto**:
- `Continuum_Mechanics_Solution/Scripts/` — contínuo Al-Cu-Nb-Fe e Al-Si-Mg-Fe, SDAS, dados térmicos.
- `Continuum_Mechanics_Solution/Sim_Luane_Final/` — cópia byte a byte (2026-10-04, 61 arquivos, sem `__pycache__`)
  de `C:\MacbookPro\ufpa\PPGEM\SDAS Al_Si_Mg_Fe\Sim Luane Final`: tabela de propriedades da Luane
  (`Tabela 1 - Propriedades simulação Luane Final.docx`), scripts AlSiMgFe / AlCuMgFe (MR 2025) / H₂O,
  tabelas `tablefull_paper_*.csv`, subpastas `métrica/` e `scientific notation/`.
  `scripts/tablefull_alsimgfe.py` já usa esta cópia (resultado idêntico ao da pasta original).
- `Continuum_Mechanics_Solution/Paper_and_presentation/` — artigos (Physica B 2024, MR 2025, JMEP 2026),
  apresentação MRS 2026 e o manuscrito `Full_Manuscript_Phase_Nucleation_edited_R2.docx` (formato das Tabelas 1–2).

Bases do CALPHAD: `C:\Users\ileao\OneDrive\Documentos\OpenCalphad\OC6` — fora do projeto, mas **espelhadas pelo
OneDrive** (não precisam ser copiadas). O OC6P está instalado no WSL da distro `Ubuntu` da máquina antiga
(`~/OC6HOME/oc6P`); na máquina nova, reinstalar se for usar as macros de `calphad_oc/`.

---

## 6. Ambiente na máquina nova

### 6.0 Estado instalado na máquina GPU (2026-10-04) — projeto agora em `Documentos\Nucleation_MD_GPU`
i7-13650HX (20 threads), 16 GB (WSL vê ~7 GB), RTX 4050 Laptop 6 GB, driver 596.08. Instalado **sem sudo**:

| Item | Onde | Verificação |
|---|---|---|
| Python Windows 3.12.6 + pycalphad 0.11.2 | sistema | `python -m pytest`: 77 passam, 1 pulado |
| LAMMPS 10Sep2025 + PLUMED 2.9.4 + GPU (CUDA 12.6, mixed, sm_89) | WSL `Ubuntu-22.04`, `~/opt/lammps-10Sep2025-cuda`; cmake portátil em `~/opt/cmake` | smoke test CPU e PLUMED **bit-idênticos** a `log.wsl`/`log.wsl_plumed_1np`; `eam/fs/gpu` roda na RTX 4050; módulo Python `lammps` no venv |
| OC6P (OpenCalphad, git `cae3000`, 2026-09-28) | WSL `Ubuntu` (usuário `ileaohpc`), `~/OC6HOME/oc6P` | `calphad_oc/alnb_wt.OCM` reproduz as 9 temperaturas da máquina antiga (1539,27 / 934,51 / 1059,66 K …) |

| GPUMD (git `eb7fbbf`, 2026-10-03, sm_89) | WSL `Ubuntu-22.04`, `~/opt/src/GPUMD/src/{gpumd,nep}` (`install/build_gpumd_wsl.sh`); no PATH via `nucmd_cuda_env.sh` | NEP89, Al CFC 32 000 át., NPT 300 K: estável, a(300 K) = 4,01 Å, **5,6 Matom·passo/s** |
| calorine 4.0 + ase 3.29 | venv `~/opt/venvs/nucmd` | `md/e2_nep89/screen_nep89.py` reproduz os números da máquina antiga |

Uso: `MSYS_NO_PATHCONV=1 wsl -d Ubuntu-22.04 -- bash script.sh` (com `source ~/opt/nucmd_cuda_env.sh` dentro do script).
**Não passar laços com `$var` em `wsl -- bash -c '…'`**: o `wsl` reinterpreta a linha e expande as variáveis antes; usar script em arquivo.
Velocidade (EAM Borovikov, `md/gpu_benchmark/README.md`): GPU 8,7–11 Matom·passo/s (8,8 k–1 M át.); CPU 8 proc. 2,5–2,8 → GPU ≈ 3–4,6×.
Não instalado: build só-CPU (desnecessário: o build CUDA roda em CPU sem `-sf gpu`), LAMMPS Windows, NEP_CPU no LAMMPS.

1. Python: `pip install numpy scipy pandas matplotlib pytest pycalphad python-docx` → `python -m pytest` na raiz (77 passam, 1 pulado).
2. GPU rápido sem compilar: `md/gpu_benchmark/run_gpu_windows.bat` (LAMMPS Windows tem pacote GPU via OpenCL, precisão mista).
3. WSL: `bash install/build_lammps_plumed_cuda_wsl.sh` (pré-requisitos no cabeçalho: driver NVIDIA com WSL, `nvcc`, build-essential, cmake, openmpi — exigem sudo do autor). CPU: `install/build_lammps_plumed_wsl.sh`.
   Ativar: `source ~/opt/nucmd_cuda_env.sh` (ou `nucmd_env.sh`). Venv em `~/opt/venvs/nucmd` (sem ensurepip: script usa get-pip).
4. NEP89 na GPU: GPUMD (https://github.com/brucefan1983/GPUMD) ou LAMMPS com NEP_CPU; para cálculos estáticos, `pip install calorine` (no Windows falha por falta de compilador; no WSL funciona).
5. Caminhos fixos a revisar se o usuário/OneDrive mudar: `scripts/calphad_liquidus.py` e `calphad_oc/alnb_wt.OCM`
   (bases em `Documentos\OpenCalphad\OC6`, espelhadas pelo OneDrive) e os scripts `.sh` com
   `/mnt/c/Users/ileao/OneDrive/Documentos/Nucleation_MD_GPU/...` (já atualizados). Os scripts do autor já estão todos dentro do projeto.

---

## 7. Armadilhas já encontradas (não repetir)

**LAMMPS / protocolo**
- `fix npt` termostata a temperatura de **todos** os átomos (compute em `all`), não só do grupo do fix:
  com metade congelada o líquido fica 2× mais quente. Usar `fix nvt` no grupo, ou `fix_modify temp`.
- `fix npt … dilate <grupo>` (não `partial`).
- Coexistência em **NPH com x-y travados** deu T_m enviesado (~980 K); usar **enquadramento isotérmico NPT** (df_fcc/dt = 0).
  Excluir do ajuste temperaturas em que a fase some (saturação).
- `create_atoms` numera átomos conforme a decomposição MPI → mesma semente dá trajetórias diferentes com N processos diferentes. Gerar configurações uma vez (`write_data`) e reusar.
- **CNA a T_m subconta o sólido** (~20 % × ~⅔ real) → piso de ruído no espectro de flutuação. Usar o parâmetro de Morris (`cfm_morris.py`).
- **Fita (111) para flutuação capilar precisa de múltiplos de 3 células em y = [11-2]**; com 1 célula o cristal nasce defeituoso e funde. Conferir sempre `run 0` com fração FCC = 1 e E/átomo = coesiva.
- MEAM Jelinek perde átomos fundindo a 2500 K com volume do sólido; fundir a ~1700 K.
- Desempenho no i7 antigo: 8 processos MPI ótimo; MEAM ~2e5 e EAM ~2–8e5 átomo·passo/s (4 proc.).

**CALPHAD**
- COST507: a fase **ALTI** (L1₀ γ-TiAl) tem extremo rico em Al = CFC + 2 J/mol e L(Al,Nb) = −80800 + 30T → sem Ti ela substitui a FCC_A1 e dissolve o Al₃Nb. Excluir (`EXCLUDE` em `calphad_alloys.py`).

**Ambiente Windows/WSL**
- Duas distros (também na máquina GPU): `Ubuntu-22.04` (usuário `ileao`, LAMMPS) e `Ubuntu` (padrão, usuário `ileaohpc`, OC6P em `~/OC6HOME/oc6P`).
- Git Bash converte caminhos `/mnt/c/...` → usar `MSYS_NO_PATHCONV=1 wsl -d Ubuntu-22.04 -- …`.
- Corridas longas: lançar desacopladas (`setsid nohup bash script.sh … & sleep N`); o `sleep` é necessário senão o WSL encerra o processo. Esperas no terminal do Claude expiram em 30 min/2 h — as corridas continuam.
- `wsl --shutdown` resolveu um erro de serviço (E_FAIL) do WSL.
- pycalphad 0.11: `filter_phases(db, [v.Species(c) …])`; `COST507-modified.tdb` carrega, `COST507.tdb`/`.pycalphad`/`mc_al_v2036` não.
- OC6P em macro: linha vazia após a lista de elementos + outra linha vazia para o aviso da base; `set status phase X=fix 0` + `set cond t=none` dá a temperatura de aparecimento de X; `l r 4` lista em fração mássica. O autor chama o programa de **OC6P**.
- Word: conversão docx→pdf via COM do PowerShell funciona para conferência visual.

**Interpretação do contínuo (fatos verificados, não "erros")**
- ΔT tabelado = Gibbs–Thomson ΔT = 2Γ/r com Γ = −γ(r)/ΔS_V; ∇T = Γ/A.
- ∂ΔS_V/∂r (L902) é obtido da condição r_c² = r; no modo `closure` o pacote reproduz r_c² = r.
- ΔG* usa o ΔS_V da L222; a coluna DSv_hom (L937) entra em GB e nas figuras.
- Γ²ⁿᵈ (slide 10): o ∂ΔS_V/∂r do denominador é a derivada total a/(f·ΔT).

---

## 8. Próximos passos sugeridos (ordem)

1. Obter a escolha do autor no §4 (NEP89 / γ₀ robusto / potencial próprio).
2. Na GPU: benchmark (`md/gpu_benchmark`), build CUDA, e as corridas escolhidas (protocolos prontos em `md/e2_melting`).
3. Fechar o critério do potencial; atualizar `results/md/properties_md.json` e regenerar as tabelas.
4. E4: seeding homogêneo (seeds de 5–20 nm) e NEMD com `fix ehex` para Γ = A·∇T; gravar
   `results/md/<liga>/gamma_md.csv` (colunas `r_c_m, DT_K, gradT_K_per_m, Gamma_mK`) → `scripts/plot_scales.py`.
5. E5: substrato Al₃Nb congelado (rede experimental) para θ(r).
6. E7: PLUMED (Q6 → LOCAL_Q6 → MFILTER_MORE → CONTACT_MATRIX → DFSCLUSTERING → CLUSTER_NATOMS; testado no smoke test) e FFS.
