# Nucleation MD — modelo de nucleação de Ferreira (2024): contínuo × dinâmica molecular

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)
![LAMMPS](https://img.shields.io/badge/LAMMPS-10Sep2025-orange)
![PLUMED](https://img.shields.io/badge/PLUMED-2.9.4-orange)
![GPUMD](https://img.shields.io/badge/GPUMD-NEP89-orange)

> *English summary:* atomistic (LAMMPS + PLUMED + GPUMD) and Python evaluation of the solid–liquid nucleation
> model of Ferreira (2024, *Physica B*, [doi:10.1016/j.physb.2024.416494](https://doi.org/10.1016/j.physb.2024.416494)),
> in which the gradient of the thermal field **T** displaces the Gibbs–Thomson equilibrium through the thermal-field
> tensor **Γ** = **A**·∇**T**. The goal is to compute the model quantities by molecular dynamics and validate the continuum solution
> for two Al alloys (Al-3Cu-5Nb-0.1Fe and Al-0.8Si-0.6Mg-0.2Fe). Documentation is in Portuguese.

Este repositório reúne o **pacote Python** que reproduz a solução do contínuo, as **entradas e análises de dinâmica
molecular** (MD) que calculam as mesmas grandezas na escala atômica, e os cálculos **CALPHAD** que dão as referências
termodinâmicas das ligas.

---

## Sumário

- [O modelo](#o-modelo)
- [Estado do projeto](#estado-do-projeto)
- [Resultados principais](#resultados-principais)
- [Início rápido (Python)](#início-rápido-python)
- [Reproduzir os resultados](#reproduzir-os-resultados)
- [Simulações de MD](#simulações-de-md)
- [Estrutura do repositório](#estrutura-do-repositório)
- [Documentação](#documentação)
- [Como citar](#como-citar)
- [Licença e conteúdo de terceiros](#licença-e-conteúdo-de-terceiros)

---

## O modelo

Formulação de Ferreira (2024) na forma da apresentação no MRS Meeting 2026
([`Continuum_Mechanics_Solution/Paper_and_presentation/`](Continuum_Mechanics_Solution/Paper_and_presentation), slides 4–11).

**Notação.** $\mathbf{T}$ (negrito) é o **campo térmico**, em K, e $\nabla\mathbf{T}$ o seu gradiente. $T$ (sem negrito) é
a temperatura e $\Delta T$ o sub-resfriamento. $\boldsymbol{\Gamma}$ é o tensor de campo térmico e $\mathbf{A}$ a área da
nova superfície, normal ao campo $\mathbf{T}$.

### 1. Campo térmico a partir da 1ª lei

Pela 1ª lei, a energia $E$ é função de ponto, e $Q$ e $W$ dependem do caminho:

```math
\delta Q - \delta W = dE
```

Tomando $dE \equiv \partial E_{,i}$ em relação à variável primitiva $T$, a 1ª lei escreve-se como
$\nabla\left[\frac{\partial T}{\delta Q}(Q-W)\right] = \nabla\left[\frac{\partial T}{\delta Q}E\right]$, o que define dois
gradientes de campo térmico:

```math
\nabla\mathbf{T}_{\mathrm{pi}} = \nabla\left[\frac{\partial T}{\delta Q}\,(Q-W)\right] = \nabla\left[\frac{\partial T}{\delta Q}\,E\right]
\qquad\text{(independente do caminho)}
```

```math
\nabla\mathbf{T}_{\mathrm{pd}} = \nabla\left[\frac{\partial T}{\delta Q}\,Q\right]
= \nabla\left[\frac{\partial T}{\delta Q}\,E\right] + \nabla\left[\frac{\partial T}{\delta Q}\,W\right]
= \nabla\mathbf{T}_{\mathrm{pi}} + \nabla\left[\frac{\partial T}{\delta Q}\,W\right]
\qquad\text{(dependente do caminho)}
```

Como $E$ e $W$ são variáveis independentes, os termos cruzados $\nabla\left[\frac{\partial T}{\delta W}E\right]$ e
$\nabla\left[\frac{\partial T}{\partial E}W\right]$ podem ser desprezados, e o campo de Ferreira (2024) fica

```math
\nabla\mathbf{T} = \nabla\left[\frac{\partial T}{\delta Q}\,(E+W)\right]
\cong \nabla\left[\frac{\partial T}{\partial E}\,E\right] + \nabla\left[\frac{\partial T}{\delta W}\,W\right]
```

### 2. Energia aberta em suas formas

Com $E = U + K + E_P + \Sigma_{other}$ (simetria de Emmy Noether) e
$\delta Q - \delta W = dE = dU + dK + dE_P + d\Sigma_{other}$, sendo $U = f(T, V, P, n_1 \ldots n_N)$,

```math
dU = \frac{\partial U}{\partial T}\,dT + \frac{\partial U}{\partial V}\,dV + \frac{\partial U}{\partial P}\,dP
+ \sum_{i}^{N}\frac{\partial T}{\partial n_i}\,dn_i
```

```math
\nabla\mathbf{T} \cong \nabla\left[\frac{\partial T}{\partial U + \partial K + \partial E_P + \partial\Sigma_{other}}
\left(U + K + E_P + \Sigma_{other}\right)\right] + \nabla\left[\frac{\partial T}{\delta W}\,W\right]
```

e, abrindo a energia interna nas formas translacional, vibracional, rotacional, eletrônica, magnética e outras,
mais as energias cinética e potencial, outras formas de energia e o trabalho:

```math
\nabla\mathbf{T} \cong
\nabla\left[\frac{\partial T}{\partial U}\left(U_{transl} + U_{vib} + U_{rot} + U_{elect} + U_{mag} + U_{other}\right)\right]
+ \nabla\left[\frac{\partial T}{\partial K}\,K\right]
+ \nabla\left[\frac{\partial T}{\partial E_P}\,E_P\right]
+ \nabla\left[\frac{\partial T}{\partial \Sigma_{other}}\,\Sigma_{other}\right]
+ \nabla\left[\frac{\partial T}{\delta W}\,W\right]
```

### 3. Tensor de campo térmico

Toda nova superfície criada ou deformada, de área $\mathbf{A}$, é normal ao campo térmico $\mathbf{T}$:

```math
\boldsymbol{\Gamma} = \mathbf{A}\cdot\nabla\mathbf{T} = \mathbf{A}\cdot\nabla\mathbf{T}_{\mathrm{pi}}
= \mathbf{A}\cdot\nabla\left[\frac{\partial T}{\delta Q}\,(Q-W)\right] \cong \mathbf{A}\cdot\nabla\left[\frac{\partial T}{\partial E}\,E\right]
```

### 4. Nucleação

Energia livre total do núcleo, em termos das energias livres de volume $\Delta G_V = \Delta S_V\,\Delta T$ e de superfície
$\Delta G_S = \gamma_{SL}$, com $\Delta S_V$, $\Delta T$, $\gamma_{SL}$ e $\theta$ funções do raio:

```math
\Delta G = \left(\frac{1}{3}\pi r_C^3\,\Delta S_V\,\Delta T + \pi r_C^2\,\gamma_{SL}\right) f(\theta)
= \left(\frac{1}{3}\pi r_C^3\,\Delta G_V + \pi r_C^2\,\Delta G_S\right) f(\theta),
\qquad f(\theta) = 2 - 3\cos\theta + \cos^3\theta
```

O máximo em relação a $r_C$ dá a equação parabólica

```math
\left[\left(\frac{\partial \Delta S_V}{\partial r}\Delta T + \Delta S_V\frac{\partial \Delta T}{\partial r}\right) f(\theta)
+ \Delta S_V\,\Delta T\,\frac{\partial f(\theta)}{\partial r}\right] r_C^2
+ 3\left[\left(\Delta S_V\,\Delta T + \frac{\partial \gamma_{SL}}{\partial r}\right) f(\theta)
+ \gamma_{SL}\frac{\partial f(\theta)}{\partial r}\right] r_C + 6\gamma_{SL}\, f(\theta) = 0
```

**Raio crítico de 2ª ordem**, em termos das energias livres de volume $\Delta G_V$, de superfície $\Delta G_S$ e
configuracional $\Delta G_C$:

```math
r_C^{2nd} = -\frac{3}{2}\,
\frac{\left(\Delta S_V\,\Delta T + \dfrac{\partial \gamma_{SL}}{\partial r}\right) + \dfrac{\gamma_{SL}}{f(\theta)}\dfrac{\partial f(\theta)}{\partial r}}
{\left(\dfrac{\partial \Delta S_V}{\partial r}\Delta T + \Delta S_V\dfrac{\partial \Delta T}{\partial r}\right) + \dfrac{\Delta S_V\,\Delta T}{f(\theta)}\dfrac{\partial f(\theta)}{\partial r}}
= -\frac{3}{2}\,
\frac{\left(\Delta S_V\,\Delta T + \dfrac{\partial \gamma_{SL}}{\partial r}\right) f(\theta) + \gamma_{SL}\dfrac{\partial f(\theta)}{\partial r}}
{\dfrac{\partial\left(\Delta S_V\,\Delta T\, f(\theta)\right)}{\partial r}}
= -\frac{3}{2}\,\frac{\Delta G_V + \dfrac{\partial \Delta G_S}{\partial r} + \Delta G_C}{\dfrac{\partial \Delta G_V}{\partial r}}
```

**Raio crítico de 1ª ordem** (aproximação) e comparação com a teoria clássica (CNT):

```math
r_C^{1st} = -\frac{2\gamma_{SL}\, f(\theta)}
{\left(\Delta S_V\,\Delta T + \dfrac{\partial \gamma_{SL}}{\partial r}\right) f(\theta) + \gamma_{SL}\dfrac{\partial f(\theta)}{\partial r}}
= -\frac{2\,\Delta G_S}{\Delta G_V + \dfrac{\partial \Delta G_S}{\partial r} + \Delta G_C},
\qquad r_C^{CNT} = -\frac{2\,\Delta G_S}{\Delta G_V}
```

**Tensor de campo térmico de 2ª ordem**, do máximo da energia livre de Gibbs total:

```math
\boldsymbol{\Gamma}^{2nd} = -\frac{3}{4}\,
\frac{\left[\left(\Delta S_V\,\Delta T + \dfrac{\partial \gamma_{SL}}{\partial r}\right) + \dfrac{\gamma_{SL}}{f(\theta)}\dfrac{\partial f(\theta)}{\partial r}\right]\Delta T}
{\left(\dfrac{\partial \Delta S_V}{\partial r}\Delta T + \Delta S_V\dfrac{\partial \Delta T}{\partial r}\right) + \dfrac{\Delta S_V\,\Delta T}{f(\theta)}\dfrac{\partial f(\theta)}{\partial r}}
= -\frac{3}{4}\,
\frac{\left(\Delta S_V\,\Delta T + \dfrac{\partial \gamma_{SL}}{\partial r}\right) f(\theta) + \gamma_{SL}\dfrac{\partial f(\theta)}{\partial r}}
{\dfrac{1}{\Delta T}\dfrac{\partial\left(\Delta S_V\,\Delta T\, f(\theta)\right)}{\partial r}}
= -\frac{3}{4}\,\frac{\Delta G_V + \dfrac{\partial \Delta G_S}{\partial r} + \Delta G_C}{\dfrac{\partial \Delta S_V}{\partial r}}
= \mathbf{A}\cdot\nabla\mathbf{T}
```

**Tensor de campo térmico de 1ª ordem** (aproximação), com entropias $S \equiv G/\Delta T$:

```math
\boldsymbol{\Gamma}^{1st} = -\frac{\gamma_{SL}\, f(\theta)\,\Delta T}
{\left(\Delta S_V\,\Delta T + \dfrac{\partial \gamma_{SL}}{\partial r}\right) f(\theta) + \gamma_{SL}\dfrac{\partial f(\theta)}{\partial r}}
= -\frac{\Delta G_S\,\Delta T}{\Delta G_V + \dfrac{\partial \Delta G_S}{\partial r} + \Delta G_C}
= -\frac{\Delta G_S}{\Delta S_V + \dfrac{\partial \Delta S_S}{\partial r} + \Delta S_C}
```

Equivalentemente, $\boldsymbol{\Gamma}^{1st} = \tfrac{1}{2}\,\Delta T\, r_C^{1st}$ e
$\boldsymbol{\Gamma}^{2nd} = \tfrac{1}{2}\,\Delta T\, r_C^{2nd}$ (forma de Gibbs–Thomson).

**Taxa de nucleação:**

```math
I = \frac{D\,A}{\lambda^4}\,\frac{N}{V}\,\exp\left(\frac{\Delta G_C}{\Delta G_{C,Eq}}\right)
```

em que $D$ é o coeficiente de autodifusão (m²/s), $N/V$ a densidade de modos da nova fase (modos/m³), $\lambda$ a distância
que um átomo salta até a nova superfície (m) e $A$ a área da nova fase (m²).

### 5. Ligação com o experimento e com a MD

Só o gradiente do campo térmico desloca o equilíbrio. No experimento de solidificação transiente,
$\nabla\mathbf{T}$ é obtido como $G = \dot T / V_L$. A partir de $\boldsymbol{\Gamma}$, o modelo MRB prevê o espaçamento
dendrítico secundário (SDAS), que concorda com os dados experimentais das duas ligas.

**Hipótese testada pela MD:** $\boldsymbol{\Gamma}$ vale nas duas escalas. Entre o contínuo ($r_C$ ~ μm,
$\nabla\mathbf{T}$ ~ 10³ K/m) e a MD ($r_C$ ~ nm, $\nabla\mathbf{T}$ ~ 10⁸–10¹⁰ K/m) mudam $\nabla\mathbf{T}$, $\Delta T$ e
$r_C$, mas $\boldsymbol{\Gamma} = \tfrac{1}{2}\Delta T\, r_C = \mathbf{A}\cdot\nabla\mathbf{T}$ deve se manter. O teste são
núcleos (*seeds*) sob $\nabla\mathbf{T}$ imposto por NEMD, com $r_C$ obtido por *seeding*.

<p align="center">
  <img src="results/figures/fig_gamma_gradT.png" width="820"
       alt="Γ × ∇T e ΔT × r_c para as duas ligas: solução do contínuo e previsão atomística com Γ invariante entre escalas">
</p>

*Γ × ∇**T** (a, c) e ΔT × r_c (b, d): solução do contínuo e previsão na escala atômica com o mesmo Γ. Os marcadores vazados
serão substituídos pelos resultados de MD (`results/md/<liga>/gamma_md.csv`).*

## Estado do projeto

| Etapa | Estado | Onde |
|---|---|---|
| **E1** — pacote Python (CNT, Tolman, `FerreiraModel`, fontes do contínuo, SDAS) | ✅ concluída, 77 testes | [`nucleation_md/`](nucleation_md), [`tests/`](tests) |
| Verificação do contínuo (2 ligas) | ✅ desvio ≤ 1×10⁻⁵ | [`docs/CONTINUUM_VERIFICATION.md`](docs/CONTINUUM_VERIFICATION.md) |
| Cadeia ∇**T** → Γ → SDAS | ✅ reproduz os scripts originais | [`nucleation_md/sdas.py`](nucleation_md/sdas.py), [`results/sdas/`](results/sdas) |
| CALPHAD (Al-Nb e ligas completas) | ✅ | [`scripts/calphad_*.py`](scripts), [`results/calphad/`](results/calphad) |
| Instalação LAMMPS + PLUMED (CPU e CUDA) e GPUMD | ✅ | [`install/`](install) |
| **E2** — validação do potencial interatômico (Al) | ✅ Borovikov 2024 e NEP89; veredito pendente | [`md/e2_melting/`](md/e2_melting), [`md/e2_nep89/`](md/e2_nep89), [`md/e2_cfm_thick/`](md/e2_cfm_thick) |
| **E4** — *seeding* e NEMD (Γ = A·∇**T**) | ⏳ a iniciar | — |
| **E5** — substrato Al₃Nb, θ(r) | ⏳ | — |
| **E7** — PLUMED (Q6, clusters) e *forward flux sampling* | ⏳ (protocolo testado no *smoke test*) | [`install/smoke_test/`](install/smoke_test) |

## Resultados principais

### Contínuo (reproduzido pelo pacote)

| Liga | ∇**T** experimental | Γ¹ˢᵗ [m·K] | Γ²ⁿᵈ [m·K] | r_c² | ΔT |
|---|---|---|---|---|---|
| Al-3Cu-5Nb-0,1Fe (Mendes et al. 2023) | 3294 K/m | 3,887×10⁻⁷ | 1,0685×10⁻⁶ | 5,08 μm | 0,4206 K |
| Al-0,8Si-0,6Mg-0,2Fe (Marques et al. 2025) | 7902,38 K/m | 7,637×10⁻⁷ | 1,5745×10⁻⁶ | 3,98 μm | 0,7908 K |

<p align="center">
  <img src="results/sdas/sdas_Al3Cu5Nb01Fe.png" width="430" alt="SDAS × tempo local de solidificação da liga Al-3Cu-5Nb-0,1Fe: modelos MRB e RB contra Mendes et al. (2023)">
  <img src="results/sdas/sdas_Al08Si06Mg02Fe.png" width="430" alt="SDAS × tempo local de solidificação da liga Al-0,8Si-0,6Mg-0,2Fe: modelos MRB e RB contra Marques et al. (2025)">
</p>

*SDAS × tempo local de solidificação t_SL. À esquerda, Al-3Cu-5Nb-0,1Fe (Mendes et al. 2023); à direita,
Al-0,8Si-0,6Mg-0,2Fe (Marques et al. 2025). Nas duas ligas o modelo MRB, com Γ¹ˢᵗ ou Γ²ⁿᵈ, fica dentro da dispersão
experimental, e o RB superestima o SDAS. Ajustes MRB: 7,67·t^0,377 e 9,38·t^0,379 μm (Al-Cu-Nb-Fe); 5,81·t^0,387 e
6,71·t^0,387 μm (Al-Si-Mg-Fe). Métricas em [`results/sdas/report_sdas.md`](results/sdas/report_sdas.md).*

### CALPHAD

- **Al-5 wt.% Nb** (base ALFENB): Al₃Nb primário; peritético L + Al₃Nb → CFC a **934,51 K**, com líquido de apenas
  0,0062 wt.% Nb. O α-Al nasce de líquido praticamente sem Nb, e o Al₃Nb atua como substrato.
- **Ligas completas** (COST507): T_L(CFC) = 925,14 K (Al-Si-Mg-Fe) e 925,36 K (Al-Cu-Nb-Fe). A entalpia de fusão
  H_LIQ − H_CFC na composição da liga é 385 kJ/kg (397 kJ/kg para Al puro).
  Detalhes em [`results/calphad/alloys_TL_dH.md`](results/calphad/alloys_TL_dH.md).

### Validação do potencial (E2, Al puro)

| Grandeza | Borovikov et al. 2024 (EAM/FS) | NEP89 (GPUMD) | Referência |
|---|---|---|---|
| T_m (coexistência) | 936 ± 3 K | 919 ± 3 K | DTA 936,58 K; Al 933,47 K |
| ΔH_m | 325,5 kJ/kg | 377,5 kJ/kg | Al exp. 397 kJ/kg |
| ρ_s / ρ_l a T_m | 2552 / 2398 kg/m³ | 2637 / 2473 kg/m³ | Al ≈ 2550 / 2375 kg/m³ |
| γ₀ (flutuação capilar)\* | 0,11–0,12 J/m² | 0,115–0,125 J/m² | Al exp. ≈ 0,13–0,17 J/m² |
| Al₃Nb D0₂₂ | estável ✔ | L1₂ < D0₂₂ ✘ | D0₂₂ |

\* Faixa entre janelas de ajuste, com piso de ruído, sobre fitas finas e espessas ([`md/e2_cfm_thick/E2_thick_results.md`](md/e2_cfm_thick/E2_thick_results.md)).
Potenciais rejeitados e motivos em [`HANDOFF.md`](HANDOFF.md) §3.3.

## Início rápido (Python)

```bash
git clone https://github.com/ileaof/nucleation-molecular-dynamics.git
cd nucleation-molecular-dynamics
pip install numpy scipy pandas matplotlib pytest      # + pycalphad python-docx para CALPHAD e tabelas
python -m pytest                                      # 77 passam, 1 pulado
```

Raios críticos e tensor Γ da liga Al-3Cu-5Nb-0,1Fe no gradiente experimental:

```python
from nucleation_md import FerreiraModel
from nucleation_md.sources.continuum import ContinuumSource

src = ContinuumSource(alloy="Al3Cu5Nb01Fe")          # ou "Al08Si06Mg02Fe"
state = src.state_at_gradient(src.grad_T_exp)        # 3294 K/m
model = FerreiraModel()
for order in (1, 2):
    print(order, model.r_critical(state, order=order), model.gamma_tensor(state, order=order))
# 1  1.848e-06 m   3.887e-07 m·K
# 2  5.081e-06 m   1.068e-06 m·K
```

Convenções de sinal, nomes e opções de cada equação estão em [`docs/CONVENTIONS.md`](docs/CONVENTIONS.md).

## Reproduzir os resultados

| Resultado | Comando | Saída |
|---|---|---|
| Relatório E1 (raízes, CNT × Ferreira, barreiras) | `python scripts/report_e1.py` | `results/e1/` |
| Verificação do contínuo contra os scripts originais | `python scripts/verify_continuum.py` | `docs/CONTINUUM_VERIFICATION.md` |
| SDAS (MRB/RB) × dados experimentais | `python scripts/report_sdas.py` | `results/sdas/` |
| Figuras Γ × ∇**T** e ΔT × r_c | `python scripts/plot_scales.py` | `results/figures/` |
| CALPHAD Al-Nb | `python scripts/calphad_liquidus.py` | `results/calphad/liquidus.md` |
| CALPHAD das ligas (T_L, T_S, ΔH) | `python scripts/calphad_alloys.py` | `results/calphad/alloys_TL_dH.md` |
| Tabelas contínuo × MD (Word) | `python scripts/make_tables_docx.py` | `results/tables/` |

Os scripts de CALPHAD leem as bases TDB de uma pasta local do OpenCalphad (`OCDIR`/`TDB` no início de cada script);
ajuste o caminho para a sua máquina. A base COST507 é usada **sem a fase ALTI** (ver o docstring de `calphad_alloys.py`).

## Simulações de MD

**Software:** LAMMPS 10Sep2025 com PLUMED 2.9.4 (CPU e pacote GPU/CUDA, precisão mista) e GPUMD (para o NEP89),
em WSL2/Ubuntu 22.04. Os scripts de [`install/`](install) compilam tudo em `~/opt`, sem sudo:

```bash
bash install/build_lammps_plumed_cuda_wsl.sh     # GPU (CUDA); CPU: build_lammps_plumed_wsl.sh
source ~/opt/nucmd_cuda_env.sh
lmp -sf gpu -pk gpu 1 -in in.xxx
```

**Protocolos (E2):**

| Diretório | Conteúdo |
|---|---|
| [`md/e2_melting/`](md/e2_melting) | T_m por coexistência isotérmica (`in.coexist`), ΔH_m e densidades (`in.enthalpy`), γ₀ por flutuação capilar (`in.cfm`, `cfm_morris.py`) — Borovikov 2024 |
| [`md/e2_nep89/`](md/e2_nep89) | mesmos protocolos portados para GPUMD com o NEP89 |
| [`md/e2_cfm_thick/`](md/e2_cfm_thick) | γ₀ com fitas espessas (1 ns) e reanálise comum das 9 fitas com piso de ruído (`fit_floor.py`) |
| [`md/e2_substrate/`](md/e2_substrate) | estabilidade do Al₃Nb (D0₂₂ × L1₂) por potencial |
| [`md/gpu_benchmark/`](md/gpu_benchmark) | desempenho na RTX 4050: EAM ≈ 8,7–11 Matom·passo/s (3–4,6× 8 processos de CPU) |

As trajetórias (`*.dump`, `prod.xyz`, ~830 MB) **não estão versionadas**; os scripts de cada diretório as regeneram.
Armadilhas de protocolo já encontradas (termostatos com grupos congelados, CNA a T_m, periodicidade da fita (111),
espessura mínima no GPUMD…) estão listadas em [`HANDOFF.md`](HANDOFF.md) §7.

## Estrutura do repositório

```
├── nucleation_md/                 pacote Python
│   ├── models.py                  NucleationModel, CNT, Tolman, FerreiraModel
│   ├── sources/continuum.py       ContinuumSource: estados a partir das tabelas do contínuo
│   ├── sdas.py                    modelos MRB/RB de SDAS
│   ├── surface.py                 ΔG_S (métrica/Ricci e Gurtin–Murdoch)
│   └── validation.py              critério de aceitação do potencial
├── tests/                         pytest (77 testes)
├── scripts/                       relatórios, figuras, CALPHAD, tabelas
├── md/                            entradas e análises de MD por etapa
├── potentials/                    potenciais interatômicos (terceiros)
├── calphad_oc/                    macro do OpenCalphad (OC6P)
├── install/                       build de LAMMPS + PLUMED (CPU/CUDA) e GPUMD no WSL
├── data/                          tabelas do contínuo e pontos experimentais de SDAS
├── results/                       saídas (E1, SDAS, CALPHAD, figuras, tabelas, propriedades de MD)
├── docs/                          convenções e verificação do contínuo
└── Continuum_Mechanics_Solution/  scripts originais do contínuo, artigos e apresentação (não modificados)
```

## Documentação

- [`HANDOFF.md`](HANDOFF.md) — estado detalhado, decisões, questões em aberto e próximos passos.
- [`docs/CONVENTIONS.md`](docs/CONVENTIONS.md) — sinais, nomes e equações implementadas.
- [`docs/CONTINUUM_VERIFICATION.md`](docs/CONTINUUM_VERIFICATION.md) — reprodução do contínuo coluna a coluna.
- [`md/e2_melting/E2_results.md`](md/e2_melting/E2_results.md), [`md/e2_nep89/E2_nep89_results.md`](md/e2_nep89/E2_nep89_results.md),
  [`md/e2_cfm_thick/E2_thick_results.md`](md/e2_cfm_thick/E2_thick_results.md) — validação dos potenciais.

## Como citar

Modelo:

> I. L. Ferreira, *Physica B: Condensed Matter* (2024), [doi:10.1016/j.physb.2024.416494](https://doi.org/10.1016/j.physb.2024.416494).

Se usar os potenciais, cite também as fontes originais (ver abaixo).

## Licença e conteúdo de terceiros

O código deste repositório está sob a licença [MIT](LICENSE). Não estão cobertos por ela:

- `potentials/`: potenciais interatômicos dos respectivos autores — Borovikov, Mendelev et al., *Int. J. Plast.* 178
  (2024) 104004; Mendelev et al. 2008; Jelinek et al. 2012; Fereidonnejad et al. 2022 (via
  [NIST IPR](https://www.ctcms.nist.gov/potentials/)); NEP89 — Liang et al. 2025 (arXiv:2504.21286). Valem os termos originais.
- `Continuum_Mechanics_Solution/Paper_and_presentation/`: artigos e apresentação, sob os direitos das editoras e dos autores.
