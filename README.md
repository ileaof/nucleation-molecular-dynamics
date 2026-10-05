# Nucleation_MD — modelo de nucleação de Ferreira (2024): contínuo × dinâmica molecular

Projeto para calcular atomisticamente (LAMMPS + PLUMED + Python) as grandezas do modelo de
Ferreira (2024, Physica B, doi:10.1016/j.physb.2024.416494) e validar a solução do contínuo.

> **Trocando de máquina? Leia primeiro [`HANDOFF.md`](HANDOFF.md)** (§6.0: o que já está instalado na máquina GPU).

## Estrutura

| Pasta | Conteúdo |
|---|---|
| `Continuum_Mechanics_Solution/` | Código e tabelas do contínuo (do autor; **não é modificado**) |
| `nucleation_md/` | Pacote Python: `CNT`, `Tolman`, `FerreiraModel`, fontes contínuo/MD, ΔG_S, SDAS |
| `tests/` | `python -m pytest` (77 testes) |
| `scripts/` | relatórios (`report_e1.py`, `report_sdas.py`), figuras (`plot_scales.py`), CALPHAD (`calphad_liquidus.py`), contínuo (`run_continuum.py`, `verify_continuum.py`) |
| `docs/` | `CONVENTIONS.md`, `CONTINUUM_VERIFICATION.md` |
| `data/` | tabelas do contínuo da 2ª liga e pontos experimentais de SDAS |
| `results/` | saídas (E1, SDAS, figuras, CALPHAD) |
| `potentials/` | potenciais interatômicos usados |
| `md/` | entradas e análises do LAMMPS por etapa |
| `calphad_oc/` | macro do OC6P (Al-Nb em wt.%) |
| `install/` | scripts de instalação do LAMMPS + PLUMED no WSL (CPU e CUDA) |

## Estado (2026-10-03)

- **E1 (Python)**: concluída. Reproduz o contínuo das duas ligas (Γ¹ˢᵗ, Γ²ⁿᵈ, ∇T, SDAS dos slides 25–30).
- **CALPHAD**: o α-Al nasce no peritético a 934,51 K de um líquido com 0,0062 wt.% Nb; o Nb está no Al₃Nb.
  Decisão: líquido de Al, núcleo Al-0,5 wt.% Nb, substrato Al₃Nb (opção a).
- **Potencial escolhido**: `potentials/NiAlNb_Borovikov2024/Ni-Al-Nb.eam.fs` (EAM/FS).
  Al₃Nb D0₂₂ estável (ΔH_f = −0,31 eV/átomo); substrato a ser congelado na rede experimental.
  Rejeitados: Al-Nb MEAM 2022 (Al₃Nb instável), Jelinek AlSiMgCuFe (T_m < 900 K).
- **E2 (validação do potencial)**: T_m(Al) = 935,4 K; ΔH_m = 325,5 kJ/kg; ρ_s = 2551,5 e ρ_l = 2397,9 kg/m³
  e γ₀ = 0,106 ± 0,004 J/m² (flutuação capilar, `cfm_morris.py`) (`md/e2_melting/E2_results.md`).
  Pendentes: as tolerâncias do critério (a decidir; γ₀ de referência 0,296 J/m² é do contínuo em r₀ = 6,3 μm).

## Na máquina nova (com a RTX 4050)

1. **Python (Windows)**: `pip install numpy scipy pandas matplotlib pytest pycalphad`; depois `python -m pytest` na raiz.
2. **Teste de GPU sem compilar**: `md/gpu_benchmark/run_gpu_windows.bat` (LAMMPS para Windows, pacote GPU via OpenCL).
3. **LAMMPS + PLUMED no WSL**:
   - CPU: `bash install/build_lammps_plumed_wsl.sh` → `source ~/opt/nucmd_env.sh`
   - GPU (CUDA): `bash install/build_lammps_plumed_cuda_wsl.sh` → `source ~/opt/nucmd_cuda_env.sh`
   As instalações ficam em `~/opt` do WSL e **não** vão junto com esta pasta.
4. **Caminhos fixos a conferir** se o usuário/OneDrive mudar:
   - `scripts/calphad_liquidus.py` e `calphad_oc/alnb_wt.OCM` apontam para `C:\Users\ileao\OneDrive\Documentos\OpenCalphad\OC6` (bases TDB fora do projeto).
   - Os scripts `.sh` usam `/mnt/c/Users/ileao/OneDrive/Documentos/Nucleation_MD_GPU/...`.
5. Desempenho: neste laptop (i7-1360P), 8 processos MPI foram o ótimo por simulação; EAM ≈ 2×10⁵ átomo·passo/s com 4 processos.

## Licença e conteúdo de terceiros

O código deste repositório está sob a licença MIT (`LICENSE`). Não estão cobertos por ela:
- `potentials/`: potenciais interatômicos dos respectivos autores (Borovikov et al. 2024, Mendelev 2008, Jelinek 2012,
  Fereidonnejad 2022 — via NIST IPR; NEP89 — Liang et al. 2025), sob os termos originais; citar as fontes ao usar.
- `Continuum_Mechanics_Solution/Paper_and_presentation/`: artigos e apresentação, sob os direitos das editoras/autores.

Não versionado (ver `.gitignore`): trajetórias de MD (`*.dump`, `prod.xyz`), regeneráveis pelos scripts de `md/`.
