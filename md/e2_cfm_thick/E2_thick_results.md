# γ₀ do Borovikov com fitas espessas (HANDOFF §4, opção 2) e reanálise comum das 9 fitas

Protocolo de `md/e2_melting/in.cfm` (cópia com `ny` livre: `in.cfm`), LAMMPS GPU, produção NPH 1 ns (antes 400 ps).
Fitas: (100) W = 33,0 Å (ny = 8), (110) 29,2 Å (ny = 5), (111) 40,4 Å (ny = 6); antes 12,4 / 11,7 / 20,2 Å.
T média da produção 928 K nas três; fração sólida e L_z constantes em 1 ns (interfaces estáveis).
Dumps em `~/runs/cfm_thick` (Ubuntu-22.04). Espectros completos das 9 fitas em `spectra/<conjunto>/` (`spectra_all.sh`).

## 1. Método da E2 (janela 4–10, sem piso)

| fita espessa | inclinação log-log | γ̃ [J/m²] | E2 (fina) |
|---|---|---|---|
| (100) | −2,28 | 0,0765 ± 0,0025 | 0,094 |
| (110) | −2,41 | 0,0992 ± 0,0051 | 0,096 |
| (111) | −1,20 | 0,0954 ± 0,0088 | 0,116 |

γ₀ = 0,095 ± 0,005 J/m² — mas a (111) com inclinação −1,2 não está no regime capilar; não usar.

## 2. Ajuste com piso de ruído, ⟨|A|²⟩ = k_BT/(W L γ̃ k²) + c (`fit_floor.py`), mesmos critérios para as 9 fitas

| janela de modos | Borovikov fina (E2) | Borovikov espessa | NEP89 |
|---|---|---|---|
| 2–15 | 0,112 ± 0,004 | 0,116 ± 0,004 | 0,119 ± 0,004 |
| 2–20 | 0,116 ± 0,004 | 0,120 ± 0,004 | 0,125 ± 0,004 |
| 3–12 | 0,105 ± 0,004 | 0,112 ± 0,004 | 0,115 ± 0,004 |
| 1–10 (modo 1 não relaxado; χ² alto) | 0,113 ± 0,006 | 0,134 ± 0,011 | 0,169 ± 0,019 |

Tabelas por orientação: `gamma0_floor_k<kmin>-<kmax>.md`.
Observações: o piso c sai **negativo** em (100) e (110) (espectro mais íngreme que k⁻², não ruído branco) e positivo em (111);
χ²_red de (110) é 4–9. O modelo de piso não descreve todos os espectros — a dispersão entre janelas é o melhor indicador do sistemático.

## Leitura

- Engrossar as fitas e alongar a produção **não muda o γ₀ do Borovikov além do sistemático da análise**:
  0,105–0,120 J/m² (fina) × 0,112–0,120 (espessa), janelas sem o modo 1. O 0,106 da E2 se sustenta dentro de ~±10 %.
- Com a mesma análise o NEP89 dá 0,115–0,125 J/m²: **a diferença entre os potenciais (~5 %) é menor que o sistemático**.
  γ₀ não discrimina os dois potenciais.
- Ambos ficam abaixo do Al experimental (0,13–0,17) e dos valores do contínuo (0,183; 0,296 em r₀ = 6,3 μm).
