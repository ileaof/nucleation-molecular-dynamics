# Teste de velocidade — RTX 4050 (outro computador)

Átomos: n = 40 → 256 000; n = 63 → 1 000 188; n = 100 → 4 000 000 (limite de memória de 6 GB).

1. **Sem compilar nada (Windows, OpenCL):** dê duplo clique em `run_gpu_windows.bat`
   (o LAMMPS 10Sep2025 para Windows já tem o pacote GPU com OpenCL).
2. **Build CUDA + PLUMED no WSL2:** `install/build_lammps_plumed_cuda_wsl.sh`
   (pré-requisitos no cabeçalho do script) e depois
   `source ~/opt/nucmd_cuda_env.sh; lmp -sf gpu -pk gpu 1 -in in.bench -var n 63`.

Linha a comparar em cada log: `Performance: ... katom-step/s`.
Referência neste laptop (CPU, 4 processos, EAM): ~200 katom-step/s.

## Resultados na máquina GPU (2026-10-04) — Borovikov 2024, Al líquido 1200 K, NPT, Δt = 2 fs

Build CUDA do WSL (`run_gpu_wsl.sh`, `run_small_wsl.sh`), precisão mista; CPU = i7-13650HX, 8 processos MPI.
Valores da 2ª corrida (`run ${steps}`) de cada log.

| átomos | GPU [Matom·passo/s] | GPU [ns/dia] | CPU 8 proc. [Matom·passo/s] | GPU/CPU |
|---|---|---|---|---|
| 8 788 (n = 13) | 8,7 | 170 | 2,8 | 3,1 |
| 32 000 (n = 20) | 10,5 | 56 | 2,7 | 3,9 |
| 256 000 (n = 40) | 11,2 | 7,6 | 2,5 | 4,6 |
| 1 000 188 (n = 63) | 10,6 | 1,8 | — | — |
| 4 000 000 (n = 100) | 7,5–10,3 | 0,3–0,45 | — | — |

- A GPU ocupa 1 núcleo de CPU: dá para rodar uma corrida GPU e outra CPU (≤ 8 proc.) ao mesmo tempo.
- 4 M átomos cabem nos 6 GB, mas a 2ª corrida caiu 27 % (provável limitação térmica/potência do laptop).
- Máquina antiga (i7-1360P, 4 proc.): ~0,2–0,8 Matom·passo/s.
