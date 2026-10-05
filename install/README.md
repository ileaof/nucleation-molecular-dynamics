# Instalação LAMMPS + PLUMED (WSL Ubuntu-22.04)

Construído por `build_lammps_plumed_wsl.sh` (sem sudo, tudo em `~/opt`).

| Item | Versão / local |
|---|---|
| LAMMPS | `patch_10Sep2025` — mesma versão do executável Windows |
| PLUMED | 2.9.4 (versão fixada pelo LAMMPS 10Sep2025, MD5 conferido), `--enable-modules=all` (crystallization, multicolvar, adjmat) |
| Modo PLUMED | `runtime` (LAMMPS carrega `libplumedKernel.so` em tempo de execução) |
| MPI / OpenMP | OpenMPI 4.1.2 (sistema) / sim |
| Pacotes | COLVARS EXTRA-COMPUTE EXTRA-DUMP EXTRA-FIX EXTRA-PAIR MANYBODY MC MEAM MISC MOLECULE OPENMP PLUMED REPLICA |
| Python | venv `~/opt/venvs/nucmd` (numpy, scipy, matplotlib, pandas, pytest, módulo `lammps`) |

Uso:

```bash
wsl -d Ubuntu-22.04
source ~/opt/nucmd_env.sh
mpirun -np 8 lmp -in in.xxx          # LAMMPS com fix plumed
plumed sum_hills --hills HILLS ...   # CLI do PLUMED
```

## Teste de fumaça (`smoke_test/`)

- `in.smoke` (Al FCC, `Al_mm.eam.fs`, 500 átomos, 100 passos NVE — só teste de instalação, não é escolha de potencial).
- Windows × WSL (1 processo): saída termodinâmica **bit-idêntica**.
- WSL com `fix plumed` (Q6 → LOCAL_Q6 → MFILTER_MORE → CONTACT_MATRIX → DFSCLUSTERING → CLUSTER_NATOMS):
  maior cluster sólido = 500 átomos (cristal todo), e a dinâmica é idêntica à sem PLUMED (só PRINT, sem viés).

## Observação de reprodutibilidade

`create_atoms` atribui IDs conforme a decomposição MPI. Com N processos diferentes, a mesma semente de
`velocity create` cai em sítios diferentes → trajetórias diferentes (fisicamente equivalentes).
**Regra do projeto:** gerar cada configuração inicial uma vez (`write_data`) e reutilizá-la com `read_data`.

## Notas de build

- `python3 -m venv` do Ubuntu falha sem `python3.10-venv` (sudo) → venv criado com `--without-pip` + `get-pip.py`.
- `PLUMED_MODE=shared` falhou no link (`libplumedKernel.so` não encontrado); `static` exige BLAS/LAPACK/GSL do sistema; `runtime` funciona sem dependências.
