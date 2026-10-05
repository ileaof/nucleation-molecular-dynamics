@echo off
rem Benchmark on the RTX 4050 with the Windows LAMMPS binary (GPU package, OpenCL, mixed precision).
cd /d %~dp0
for %%N in (40 63 100) do (
  lmp -sf gpu -pk gpu 1 -in in.bench -var n %%N -log bench_gpu_n%%N.log
)
findstr /C:"Performance" bench_gpu_n*.log
pause
