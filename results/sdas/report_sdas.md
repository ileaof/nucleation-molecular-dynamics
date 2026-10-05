# Cadeia de validação contínuo → SDAS

∇T experimental → Γ (FerreiraModel, Γ = A·∇T) → SDAS (MRB, Ferreira 2024) → pontos experimentais.

Métricas: rms_log = √⟨ln²(pred/exp)⟩; bias_log = ⟨ln(pred/exp)⟩; within = fração dos pontos cuja barra de erro contém a previsão.


## Al3Cu5Nb01Fe — ∇T = 3294.0 K/m (Mendes et al. (2023))

r_c² = 5.0809 μm, ΔT = 0.42059 K, Γ¹ˢᵗ = 3.88679e-07 m·K, Γ²ⁿᵈ = 1.06849e-06 m·K

| modelo | Γ [m K] | SDAS = a·t^b | rms_log | bias_log | within |
|---|---|---|---|---|---|
| MRB Γ¹ˢᵗ | 3.8868e-07 | 7.67·t^0.377 | 0.054 | -0.000 | 100% |
| MRB Γ²ⁿᵈ | 1.0685e-06 | 9.38·t^0.379 | 0.219 | +0.213 | 94% |
| RB Γ²ⁿᵈ | 1.0685e-06 | 27.59·t^0.333 | 1.078 | +1.074 | 0% |

## Al08Si06Mg02Fe — ∇T = 7902.38 K/m (Marques et al. (2025))

r_c² = 3.9821 μm, ΔT = 0.79077 K, Γ¹ˢᵗ = 7.63677e-07 m·K, Γ²ⁿᵈ = 1.57447e-06 m·K

| modelo | Γ [m K] | SDAS = a·t^b | rms_log | bias_log | within |
|---|---|---|---|---|---|
| MRB Γ¹ˢᵗ | 7.6368e-07 | 5.81·t^0.387 | 0.054 | +0.027 | 100% |
| MRB Γ²ⁿᵈ | 1.5745e-06 | 6.71·t^0.387 | 0.180 | +0.173 | 75% |
| RB Γ²ⁿᵈ | 1.5745e-06 | 27.03·t^0.333 | 1.392 | +1.391 | 0% |