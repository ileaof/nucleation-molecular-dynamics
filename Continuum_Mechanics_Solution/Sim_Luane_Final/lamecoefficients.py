def calcular_coeficientes_lame(E, nu):
    """
    Calcula os Coeficientes de Lamè (lambda e mu) a partir do 
    Módulo de Young (E) e Coeficiente de Poisson (nu).
    
    Args:
        E (float): Módulo de Young (unidade de tensão, ex: Pa)
        nu (float): Coeficiente de Poisson (adimensional)
        
    Returns:
        tuple: (lambda, mu)
    """
    # Segundo Coeficiente de Lamè (mu - Módulo de Cisalhamento)
    mu = E / (2 * (1 + nu))
    
    # Primeiro Coeficiente de Lamè (lambda)
    lame_lambda = (E * nu) / ((1 + nu) * (1 - 2 * nu))
    
    return lame_lambda, mu

# --- Exemplo de Uso ---
# Exemplo: Aço (E = 200 GPa, nu = 0.3), Alumínio (E = 72.4 GPa, nu = 0.33), gelo (E = 9 GPa, nu = 0.33), Nb (E = 105 GPa, nu = 0.4))
modulo_young = 105e9 #9e9 #72.4e9  # Pa
poisson = 0.4  #0.33

lmbda, mu = calcular_coeficientes_lame(modulo_young, poisson)

print(f"Módulo de Young (E): {modulo_young:.2e} Pa")
print(f"Coeficiente de Poisson (nu): {poisson}")
print("-" * 30)
print(f"1º Coeficiente de Lamè (λ): {lmbda:.2e} Pa")
print(f"2º Coeficiente de Lamè (μ/Cisalhamento): {mu:.2e} Pa")
