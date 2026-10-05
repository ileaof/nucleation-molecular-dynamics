def midpoint(f, a, b, n, rhos, cps, TL, TF, beta, k):
    """This function integrates only the Energy as a function of solid fraction."""
    h = (b-a)/n
    result = 0.0
    for i in range(n):
        result += f((a + h/2.0) + i*h, rhos, cps, TL, TF, beta, k)
    result *= h
    return result

