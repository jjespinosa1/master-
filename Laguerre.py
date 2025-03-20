import numpy as np
from scipy.special import eval_genlaguerre

def laguerre_function(n, alpha, x):
    """
    Evalúa el polinomio de Laguerre generalizado Lf(n, alpha, x).

    Parameters:
    n (int): El orden del polinomio.
    alpha (float): El parámetro alpha.
    x (float): El punto donde se evalúa el polinomio.

    Returns:
    float: El valor del polinomio de Laguerre generalizado en x.
    """
    return eval_genlaguerre(n, alpha, x)

# Parámetros evaluados en Laguerre.f
n = 5
alpha = 0.5
x = 0.25

# Evaluar el polinomio de Laguerre
fx = laguerre_function(n, alpha, x)
print(f"Lf({n},{alpha},{x}) = {fx:.6f}")