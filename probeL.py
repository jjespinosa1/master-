import matplotlib.pyplot as plt
import numpy as np
from scipy.special import genlaguerre
x = np.arange(-4.0, 12.0, 0.01)
fig, ax = plt.subplots()

ax.axhline(0, color='red', linestyle='--')
ax.set_ylim(-5.0, 10.0)
ax.set_title(r'Generalized Laguerre polynomials $L_3^{\alpha}$')
for alpha in np.arange(0, 1):
    ax.plot(x, genlaguerre(11, alpha)(x), label=rf'$L_3^{(alpha)}$')
plt.legend(loc='best')
plt.show()