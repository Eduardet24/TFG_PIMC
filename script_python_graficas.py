import numpy as np
import matplotlib.pyplot as plt

# Parámetros
D = 1.0
T = 1.0 
beta = 1.0 / T

# 1. Leer la distribución radial 1D. El programa añade un bloque por 
# cada bloque de Monte Carlo, por lo que el archivo contiene varios histogramas.
data = np.loadtxt('outrd.dat')
saltos = np.where(np.diff(data[:, 0]) < 0)[0] + 1
bloques = np.split(data, saltos)

if any(len(bloque) != len(bloques[0]) for bloque in bloques):
    raise ValueError('Los bloques de outrd.dat no tienen el mismo tamaño')

x_sim = bloques[0][:, 0]
if not all(np.allclose(bloque[:, 0], x_sim) for bloque in bloques[1:]):
    raise ValueError('Los centros de bin cambian entre bloques')

# Promediar los bloques reduce el ruido estadístico del último bloque aislado.
n_sim = np.mean([bloque[:, 1] for bloque in bloques], axis=0)

# 3. Fórmula teórica original del profesor (Campana centrada en 0)
x_teo = np.linspace(x_sim.min(), x_sim.max(), 500)
n_teo = (np.tanh(beta / 2.0) / np.pi)**(D / 2.0) * np.exp(-(x_teo**2) * np.tanh(beta / 2.0))

# 4. Graficar
plt.plot(x_teo, n_teo, 'r-', linewidth=2, label='Teoría (Fórmula original)')
plt.plot(x_sim, n_sim, 'b.', markersize=5, label=f'Simulación PIMC ({len(bloques)} bloques)')
plt.title('Densidad Cartesiana del Oscilador Armónico 1D')
plt.xlabel('Posición (z)')
plt.ylabel('Densidad n(z)')
plt.legend()
plt.grid(True)
plt.show()