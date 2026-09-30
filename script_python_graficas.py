## @file script_python_graficas.py
#  @brief Analiza la densidad de PIMC del oscilador armónico 1D.
#
#  Este script compara la distribución de posiciones obtenida mediante PIMC
#  con la solución analítica del oscilador armónico cuántico a temperatura
#  finita. También calcula el radio cuadrático medio y estima las energías
#  cinética, potencial y total usando el teorema del virial.
#
#  @details
#  El programa MC guarda la distribución radial en `outrd.dat`. En 1D, esta
#  salida representa la rama x >= 0 de una distribución simétrica y contiene
#  varios bloques de Monte Carlo consecutivos. El script separa los bloques,
#  comprueba que tienen la misma malla y calcula su promedio para reducir el
#  ruido estadístico.
#
#  Para una partícula en un oscilador armónico con hbar = m = omega = 1, la
#  predicción analítica usada es:
#  @f[
#  n(x) = \left(\frac{\tanh(\beta/2)}{\pi}\right)^{D/2}
#         \exp\left[-x^2\tanh(\beta/2)\right].
#  @f]
#
#  El valor exacto del radio cuadrático medio y de la energía total es:
#  @f[
#  \langle x^2\rangle = E = \frac{D}{2\tanh(\beta/2)}.
#  @f]
#
#  En este problema se cumple Ekin = Epot = <x^2>/2. La salida `oute.dat`
#  no se utiliza para separar estas contribuciones, ya que el estimador
#  empleado por esta configuración de PIMC no guarda una cinética independiente.
#
#  @pre `outrd.dat` debe estar en el directorio de ejecución.
#  @post Se imprime la comparación entre simulación y teoría, se muestra una
#        gráfica de la densidad y se muestra una gráfica de residuos con barras
#        de error estándar.

import numpy as np
import matplotlib.pyplot as plt

## Número de dimensiones del sistema.
D = 1.0
## Temperatura del sistema en unidades reducidas.
T = 1.0
## Temperatura inversa, beta = 1/T.
beta = 1.0 / T

## Leer la distribución radial 1D generada por PIMC.
##
## Cada bloque de Monte Carlo contiene dos columnas: posición y densidad.
## El programa concatena los bloques en el mismo archivo.
data = np.loadtxt('outrd.dat')
## Índices donde empieza un bloque nuevo; la posición vuelve a disminuir.
saltos = np.where(np.diff(data[:, 0]) < 0)[0] + 1
## Lista de histogramas individuales, uno por bloque de Monte Carlo.
bloques = np.split(data, saltos)

## Comprobar que todos los bloques tienen el mismo número de bins.
if any(len(bloque) != len(bloques[0]) for bloque in bloques):
    raise ValueError('Los bloques de outrd.dat no tienen el mismo tamaño')

## Coordenadas de los centros de bin de la distribución simulada.
x_sim = bloques[0][:, 0]
## Comprobar que todos los bloques usan los mismos centros de bin.
if not all(np.allclose(bloque[:, 0], x_sim) for bloque in bloques[1:]):
    raise ValueError('Los centros de bin cambian entre bloques')

## Densidades de cada bloque, usadas para calcular el promedio y su error.
n_bloques = np.array([bloque[:, 1] for bloque in bloques])
## Densidad simulada promedio sobre todos los bloques.
n_sim = np.mean(n_bloques, axis=0)
## Error estándar del promedio en cada centro de bin.
n_error = np.std(n_bloques, axis=0, ddof=1) / np.sqrt(len(bloques))

## Coordenadas empleadas para dibujar la solución analítica.
x_teo = np.linspace(x_sim.min(), x_sim.max(), 500)
## Perfil de densidad analítico del oscilador armónico a temperatura T.
n_teo = (np.tanh(beta / 2.0) / np.pi)**(D / 2.0) * np.exp(-(x_teo**2) * np.tanh(beta / 2.0))

## Radio cuadrático medio teórico, <x^2>.
r2_teo = D / (2.0 * np.tanh(beta / 2.0))
## Energía cinética teórica; en este sistema coincide con la potencial.
e_kin_teo = r2_teo / 2.0
## Energía potencial teórica del potencial V(x) = x^2/2.
e_pot_teo = r2_teo / 2.0
## Energía total teórica, E = Ekin + Epot.
e_teo = e_kin_teo + e_pot_teo

## Anchura de los bins de la distribución simulada.
dx = x_sim[1] - x_sim[0]
## Estimación de <x^2>; se multiplica por 2 porque solo se guarda x >= 0.
r2_sim = 2.0 * np.sum(x_sim**2 * n_sim) * dx

## Energía potencial estimada a partir de la densidad simulada.
e_pot_sim = r2_sim / 2.0
## Energía cinética estimada mediante el teorema del virial.
e_kin_sim = r2_sim / 2.0
## Energía total estimada, suma de cinética y potencial.
e_sim = e_kin_sim + e_pot_sim

## Mostrar resultados numéricos de la simulación y de la teoría.
print(f'Ekin simulacion = {e_kin_sim:.6f}')
print(f'Epot simulacion = {e_pot_sim:.6f}')
print(f'E total simulacion = {e_sim:.6f}')
print(f'Ekin teoria = {e_kin_teo:.6f}')
print(f'Epot teoria = {e_pot_teo:.6f}')
print(f'E total teoria = {e_teo:.6f}')
print(f'<x^2> simulacion = {r2_sim:.6f}')
print(f'<x^2> teoria = {r2_teo:.6f}')

## Residuo entre la densidad simulada y la solución teórica en los mismos bins.
n_teo_sim = (np.tanh(beta / 2.0) / np.pi)**(D / 2.0) * np.exp(
    -(x_sim**2) * np.tanh(beta / 2.0)
)
residuo = n_sim - n_teo_sim
## RMS del residuo como medida global de la discrepancia simulación-teoría.
residuo_rms = np.sqrt(np.mean(residuo**2))
print(f'Residuo RMS = {residuo_rms:.6f}')

## Dibujar la distribución analítica y la simulada.
plt.plot(x_teo, n_teo, 'r-', linewidth=2, label='Teoría (Fórmula original)')
## Cada punto azul representa un centro de bin de la densidad promediada.
plt.plot(x_sim, n_sim, 'b.', markersize=5, label=f'Simulación PIMC ({len(bloques)} bloques)')
## Etiquetas y formato de la figura.
plt.title('Densidad Cartesiana del Oscilador Armónico 1D')
plt.xlabel('Posición (z)')
plt.ylabel('Densidad n(z)')
plt.legend()
plt.grid(True)

## Dibujar el error de la simulación respecto a la teoría.
## Las barras representan el error estándar del promedio de los bloques.
plt.figure()
plt.axhline(0.0, color='black', linewidth=1)
plt.errorbar(x_sim, residuo, yerr=n_error, fmt='g.', markersize=5,
             capsize=2, label='Residuo +/- error estándar')
plt.title('Residuo de la densidad PIMC respecto a la teoría')
plt.xlabel('Posición (z)')
plt.ylabel(r'$n_{PIMC}(z) - n_{teoria}(z)$')
plt.legend()
plt.grid(True)

plt.show()