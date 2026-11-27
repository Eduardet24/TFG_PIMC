# Simulación PIMC de Osciladores Armónicos 1D y 2D

Entorno de trabajo para simulaciones de **Monte Carlo Cuántico (QMC)**, centrado en **Path Integral Monte Carlo (PIMC)** para osciladores armónicos en 1D y 2D. Incluye un **algoritmo Worm** con muestreo de permutaciones bosónicas, estimadores de densidad, energía, OBDM y superfluidez, y un benchmark canónico 2D para $N > 1$.

El código base en C fue desarrollado originalmente por G. E. Astrakharchik (2001–2017) y ha sido adaptado y ampliado para el Trabajo de Fin de Grado (TFG).

---

## 1. Estructura del Repositorio

```text
Grigori_codigo/
├── MC                      # Binario compilado
├── MC.cfg                  # Configuración principal de la simulación
├── *.c, *.h                # Código fuente en C (MC, PIMC, DMC, VMC, Worm, etc.)
├── Makefile.txt / .AMD     # Makefiles para compilación (icc / AMD)
├── scritpt_Rdz_grigori     # Script de graficación rápida con gnuplot
├── scripts/                # Scripts de análisis, verificación y campañas en Python
└── Resultados/             # Salidas organizadas por estudio/campaña
```

---

## 2. Requisitos y Compilación

- **Compilador C:** `gcc` con OpenMP (recomendado) o Intel C Compiler (`icc`).
- **Python 3.8+:** `pip install numpy matplotlib scipy`
- **Opcionales:** `gnuplot`, `doxygen`.

Compilación desde la raíz (`Grigori_codigo/`):

```bash
# GCC con OpenMP (recomendada)
gcc -O2 -fopenmp *.c -lm -o MC

# Alternativa con Intel (icc)
make -f Makefile.txt
```

---

## 3. Configuración (`MC.cfg`) y Salidas

> [!IMPORTANT]
> **Sintaxis estricta:** los parámetros deben escribirse sin espacio antes del signo `=` (ej. `T= 0.5`, **no** `T = 0.5`).

| Parámetro | Valor típico | Descripción |
| :--- | :--- | :--- |
| `MC=` | `3` | Método: `0` (VMC), `1` (DMC), `3` (PIMC). |
| `Npop=` | `5` .. `30` | Número de *beads* o réplicas temporales ($N_{\text{beads}}$). |
| `T=` | `0.5` .. `4.0` | Temperatura en unidades reducidas ($\hbar = m = \omega = k_B = 1$). |
| `N=` | `1` .. `5` | Número de partículas. |
| `omega_z=` | `1` | Frecuencia del potencial armónico. |
| `acceptance_rate=` | `50` | Tasa de aceptación objetivo (%). Con `0` usa paso fijo `dt_vmc`. |
| `blck_heating=` / `blck=` | `100` / `1000` | Bloques de termalización y de medida. |
| `Niter=` / `Nmeasure=` | `1000` / `1` | Iteraciones por bloque y frecuencia de muestreo. |
| `measure_energy=` / `measure_RadDistr=` | `1` / `1` | Activa el guardado de energía (`oute.dat`) y densidad (`outrd.dat`). |

**Salidas principales de `./MC`:**
- `oute.dat`: energías por bloque (`# Etot Epot Ekin`).
- `outrd.dat`: distribución espacial/radial $n(r)$.
- `MC.log`: registro de ejecución y tasas de aceptación.
- `outworm_*.dat`: estadísticas del algoritmo Worm (ciclos, permutaciones, *gap*, OBDM, *snapshots*).

---

## 4. Uso Rápido

Ejecutar siempre los comandos desde la raíz del repositorio (`Grigori_codigo/`).

### Simulación individual
```bash
./MC
python3 scripts/Script_Armonico1D_basico.py
```

### Campaña automatizada por lotes (1D)
```bash
# Previsualizar el plan sin ejecutar simulaciones (dry-run)
python3 scripts/Script_Lote_PIMC.py

# Ejecutar matriz por defecto o personalizada
python3 scripts/Script_Lote_PIMC.py --run
python3 scripts/Script_Lote_PIMC.py --run --temperaturas 0.5,1.0,2.0 --beads 5,10,15,20,25,30 --campaign lote_modos_normales_blck1000
```
Las ejecuciones completadas se guardan en `Resultados/1DArmonico/AR<rate>/<campaña>/Temp<T>/N<beads>/` con un marcador `.complete` para evitar repetir simulaciones ya terminadas.

### Benchmark 2D ($N > 1$)
```bash
python3 scripts/Script_Comparacion_N_Particulas_2D.py --run --particles 5
```

---

## 5. Scripts de Análisis y Verificación (`scripts/`)

| Script | Propósito |
| :--- | :--- |
| `Script_Lote_PIMC.py` | Automatización de barridos en $(T, N_{\text{beads}})$. |
| `Script_Armonico1D_basico.py` | Densidad $n(x)$, residuos y energías por virial de la corrida actual. |
| `Script_EPvsBeads` | Convergencia de $E_{\text{pot}}$ frente a $1/N_{\text{beads}}$ y extrapolación al continuo. |
| `Script_EPvsT` | Comparación de $E_{\text{pot}}(T)$ frente a los modelos clásico y cuántico exacto. |
| `Script_Convergencia_AR0.py` | Ajustes cuadráticos de convergencia para series `AR0` / `AR50`. |
| `Script_Caminos_Armonico_2D.py` | Visualización de caminos cerrados y fluctuaciones del centroide en 2D. |
| `Script_Comparacion_N_Particulas_2D.py` | Benchmark $N>1$ en 2D frente a curvas analíticas distinguible y Bose canónica. |
| `ejecutar_verificacion_fase2.py` / `fase3.py` / `fase4.py` | Baterías de verificación del algoritmo Worm (energía, densidad, permutaciones). |
| `ejecutar_verificacion_obdm.py` | Verificación de la matriz densidad de un cuerpo (OBDM). |
| `generar_densidad_comparativa.py` | Gráficos comparativos de densidad radial (con/sin *swap*, bandas $3\sigma$). |
| `generar_graficos_energia.py` | Gráficos de energía y tablas de verificación. |
| `generar_graficos_obdm.py` | Mapas de calor 2D y cortes de la OBDM. |
| `generar_graficos_superfluidez.py` | Análisis y figuras de fracción superfluida. |
| `generar_graficos_dornheim_worm.py` | Figuras de *snapshots* del algoritmo Worm. |
| `worm_config.py` | Utilidades compartidas de configuración para los scripts Worm. |

---

## 6. Referencias Analíticas (Oscilador Armónico 1D)

En unidades reducidas ($\hbar = m = \omega = k_B = 1$) con $\beta = 1/T$ y $\alpha = \tanh(\beta/2)$:

- **Densidad:** $n(x) = \sqrt{\frac{\alpha}{\pi}} e^{-\alpha x^2}$
- **Dispersión cuadrática media:** $\langle x^2 \rangle = \frac{1}{2\alpha}$
- **Energías (virial):** $E_{\text{pot}} = E_{\text{kin}} = \frac{1}{2}\langle x^2 \rangle = \frac{1}{4}\coth\left(\frac{1}{2T}\right)$, $\quad E_{\text{tot}} = 2E_{\text{pot}}$
- **Límites de $E_{\text{pot}}$:** $0.25$ cuando $T \to 0$ (estado fundamental); $\frac{1}{2}T$ cuando $T \to \infty$ (equipartición clásica).
