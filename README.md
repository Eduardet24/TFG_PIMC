# Simulación PIMC de Osciladores Armónicos 1D y 2D

Este repositorio contiene la implementación y el entorno de trabajo para simulaciones de **Monte Carlo Cuántico** (Quantum Monte Carlo, QMC), principalmente **Path Integral Monte Carlo (PIMC)** para osciladores armónicos 1D y 2D. Incluye un benchmark 2D de partículas ideales con `N > 1`; el muestreador actual conserva caminos etiquetados y aún no incluye permutaciones bosónicas. El código base en C está desarrollado originalmente por G. E. Astrakharchik (2001–2017) y ha sido configurado y complementado para el Trabajo de Fin de Grado (TFG).

---

## Índice

- [1. Estructura del Repositorio](#1-estructura-del-repositorio)
- [2. Requisitos y Dependencias](#2-requisitos-y-dependencias)
- [3. Compilación](#3-compilación)
- [4. Configuración (`MC.cfg`)](#4-configuración-mccfg)
- [5. Flujo de Simulación](#5-flujo-de-simulación)
  - [5.1 Simulación Individual](#51-simulación-individual)
  - [5.2 Campañas Automatizadas por Lotes (`Script_Lote_PIMC.py`)](#52-campañas-automatizadas-por-lotes-script_lote_pimcpy)
- [6. Scripts de Análisis y Gráficas](#6-scripts-de-análisis-y-gráficas)
  - [6.1 Densidad y Virial (`Script_Armonico1D_basico.py`)](#61-densidad-y-virial-script_armonico1d_basicopy)
  - [6.2 Convergencia en Beads (`Script_EPvsBeads`)](#62-convergencia-en-beads-script_epvsbeads)
  - [6.3 Energía Potencial vs Temperatura (`Script_EPvsT`)](#63-energía-potencial-vs-temperatura-script_epvst)
  - [6.4 Ajustes Cuadráticos (`Script_Convergencia_AR0.py`)](#64-ajustes-cuadráticos-script_convergencia_ar0py)
  - [6.5 Script Gnuplot (`scritpt_Rdz_grigori`)](#65-script-gnuplot-scritpt_rdz_grigori)
  - [6.6 Benchmark de N partículas en 2D](#66-benchmark-de-n-partículas-en-2d)
- [7. Organización de Resultados (`Resultados/`)](#7-organización-de-resultados-resultados)
- [8. Modelos Físicos y Referencias Analíticas](#8-modelos-físicos-y-referencias-analíticas)

---

## 1. Estructura del Repositorio

```text
Grigori_codigo/
├── MC                  # Ejecutable binario compilado
├── MC.cfg              # Archivo principal de configuración y parámetros
├── main.c, pimc.c, ... # Código fuente en C (módulos MC, PIMC, DMC, VMC, etc.)
├── Makefile.txt        # Makefile para compilación con Intel (icc)
├── Makefile.AMD        # Makefile alternativo
│
├── Script_Lote_PIMC.py          # Automatización de barridos (T, N_beads)
├── Script_Armonico1D_basico.py  # Análisis de densidad n(x) y virial para una corrida
├── Script_EPvsBeads             # Estudio de convergencia Epot vs 1/N_beads
├── Script_EPvsT                 # Comparación Epot vs T (analítico, clásico y PIMC)
├── Script_Convergencia_AR0.py   # Extrapolaciones parabólicas de convergencia
├── Script_Comparacion_N_Particulas_2D.py # Benchmark N>1 y referencias ideales 2D
├── scritpt_Rdz_grigori          # Script de graficación rápida con gnuplot
│
├── Resultados/
│   └── 1DArmonico/
│       ├── AR0/                 # Series con acceptance_rate = 0 (paso fijo dt_vmc)
│       ├── AR50/                # Series con acceptance_rate = 50%
│       │   └── lote_modos_normales_blck1000/
│       │       └── Temp<T>/
│       │           └── N<beads>/
│       └── ejecucion_actual/    # Salidas temporales del análisis individual
│
└── Doxyfile            # Configuración para documentación con Doxygen
```

---

## 2. Requisitos y Dependencias

- **Compilador de C:** `gcc` con soporte OpenMP o Intel C Compiler (`icc`).
- **Python 3:** versión 3.8 o superior.
- **Librerías Python:**
  ```bash
  pip install numpy matplotlib scipy
  ```
- **Opcionales:**
  - `gnuplot` (para ejecutar `scritpt_Rdz_grigori`).
  - `doxygen` (para generar documentación del código C).

---

## 3. Compilación

### Compilación estándar con GCC (recomendada)

Desde el directorio `Grigori_codigo/`:

```bash
gcc -O2 -fopenmp *.c -lm -o MC
```

### Compilación con Intel C Compiler (`icc`)

```bash
make -f Makefile.txt
```

Para limpiar archivos objeto generados:

```bash
make -f Makefile.txt clean
```

---

## 4. Configuración (`MC.cfg`)

El archivo [MC.cfg](file:///home/edugimeno/Documents/Uni/TFG/Grigori_codigo/MC.cfg) gobierna todos los parámetros de la simulación. 

> [!IMPORTANT]
> **Sintaxis estricta:** Los parámetros deben especificarse sin espacio antes del signo igual (ejemplo: `T= 0.5`, **no** `T = 0.5`).

### Parámetros principales

| Parámetro | Valor típico | Descripción |
| :--- | :--- | :--- |
| `MC=` | `3` | Tipo de simulación: `0` (VMC), `1` (DMC), `3` (PIMC). |
| `Npop=` | `5` .. `30` | Número de réplicas temporales o *beads* ($N_{\text{beads}}$). |
| `T=` | `0.5` .. `4.0` | Temperatura del sistema en unidades reducidas ($\hbar = m = \omega = k_B = 1$). |
| `N=` | `1` | Número de partículas (1 para oscilador armónico individual). |
| `omega_z=` | `1` | Frecuencia del potencial armónico $V(z) = \frac{1}{2} m \omega_z^2 z^2$. |
| `acceptance_rate=` | `50` | Tasa de aceptación objetivo (%). Si es `0`, se usa desplazamiento fijo vía `dt_vmc`. |
| `dt_vmc=` | `1e-2` | Paso temporal/desplazamiento cuando `acceptance_rate= 0`. |
| `blck_heating=` | `100` | Bloques de termalización/calentamiento descartados. |
| `blck=` | `1000` | Número de bloques para toma de medidas. |
| `Niter=` | `1000` | Iteraciones por bloque. |
| `Nmeasure=` | `1` | Frecuencia de muestreo (cada cuántas iteraciones se mide). |
| `measure_energy=` | `1` | Medir energías y guardar en `file_energy` (`oute.dat`). |
| `measure_RadDistr=`| `1` | Medir distribución radial/espacial en `file_RD` (`outrd.dat`). |

### Archivos de salida producidos por `./MC`

- **`oute.dat`**: Registro de energía por bloque/muestra. Columnas: `# Etot Epot Ekin`.
- **`outrd.dat`**: Distribución de posiciones/radial $n(r)$.
- **`MC.log`**: Registro de ejecución, estadísticas de aceptación y resumen.

---

## 5. Flujo de Simulación

### 5.1 Simulación Individual

1. Ajusta los parámetros deseados en `MC.cfg` (por ejemplo, `Npop= 20`, `T= 1.0`).
2. Ejecuta el binario:
   ```bash
   ./MC
   ```
3. Analiza los resultados inmediatamente con:
   ```bash
   python3 Script_Armonico1D_basico.py
   ```
   Las gráficas y el resumen numérico se guardarán en `Resultados/1DArmonico/ejecucion_actual/`.

### 5.2 Campañas Automatizadas por Lotes (`Script_Lote_PIMC.py`)

Para barridos sistemáticos en temperatura y número de beads sin intervención manual:

```bash
# Ver el plan de ejecución sin correr simulaciones (dry-run)
python3 Script_Lote_PIMC.py

# Lanzar la matriz predeterminada
python3 Script_Lote_PIMC.py --run

# Lanzar una campaña personalizada
python3 Script_Lote_PIMC.py --run --temperaturas 0.5,1.0,2.0 --beads 5,10,15,20,25,30 --campaign lote_modos_normales_blck1000
```

**Características del script de lotes:**
- Crea la jerarquía `Resultados/1DArmonico/AR<rate>/<campaña>/Temp<T>/N<beads>/`.
- Copia de forma aislada la configuración `MC.cfg` y el ejecutable `MC` en cada carpeta.
- Guarda el log completo en `run.log`.
- Comprueba la validez y número de filas de `oute.dat` antes de marcar la corrida con `.complete`.
- Si una corrida ya fue completada exitosamente, la omite automáticamente.

### 5.3 Benchmark de partículas ideales en 2D

El benchmark usa `INTERACTION_ABSENT`, una trampa armónica 2D y el muestreador exacto de modos normales para cada camino. Se ejecuta con `N > 1`, pero sin intercambios: representa partículas distinguibles. Las curvas analíticas distinguible y Bose canónica permiten comprobar la energía y el perfil de densidad antes de implementar el Worm Algorithm.

```bash
python3 Script_Comparacion_N_Particulas_2D.py --run --particles 5
```

Con `--run`, el script compila los fuentes C actuales en un directorio temporal y usa ese ejecutable para las simulaciones, evitando depender del `MC` precompilado del repositorio. Los datos se guardan bajo `Resultados/NParticulas2D/`. Para evaluar convergencia de discretización, repetir la corrida con distintos valores de `--beads`; las referencias analíticas corresponden al límite continuo.

---

## 6. Scripts de Análisis y Gráficas

### 6.1 Densidad y Virial (`Script_Armonico1D_basico.py`)
- **Propósito:** Procesa `outrd.dat` y `oute.dat` de la simulación actual o local.
- **Salida:**
  - Compara la distribución $n(x)$ simulada con la curva analítica cuántica $n_{\text{exact}}(x) = \sqrt{\frac{\alpha}{\pi}} e^{-\alpha x^2}$.
  - Gráfica de residuos con barras de error estadístico.
  - Estimación de $\langle x^2 \rangle$, energía potencial $E_{\text{pot}} = \frac{1}{2} \langle x^2 \rangle$, energía cinética y energía total según el teorema del virial.
  - Guarda los archivos en `Resultados/1DArmonico/ejecucion_actual/`.

### 6.2 Convergencia en Beads (`Script_EPvsBeads`)
- **Propósito:** Evalúa la convergencia de la energía potencial en función del número de cuentas temporales ($N_{\text{beads}}$) a temperatura fija.
- **Método:**
  - Lee los `oute.dat` de todas las carpetas `N<beads>` completadas para la temperatura configurada en `MC.cfg`.
  - Promedia por bloques para estimar el valor medio y la incertidumbre estadística de $E_{\text{pot}}$.
  - Grafica $E_{\text{pot}}$ frente a $1/N_{\text{beads}}$ y ajusta una recta para extrapolar al límite continuo $N_{\text{beads}} \to \infty$ ($1/N \to 0$).
  - Compara el resultado extrapolado con el valor cuántico analítico exacto.
- **Ejecución:**
  ```bash
  python3 Script_EPvsBeads
  ```
  Guarda la gráfica en `convergencia_energia_potencialvsBeads.png` dentro del directorio de la campaña.

### 6.3 Energía Potencial vs Temperatura (`Script_EPvsT`)
- **Propósito:** Analiza la respuesta térmica del oscilador armónico comparando los regímenes cuántico y clásico.
- **Método:**
  - Escanea los resultados de una campaña, detecta los $N_{\text{beads}}$ disponibles y solicita interactivamente cuál analizar.
  - Representa $\langle U \rangle(T)$ frente a:
    1. **Modelo clásico:** $E_{\text{pot}} = \frac{1}{2} T$ (teorema de equipartición).
    2. **Modelo cuántico exacto:** $E_{\text{pot}} = \frac{1}{4} \coth\left(\frac{1}{2T}\right)$.
    3. **Puntos PIMC:** Con sus respectivas barras de error.
- **Ejecución:**
  ```bash
  python3 Script_EPvsT
  ```
  Genera `energia_potencial_vs_temperatura.png`.

### 6.4 Ajustes Cuadráticos (`Script_Convergencia_AR0.py`)
- Realiza ajustes polinómicos de orden 2 en $1/N_{\text{beads}}$ para series con tasa de aceptación cero (`AR0`) o fija (`AR50`) a partir de archivos `convergencia_beads_1D.csv`.

### 6.5 Script Gnuplot (`scritpt_Rdz_grigori`)
- Script rápido para visualizar la función de distribución radial usando gnuplot:
  ```bash
  gnuplot scritpt_Rdz_grigori
  ```

### 6.6 Benchmark de N partículas en 2D

`Script_Comparacion_N_Particulas_2D.py` compara energía por partícula y densidad areal radial con las soluciones ideales distinguible y Bose canónica a N fijo. Promedia las últimas mitades de las medidas y estima incertidumbres agrupando las energías por bloques y la densidad por bloques guardados. La rama PIMC actual no muestrea permutaciones, así que el resultado numérico debe contrastarse con la curva distinguible; la curva Bose es el objetivo de la futura implementación Worm, no una predicción de la simulación actual.

---

## 7. Organización de Resultados (`Resultados/`)

Los resultados se archivan siguiendo una taxonomía clara según la tasa de aceptación, la campaña y los parámetros físicos:

```text
Resultados/
└── 1DArmonico/
    ├── AR0/                                   # Paso fijo dt_vmc
    │   ├── N10/
    │   ├── N20/
    │   └── N30/
    └── AR50/                                  # Tasa de aceptación 50%
        └── lote_modos_normales_blck1000/      # Nombre de la campaña
            ├── Temp0.5/
            │   ├── N5/
            │   │   ├── MC.cfg                 # Parámetros exactos de esta corrida
            │   │   ├── oute.dat               # Energías por bloque
            │   │   ├── outrd.dat              # Distribución radial
            │   │   ├── run.log                # Salida estándar del ejecutable
            │   │   └── .complete              # Marcador de corrida finalizada
            │   ├── N10/
            │   └── ...
            ├── Temp1/
            ├── Temp2/
            ├── convergencia_energia_potencialvsBeads.png
            └── energia_potencial_vs_temperatura.png
```

---

## 8. Modelos Físicos y Referencias Analíticas

Para una partícula en un potencial armónico $V(x) = \frac{1}{2} m \omega^2 x^2$ en unidades reducidas ($\hbar = m = \omega = k_B = 1$):

### 1. Parámetro térmico
$$\beta = \frac{1}{T}, \quad \alpha = \tanh\left(\frac{\beta}{2}\right) = \tanh\left(\frac{1}{2T}\right)$$

### 2. Distribución de densidad analítica
$$n(x) = \left(\frac{\alpha}{\pi}\right)^{1/2} \exp\left(-\alpha x^2\right)$$

### 3. Dispersión cuadrática media
$$\langle x^2 \rangle = \frac{1}{2\alpha} = \frac{1}{2\tanh\left(\frac{1}{2T}\right)}$$

### 4. Energías por el Teorema del Virial
Para el oscilador armónico unidimensional:
$$E_{\text{pot}} = \frac{1}{2} \langle x^2 \rangle, \quad E_{\text{kin}} = E_{\text{pot}}, \quad E_{\text{tot}} = 2 E_{\text{pot}} = \frac{1}{2\tanh\left(\frac{1}{2T}\right)}$$

### 5. Límites asintóticos de la energía potencial
- **Límite de baja temperatura ($T \to 0$ / fundamental cuántico):**
  $$E_{\text{pot}} \to \frac{1}{4} \hbar \omega = 0.25$$
- **Límite de alta temperatura ($T \to \infty$ / límite clásico):**
  $$\tanh(x) \approx x \implies E_{\text{pot}} \to \frac{1}{2} T$$
