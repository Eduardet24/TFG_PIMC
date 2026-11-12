#!/usr/bin/env python3
r"""
generar_graficos_obdm.py

Generación y análisis de la Matriz de Densidad de Una Partícula (OBDM)
rho_1(r, r') = <Psi^\dagger(r) Psi(r')> medida con las cadenas abiertas
del Algoritmo Worm (Prokof'ev-Svistunov / Dornheim 2014 / Boninsegni 2006).

Incluye:
  1. Solución analítica exacta de referencia (Kernel de Mehler y expansión en ciclos).
  2. Carga y normalización de la matriz medida en Sector G (gap == 1).
  3. Gráficas de mapa de calor 2D (PIMC vs Analítico vs Residuos).
  4. Cortes unidimensionales: diagonal local rho_1(x, x) y cortes de coherencia rho_1(x, -x), rho_1(x, 0).
  5. Decaimiento de la coherencia cuántica en función de la distancia de separación s = |r - r'|.
  6. Exportación de figuras de alta resolución a Overleaf.
"""

import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from worm_config import get_n_particles

OUTPUT_DIR = os.path.join(SCRIPT_DIR, "Resultados", "Verificacion_OBDM")
OVERLEAF_IMG_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "overleaf", "imágenes")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(OVERLEAF_IMG_DIR, exist_ok=True)

# Parámetros físicos por defecto
N_PARTICLES = get_n_particles(os.path.join(SCRIPT_DIR, "MC.cfg"), "WORM_FASE4_PARTICLES")
if N_PARTICLES < 2 and "WORM_PARTICLES" not in os.environ and "WORM_FASE4_PARTICLES" not in os.environ:
    N_PARTICLES = 2
TEMP = 0.5
BETA = 1.0 / TEMP
OMEGA = 1.0


# ==============================================================================
# 1. Modelos Analíticos Exactos
# ==============================================================================

def canonical_bose_partition(n_particles, temp, omega=OMEGA):
    """Calcula la función de partición canónica Z_N recursivamente."""
    beta = 1.0 / temp
    partition = np.zeros(n_particles + 1)
    partition[0] = 1.0
    for count in range(1, n_particles + 1):
        for cycle in range(1, count + 1):
            q = np.exp(-cycle * beta * omega)
            one_cycle_partition = 1.0 / (1.0 - q)**2
            partition[count] += one_cycle_partition * partition[count - cycle] / count
    return partition


def normalized_mehler_1d(x, xp, tau, omega=OMEGA):
    """
    Kernel de Mehler unidimensional exacto normalizado a integral sobre x = 1:
    rho_HO(x, xp; tau) / z_1D(tau).
    """
    u = tau * omega
    sinh_u = np.sinh(u)
    cosh_u = np.cosh(u)
    tanh_half = np.tanh(0.5 * u)
    pref = np.sqrt(omega * tanh_half / np.pi)
    exponent = -(omega / (2.0 * sinh_u)) * ((x**2 + xp**2) * cosh_u - 2.0 * x * xp)
    return pref * np.exp(exponent)


def analytical_bose_obdm_1d(x, xp, n_particles, temp, omega=OMEGA, per_particle=False):
    r"""
    Matriz de densidad marginal 1D rho_1(x, x') = \int dy rho_1(x, y; x', y)
    para N bosones no interactuantes en una trampa armónica isótropa 2D.
    """
    z_can = canonical_bose_partition(n_particles, temp, omega)
    beta = 1.0 / temp
    
    # Manejar tanto arrays 2D para matrices como 1D para cortes
    if isinstance(x, np.ndarray) and isinstance(xp, np.ndarray) and x.ndim == 2 and xp.ndim == 2:
        val = np.zeros_like(x, dtype=float)
    elif isinstance(x, np.ndarray) and isinstance(xp, np.ndarray) and x.ndim == 1 and xp.ndim == 1 and x.shape == xp.shape:
        val = np.zeros_like(x, dtype=float)
    else:
        val = np.zeros_like(x, dtype=float)

    for k in range(1, n_particles + 1):
        q = np.exp(-k * beta * omega)
        z2d = 1.0 / (1.0 - q)**2
        # Peso del ciclo k en la matriz de densidad de una partícula
        weight = z2d * z_can[n_particles - k] / z_can[n_particles]
        if per_particle:
            weight /= n_particles
        val += weight * normalized_mehler_1d(x, xp, k * beta, omega)
        
    return val


def analytical_distinguishable_obdm_1d(x, xp, n_particles, temp, omega=OMEGA, per_particle=False):
    """Matriz de densidad marginal 1D para partículas distinguibles (sin intercambio, k=1)."""
    beta = 1.0 / temp
    norm = 1.0 if per_particle else float(n_particles)
    return norm * normalized_mehler_1d(x, xp, beta, omega)


# ==============================================================================
# 2. Carga y Normalización de Datos PIMC
# ==============================================================================

def load_worm_obdm_data(data_path):
    """
    Lee el archivo outworm_obdm.dat generado por el motor en C.
    Devuelve x_grid, matrix_2d, metadata_dict.
    """
    if not os.path.isfile(data_path):
        raise FileNotFoundError(f"No se encontró el fichero de datos: {data_path}")

    metadata = {}
    with open(data_path, "r") as f:
        for line in f:
            if line.startswith("# grid:"):
                # Formato: # grid: 100 min: -3.5 max: 3.5 step: 0.07 samples: 1234 norm_trace: ... N: 2
                parts = line.strip("# \n").split()
                for idx in range(0, len(parts)-1, 2):
                    key = parts[idx].rstrip(":")
                    val = parts[idx+1]
                    try:
                        metadata[key] = int(val) if val.isdigit() or val.startswith("-") and val[1:].isdigit() else float(val)
                    except ValueError:
                        metadata[key] = val
                break

    grid_size = int(metadata.get("grid", 100))
    min_x = float(metadata.get("min", -3.5))
    max_x = float(metadata.get("max", 3.5))
    step = float(metadata.get("step", (max_x - min_x) / grid_size))

    # Cargar datos tabulares: col 0: x, col 1: x_prime, col 2: count, col 3: rho1_norm
    raw_data = np.loadtxt(data_path, comments="#")
    x_vals = np.linspace(min_x + 0.5 * step, max_x - 0.5 * step, grid_size)
    
    # Reconstruir matriz cuadrada
    matrix_counts = np.zeros((grid_size, grid_size), dtype=float)
    if raw_data.shape[1] >= 3:
        for row in raw_data:
            ix = int(np.floor((row[0] - min_x) / step))
            ixp = int(np.floor((row[1] - min_x) / step))
            if 0 <= ix < grid_size and 0 <= ixp < grid_size:
                matrix_counts[ix, ixp] = row[2]
                matrix_counts[ixp, ix] = row[2]  # Simetrización Hermítica

    # Normalización exacta por traza: \int rho_1(x, x) dx = N
    trace_integral = np.sum(np.diag(matrix_counts)) * step
    n_part = int(metadata.get("N", N_PARTICLES))
    if trace_integral > 0:
        matrix_norm = matrix_counts * (float(n_part) / trace_integral)
    else:
        matrix_norm = matrix_counts

    metadata["trace_integral"] = trace_integral
    return x_vals, matrix_norm, matrix_counts, metadata


def load_worm_dist_data(dist_path):
    """Lee el archivo outworm_obdm_dist.dat."""
    if not os.path.isfile(dist_path):
        return None, None
    raw = np.loadtxt(dist_path, comments="#")
    if len(raw) == 0:
        return None, None
    s_vals = raw[:, 0]
    p_s = raw[:, 2] if raw.shape[1] >= 3 else raw[:, 1]
    return s_vals, p_s


# ==============================================================================
# 3. Generación de Gráficos de Publicación
# ==============================================================================

def generar_mapa_calor_obdm(x_vals, pimc_matrix, n_particles=N_PARTICLES, temp=TEMP, save_dir=OUTPUT_DIR):
    """
    Genera el panel comparativo de Mapas de Calor 2D:
      [1] PIMC Worm (Muestreado con cadenas abiertas)
      [2] Analítico Bose (Solución canónica exacta)
      [3] Residuos / Discrepancia absoluta
    """
    X, XP = np.meshgrid(x_vals, x_vals)
    bose_exact = analytical_bose_obdm_1d(X, XP, n_particles, temp)
    diff_matrix = pimc_matrix - bose_exact

    vmax = max(np.max(pimc_matrix), np.max(bose_exact))
    vmin = 0.0

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)

    # Panel 1: PIMC
    im0 = axes[0].imshow(pimc_matrix, extent=[x_vals[0], x_vals[-1], x_vals[0], x_vals[-1]],
                         origin="lower", cmap="viridis", vmin=vmin, vmax=vmax)
    axes[0].set_title(r"$\rho_1(x, x')$ -- PIMC Algoritmo Worm", fontsize=12, fontweight="bold")
    axes[0].set_xlabel(r"Posición $x$ [$a_{\mathrm{ho}}$]", fontsize=11)
    axes[0].set_ylabel(r"Posición $x'$ [$a_{\mathrm{ho}}$]", fontsize=11)
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    # Panel 2: Analítico Bose
    im1 = axes[1].imshow(bose_exact, extent=[x_vals[0], x_vals[-1], x_vals[0], x_vals[-1]],
                         origin="lower", cmap="viridis", vmin=vmin, vmax=vmax)
    axes[1].set_title(r"$\rho_1(x, x')$ -- Analítico Bose Canónico", fontsize=12, fontweight="bold")
    axes[1].set_xlabel(r"Posición $x$ [$a_{\mathrm{ho}}$]", fontsize=11)
    axes[1].set_ylabel(r"Posición $x'$ [$a_{\mathrm{ho}}$]", fontsize=11)
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    # Panel 3: Residuos
    max_err = max(np.max(np.abs(diff_matrix)), 1e-4)
    im2 = axes[2].imshow(diff_matrix, extent=[x_vals[0], x_vals[-1], x_vals[0], x_vals[-1]],
                         origin="lower", cmap="coolwarm", vmin=-max_err, vmax=max_err)
    axes[2].set_title(r"Residuo: $\rho_1^{\mathrm{PIMC}} - \rho_1^{\mathrm{Bose}}$", fontsize=12, fontweight="bold")
    axes[2].set_xlabel(r"Posición $x$ [$a_{\mathrm{ho}}$]", fontsize=11)
    axes[2].set_ylabel(r"Posición $x'$ [$a_{\mathrm{ho}}$]", fontsize=11)
    fig.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)

    out_file = os.path.join(save_dir, "01_obdm_mapa_calor_2d_comparativo.png")
    fig.savefig(out_file, dpi=300)
    plt.close(fig)
    print(f"  -> Guardado: {out_file}")

    overleaf_file = os.path.join(OVERLEAF_IMG_DIR, "01_obdm_mapa_calor_2d_comparativo.png")
    if os.path.exists(OVERLEAF_IMG_DIR):
        import shutil
        shutil.copy2(out_file, overleaf_file)


def generar_cortes_obdm(x_vals, pimc_matrix, n_particles=N_PARTICLES, temp=TEMP, save_dir=OUTPUT_DIR):
    """
    Genera cortes 1D de la matriz de densidad:
      Panel A: Diagonal local rho_1(x, x) (perfil de densidad espacial).
      Panel B: Corte anti-diagonal rho_1(x, -x) (coherencia cuántica simétrica).
      Panel C: Coherencia espacial relativa rho_1(x, 0) respecto al origen.
    """
    diag_pimc = np.diag(pimc_matrix)
    diag_bose = analytical_bose_obdm_1d(x_vals, x_vals, n_particles, temp)
    diag_dist = analytical_distinguishable_obdm_1d(x_vals, x_vals, n_particles, temp)

    # Corte anti-diagonal: x' = -x
    # Dado que x_vals va de -L a +L, el elemento x_vals[i] tiene opuesto x_vals[-(i+1)]
    antidiag_pimc = np.array([pimc_matrix[i, -(i+1)] for i in range(len(x_vals))])
    antidiag_bose = analytical_bose_obdm_1d(x_vals, -x_vals, n_particles, temp)
    antidiag_dist = analytical_distinguishable_obdm_1d(x_vals, -x_vals, n_particles, temp)

    # Corte x' = 0 (coherencia con el centro de la trampa)
    idx_center = np.argmin(np.abs(x_vals))
    cut_center_pimc = pimc_matrix[:, idx_center]
    cut_center_bose = analytical_bose_obdm_1d(x_vals, np.zeros_like(x_vals), n_particles, temp)
    cut_center_dist = analytical_distinguishable_obdm_1d(x_vals, np.zeros_like(x_vals), n_particles, temp)

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)

    # Panel A: Diagonal Local n(x)
    axes[0].plot(x_vals, diag_bose, color="crimson", linewidth=2.2, label="Analítico Bose Canónico", zorder=2)
    axes[0].plot(x_vals, diag_dist, color="black", linestyle="--", linewidth=1.8, label="Distinguibles", zorder=1)
    axes[0].plot(x_vals, diag_pimc, "o", color="#0055aa", markersize=4.0, alpha=0.85, label="PIMC Worm Sector G", zorder=3)
    axes[0].set_title(r"Diagonal Local: $\rho_1(x, x) = n(x)$", fontsize=11, fontweight="bold")
    axes[0].set_xlabel(r"Posición $x$ [$a_{\mathrm{ho}}$]", fontsize=10)
    axes[0].set_ylabel(r"Densidad Marginal $\rho_1(x, x)$", fontsize=10)
    axes[0].legend(fontsize=9, loc="upper right")
    axes[0].grid(True, linestyle=":", alpha=0.5)

    # Panel B: Coherencia Anti-diagonal rho_1(x, -x)
    axes[1].plot(x_vals, antidiag_bose, color="crimson", linewidth=2.2, label="Analítico Bose", zorder=2)
    axes[1].plot(x_vals, antidiag_dist, color="black", linestyle="--", linewidth=1.8, label="Distinguibles", zorder=1)
    axes[1].plot(x_vals, antidiag_pimc, "s", color="darkorange", markersize=4.0, alpha=0.85, label="PIMC Worm", zorder=3)
    axes[1].set_title(r"Coherencia Anti-diagonal: $\rho_1(x, -x)$", fontsize=11, fontweight="bold")
    axes[1].set_xlabel(r"Posición $x$ [$a_{\mathrm{ho}}$]", fontsize=10)
    axes[1].set_ylabel(r"Coherencia Cuántica $\rho_1(x, -x)$", fontsize=10)
    axes[1].legend(fontsize=9, loc="upper right")
    axes[1].grid(True, linestyle=":", alpha=0.5)

    # Panel C: Coherencia con el Centro rho_1(x, 0)
    axes[2].plot(x_vals, cut_center_bose, color="crimson", linewidth=2.2, label="Analítico Bose", zorder=2)
    axes[2].plot(x_vals, cut_center_dist, color="black", linestyle="--", linewidth=1.8, label="Distinguibles", zorder=1)
    axes[2].plot(x_vals, cut_center_pimc, "^", color="purple", markersize=4.0, alpha=0.85, label="PIMC Worm", zorder=3)
    axes[2].set_title(r"Coherencia con el Centro: $\rho_1(x, 0)$", fontsize=11, fontweight="bold")
    axes[2].set_xlabel(r"Posición $x$ [$a_{\mathrm{ho}}$]", fontsize=10)
    axes[2].set_ylabel(r"Amplitud de Coherencia $\rho_1(x, 0)$", fontsize=10)
    axes[2].legend(fontsize=9, loc="upper right")
    axes[2].grid(True, linestyle=":", alpha=0.5)

    out_file = os.path.join(save_dir, "02_obdm_cortes_coherencia_cuantica.png")
    fig.savefig(out_file, dpi=300)
    plt.close(fig)
    print(f"  -> Guardado: {out_file}")

    overleaf_file = os.path.join(OVERLEAF_IMG_DIR, "02_obdm_cortes_coherencia_cuantica.png")
    if os.path.exists(OVERLEAF_IMG_DIR):
        import shutil
        shutil.copy2(out_file, overleaf_file)


def generar_grafico_distancia_extremos(dist_path, n_particles=N_PARTICLES, temp=TEMP, save_dir=OUTPUT_DIR):
    """Grafica la distribución de distancia espacial entre los extremos Ira y Masha."""
    s_vals, p_s = load_worm_dist_data(dist_path)
    if s_vals is None:
        print("  [Aviso] No se encontró outworm_obdm_dist.dat, saltando gráfica de distancia.")
        return

    fig, ax = plt.subplots(figsize=(7.5, 4.8), constrained_layout=True)
    mask = s_vals <= 5.0
    ax.plot(s_vals[mask], p_s[mask], "o-", color="teal", markersize=3.5, linewidth=1.5,
            label=r"PIMC: $p(s)$ extremos worm (gap == 1)")
    ax.set_title(r"Distribución de Separación de Extremos del Worm $s = |\mathbf{r}_{\mathrm{ira}} - \mathbf{r}_{\mathrm{masha}}|$",
                 fontsize=11, fontweight="bold")
    ax.set_xlabel(r"Distancia de separación $s$ [$a_{\mathrm{ho}}$]", fontsize=10)
    ax.set_ylabel(r"Densidad de probabilidad $p(s)$", fontsize=10)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, linestyle=":", alpha=0.5)

    out_file = os.path.join(save_dir, "03_obdm_distribucion_distancia_extremos.png")
    fig.savefig(out_file, dpi=300)
    plt.close(fig)
    print(f"  -> Guardado: {out_file}")

    overleaf_file = os.path.join(OVERLEAF_IMG_DIR, "03_obdm_distribucion_distancia_extremos.png")
    if os.path.exists(OVERLEAF_IMG_DIR):
        import shutil
        shutil.copy2(out_file, overleaf_file)


# ==============================================================================
# 4. Función Principal de Análisis
# ==============================================================================

def main():
    data_dir = sys.argv[1] if len(sys.argv) > 1 else SCRIPT_DIR
    obdm_file = os.path.join(data_dir, "outworm_obdm.dat")
    dist_file = os.path.join(data_dir, "outworm_obdm_dist.dat")

    if not os.path.isfile(obdm_file):
        print(f"[Aviso] No se encontró '{obdm_file}'. Ejecute primero una simulación PIMC con Worm.")
        return 1

    print(f"\n==================================================================")
    print(f"  ANÁLISIS DE LA MATRIZ DE DENSIDAD DE UNA PARTÍCULA (OBDM)")
    print(f"==================================================================")
    print(f"  Leyendo archivo: {obdm_file}")

    x_vals, pimc_norm, pimc_counts, meta = load_worm_obdm_data(obdm_file)
    n_part = meta.get("N", N_PARTICLES)
    samples = meta.get("samples", 0)
    step = meta.get("step", 0.0)

    print(f"  -> Puntos de cuadrícula: {len(x_vals)}")
    print(f"  -> Muestras gap == 1  : {samples}")
    print(f"  -> Rango espacial     : [{x_vals[0]:.2f}, {x_vals[-1]:.2f}] (paso dx = {step:.4f})")
    print(f"  -> Integral traza raw : {meta.get('trace_integral', 0.0):.2f}")

    # 1. Mapa de Calor 2D
    generar_mapa_calor_obdm(x_vals, pimc_norm, n_particles=n_part, temp=TEMP, save_dir=OUTPUT_DIR)

    # 2. Cortes Unidimensionales
    generar_cortes_obdm(x_vals, pimc_norm, n_particles=n_part, temp=TEMP, save_dir=OUTPUT_DIR)

    # 3. Distribución de Distancia de Extremos
    generar_grafico_distancia_extremos(dist_file, n_particles=n_part, temp=TEMP, save_dir=OUTPUT_DIR)

    print(f"==================================================================\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
