#!/usr/bin/env python3
r"""
generar_graficos_obdm.py

Generación y análisis de la Matriz de Densidad de Una Partícula (OBDM)
rho_1(r, r') = <Psi^\dagger(r) Psi(r')> medida con las cadenas abiertas
del Algoritmo Worm (Prokof'ev-Svistunov / Dornheim 2014 / Boninsegni 2006).

Incluye:
  1. Solución analítica exacta de referencia (Kernel de Mehler y expansión en ciclos).
    2. Carga de la matriz medida en Sector G (gap == 1), normalizada por la traza en C.
  3. Gráficas de mapa de calor 2D (PIMC vs Analítico vs Residuos).
  4. Cortes unidimensionales: diagonal local rho_1(x, x) y cortes de coherencia rho_1(x, -x), rho_1(x, 0).
  5. Decaimiento de la coherencia cuántica en función de la distancia de separación s = |r - r'|.
  6. Exportación de figuras de alta resolución a Overleaf.
"""

import argparse
import os
import shutil
import sys
import numpy as np
import matplotlib.pyplot as plt

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)
SCRIPT_DIR = os.path.dirname(_THIS_DIR)

from worm_config import get_n_particles
from pimc_common import canonical_bose_partition as _full_bose_partition, read_cfg_value

OUTPUT_DIR = os.path.join(SCRIPT_DIR, "Resultados", "Verificacion_OBDM")
OVERLEAF_IMG_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "overleaf", "imágenes")


def _get_cfg_number(cfg_path, key, default=None):
    if os.path.isfile(cfg_path):
        with open(cfg_path, "r", encoding="utf-8", errors="ignore") as f:
            try:
                return read_cfg_value(f.read(), key)
            except ValueError:
                return default
    return default


# Parámetros físicos por defecto (usados si el caller no los pasa ni están en metadatos)
N_PARTICLES = get_n_particles(os.path.join(SCRIPT_DIR, "MC.cfg"), "WORM_FASE4_PARTICLES")
if N_PARTICLES < 2 and "WORM_PARTICLES" not in os.environ and "WORM_FASE4_PARTICLES" not in os.environ:
    N_PARTICLES = 2
TEMP = _get_cfg_number(os.path.join(SCRIPT_DIR, "MC.cfg"), "T", 0.5)
OMEGA = _get_cfg_number(os.path.join(SCRIPT_DIR, "MC.cfg"), "omega_x", 1.0)


# ==============================================================================
# 1. Modelos Analíticos Exactos
# ==============================================================================

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
    z_can = _full_bose_partition(n_particles, temp, omega)[0]
    beta = 1.0 / temp

    x_arr, xp_arr = np.broadcast_arrays(
        np.asarray(x, dtype=float), np.asarray(xp, dtype=float)
    )
    val = np.zeros_like(x_arr, dtype=float)

    for k in range(1, n_particles + 1):
        q = np.exp(-k * beta * omega)
        z2d = 1.0 / (1.0 - q)**2
        # Peso del ciclo k en la matriz de densidad de una partícula
        weight = z2d * z_can[n_particles - k] / z_can[n_particles]
        if per_particle:
            weight /= n_particles
        val += weight * normalized_mehler_1d(x_arr, xp_arr, k * beta, omega)

    if np.ndim(x) == 0 and np.ndim(xp) == 0:
        return float(val)
    return val


def analytical_distinguishable_obdm_1d(x, xp, n_particles, temp, omega=OMEGA, per_particle=False):
    """Matriz de densidad marginal 1D para partículas distinguibles (sin intercambio, k=1)."""
    beta = 1.0 / temp
    norm = 1.0 if per_particle else float(n_particles)
    return norm * normalized_mehler_1d(x, xp, beta, omega)


def analytical_bose_endpoint_dist(s_vals, n_particles, temp, omega=OMEGA, tau=0.0):
    """
    Distribución analítica exacta de la separación espacial 2D s = |r_ira - r_masha|
    en el sector G (gap=1, tiempo efectivo k*beta - tau) para N bosones en trampa armónica 2D.
    """
    z_can = _full_bose_partition(n_particles, temp, omega)[0]
    beta = 1.0 / temp
    p_s = np.zeros_like(s_vals, dtype=float)
    w_sum = 0.0
    for k in range(1, n_particles + 1):
        q = np.exp(-k * beta * omega)
        z2d = 1.0 / (1.0 - q)**2
        u_eff = max(k * beta - tau, 1e-6) * omega
        weight = z2d * (z_can[n_particles - k] / z_can[n_particles]) * np.tanh(0.5 * u_eff)
        sigma2 = 2.0 * np.tanh(0.5 * u_eff) / omega
        w_sum += weight
        p_s += weight * (s_vals / sigma2) * np.exp(-s_vals**2 / (2.0 * sigma2))
    return p_s / w_sum if w_sum > 0 else p_s


def analytical_distinguishable_endpoint_dist(s_vals, temp, omega=OMEGA, tau=0.0):
    """Distribución de separación espacial 2D s para partículas distinguibles (k=1)."""
    beta = 1.0 / temp
    u_eff = max(beta - tau, 1e-6) * omega
    sigma2 = 2.0 * np.tanh(0.5 * u_eff) / omega
    return (s_vals / sigma2) * np.exp(-s_vals**2 / (2.0 * sigma2))


# ==============================================================================
# 2. Carga y Normalización de Datos PIMC
# ==============================================================================

def _symmetrize_if_needed(matrix):
    if not np.array_equal(matrix, matrix.T):
        return 0.5 * (matrix + matrix.T)
    return matrix


def load_worm_obdm_data(data_path):
    """
    Lee el archivo outworm_obdm.dat generado por el motor en C.
    Devuelve x_grid, matrix_norm, matrix_counts, metadata_dict.
    """
    if not os.path.isfile(data_path):
        raise FileNotFoundError(f"No se encontró el fichero de datos: {data_path}")

    metadata = {}
    with open(data_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("# grid:"):
                # Formato: # grid: 100 min: -3.5 max: 3.5 step: 0.07 samples: 1234 norm_trace: ... N: 2 T: 0.5 omega: 1 tau: ...
                parts = line.strip("# \n").split()
                for idx in range(0, len(parts) - 1, 2):
                    key = parts[idx].rstrip(":")
                    val = parts[idx + 1]
                    try:
                        if val.isdigit() or (val.startswith("-") and val[1:].isdigit()):
                            metadata[key] = int(val)
                        else:
                            metadata[key] = float(val)
                    except ValueError:
                        metadata[key] = val
                break

    grid_size = int(metadata.get("grid", 100))
    min_x = float(metadata.get("min", -3.5))
    max_x = float(metadata.get("max", 3.5))
    step = float(metadata.get("step", (max_x - min_x) / grid_size))

    # Cargar datos tabulares: col 0: x, col 1: x_prime, col 2: count, col 3: rho1_norm
    raw_data = np.loadtxt(data_path, comments="#", dtype=float, ndmin=2)
    if raw_data.size == 0:
        raise ValueError(f"El fichero {data_path} no contiene datos numéricos")

    x_vals = np.linspace(min_x + 0.5 * step, max_x - 0.5 * step, grid_size)

    if raw_data.shape[1] < 3:
        matrix_counts = np.zeros((grid_size, grid_size), dtype=float)
        matrix_norm = np.zeros((grid_size, grid_size), dtype=float)
    elif (
        raw_data.shape[0] == grid_size * grid_size
        and np.allclose(raw_data[:, 0], np.repeat(x_vals, grid_size), atol=1e-4 * max(step, 1e-6), rtol=0.0)
        and np.allclose(raw_data[:, 1], np.tile(x_vals, grid_size), atol=1e-4 * max(step, 1e-6), rtol=0.0)
    ):
        matrix_counts = _symmetrize_if_needed(raw_data[:, 2].reshape(grid_size, grid_size).copy())
        if raw_data.shape[1] >= 4:
            matrix_norm = _symmetrize_if_needed(raw_data[:, 3].reshape(grid_size, grid_size).copy())
        else:
            matrix_norm = np.zeros((grid_size, grid_size), dtype=float)
    else:
        matrix_counts = np.zeros((grid_size, grid_size), dtype=float)
        matrix_norm = np.zeros((grid_size, grid_size), dtype=float)
        ix = np.floor((raw_data[:, 0] - min_x) / step).astype(int)
        ixp = np.floor((raw_data[:, 1] - min_x) / step).astype(int)
        valid = (ix >= 0) & (ix < grid_size) & (ixp >= 0) & (ixp < grid_size)
        ix_v = ix[valid]
        ixp_v = ixp[valid]
        counts_v = raw_data[valid, 2]
        matrix_counts[ix_v, ixp_v] = counts_v
        matrix_counts[ixp_v, ix_v] = counts_v
        matrix_counts = _symmetrize_if_needed(matrix_counts)
        if raw_data.shape[1] >= 4:
            norm_v = raw_data[valid, 3]
            matrix_norm[ix_v, ixp_v] = norm_v
            matrix_norm[ixp_v, ix_v] = norm_v
            matrix_norm = _symmetrize_if_needed(matrix_norm)

    n_part = int(metadata.get("N", N_PARTICLES))
    raw_trace_integral = float(np.sum(np.diag(matrix_counts)) * step)
    if raw_data.shape[1] < 4 and raw_trace_integral > 0:
        matrix_norm = matrix_counts * (float(n_part) / raw_trace_integral)

    metadata["raw_trace_integral"] = raw_trace_integral
    metadata["trace_integral"] = float(np.sum(np.diag(matrix_norm)) * step)
    return x_vals, matrix_norm, matrix_counts, metadata


def load_worm_dist_data(dist_path):
    """Lee el archivo outworm_obdm_dist.dat."""
    if not os.path.isfile(dist_path):
        return None, None
    raw = np.loadtxt(dist_path, comments="#", dtype=float, ndmin=2)
    if raw.size == 0:
        return None, None
    s_vals = raw[:, 0]
    p_s = raw[:, 2] if raw.shape[1] >= 3 else raw[:, 1]
    return s_vals, p_s


def _find_local_param(data_dir, keys):
    """Busca un parámetro numérico en MC.cfg / MC.log / outworm_superfluid.dat del directorio."""
    candidates = [
        ("MC.cfg local", os.path.join(data_dir, "MC.cfg")),
        ("MC.log local", os.path.join(data_dir, "MC.log")),
        ("outworm_superfluid.dat", os.path.join(data_dir, "outworm_superfluid.dat")),
    ]
    if os.path.abspath(data_dir) != os.path.abspath(SCRIPT_DIR):
        candidates.extend([
            ("sim_run/MC.cfg", os.path.join(data_dir, "sim_run", "MC.cfg")),
            ("sim_run/MC.log", os.path.join(data_dir, "sim_run", "MC.log")),
        ])
    for label, path in candidates:
        for key in keys:
            val = _get_cfg_number(path, key)
            if val is not None:
                return val, f"{label} ({path})"
    return None, None


def resolve_obdm_parameters(data_dir, metadata, temp=None, n_particles=None, omega=None):
    """Resuelve (N, T, omega) priorizando: CLI -> cabecera outworm_obdm.dat -> MC.cfg/MC.log local -> MC.cfg raíz."""
    sources = {}

    if n_particles is not None:
        resolved_n = int(n_particles)
        sources["N"] = "CLI"
    elif "N" in metadata:
        resolved_n = int(metadata["N"])
        sources["N"] = "cabecera outworm_obdm.dat"
    else:
        local_n, src_n = _find_local_param(data_dir, ("N",))
        if local_n is not None:
            resolved_n = int(local_n)
            sources["N"] = src_n
        else:
            resolved_n = N_PARTICLES
            sources["N"] = "MC.cfg raíz (fallback)"

    if temp is not None:
        resolved_t = float(temp)
        sources["T"] = "CLI"
    elif "T" in metadata:
        resolved_t = float(metadata["T"])
        sources["T"] = "cabecera outworm_obdm.dat"
    else:
        local_t, src_t = _find_local_param(data_dir, ("T",))
        if local_t is not None:
            resolved_t = float(local_t)
            sources["T"] = src_t
        else:
            resolved_t = float(TEMP)
            sources["T"] = "MC.cfg raíz (fallback)"

    if omega is not None:
        resolved_omega = float(omega)
        sources["omega"] = "CLI"
    elif "omega" in metadata:
        resolved_omega = float(metadata["omega"])
        sources["omega"] = "cabecera outworm_obdm.dat"
    else:
        local_omega, src_omega = _find_local_param(data_dir, ("omega_x", "omega_z"))
        if local_omega is not None:
            resolved_omega = float(local_omega)
            sources["omega"] = src_omega
        else:
            resolved_omega = float(OMEGA)
            sources["omega"] = "MC.cfg raíz (fallback)"

    if resolved_n < 1 or resolved_t <= 0.0 or resolved_omega <= 0.0:
        raise ValueError(f"Parámetros inválidos: N={resolved_n}, T={resolved_t}, omega={resolved_omega}")

    return resolved_n, resolved_t, resolved_omega, sources


# ==============================================================================
# 3. Generación de Gráficos de Publicación
# ==============================================================================

def _save_and_sync_figure(fig, save_dir, filename, copy_overleaf=True):
    os.makedirs(save_dir, exist_ok=True)
    out_file = os.path.join(save_dir, filename)
    fig.savefig(out_file, dpi=300)
    plt.close(fig)
    print(f"  -> Guardado: {out_file}")

    if copy_overleaf and os.path.exists(OVERLEAF_IMG_DIR):
        shutil.copy2(out_file, os.path.join(OVERLEAF_IMG_DIR, filename))
    return out_file


def generar_mapa_calor_obdm(
    x_vals,
    pimc_matrix,
    n_particles=N_PARTICLES,
    temp=TEMP,
    save_dir=OUTPUT_DIR,
    omega=OMEGA,
    copy_overleaf=True,
):
    """
    Genera el panel comparativo de Mapas de Calor 2D:
      [1] PIMC Worm (Muestreado con cadenas abiertas)
      [2] Analítico Bose (Solución canónica exacta)
      [3] Residuos / Discrepancia absoluta
    """
    X, XP = np.meshgrid(x_vals, x_vals)
    bose_exact = analytical_bose_obdm_1d(X, XP, n_particles, temp, omega=omega)
    diff_matrix = pimc_matrix - bose_exact

    vmax = max(np.max(pimc_matrix), np.max(bose_exact))
    vmin = 0.0

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)

    # Panel 1: PIMC
    half_step = 0.5 * (x_vals[1] - x_vals[0])
    extent = [x_vals[0] - half_step, x_vals[-1] + half_step,
              x_vals[0] - half_step, x_vals[-1] + half_step]
    im0 = axes[0].imshow(pimc_matrix.T, extent=extent,
                         origin="lower", cmap="viridis", vmin=vmin, vmax=vmax)
    axes[0].set_title(r"$\rho_1(x, x')$ -- PIMC Algoritmo Worm", fontsize=12, fontweight="bold")
    axes[0].set_xlabel(r"Posición $x$ [$a_{\mathrm{ho}}$]", fontsize=11)
    axes[0].set_ylabel(r"Posición $x'$ [$a_{\mathrm{ho}}$]", fontsize=11)
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    # Panel 2: Analítico Bose
    im1 = axes[1].imshow(bose_exact.T, extent=extent,
                         origin="lower", cmap="viridis", vmin=vmin, vmax=vmax)
    axes[1].set_title(r"$\rho_1(x, x')$ -- Analítico Bose Canónico", fontsize=12, fontweight="bold")
    axes[1].set_xlabel(r"Posición $x$ [$a_{\mathrm{ho}}$]", fontsize=11)
    axes[1].set_ylabel(r"Posición $x'$ [$a_{\mathrm{ho}}$]", fontsize=11)
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    # Panel 3: Residuos
    max_err = max(np.max(np.abs(diff_matrix)), 1e-4)
    im2 = axes[2].imshow(diff_matrix.T, extent=extent,
                         origin="lower", cmap="coolwarm", vmin=-max_err, vmax=max_err)
    axes[2].set_title(r"Residuo: $\rho_1^{\mathrm{PIMC}} - \rho_1^{\mathrm{Bose}}$", fontsize=12, fontweight="bold")
    axes[2].set_xlabel(r"Posición $x$ [$a_{\mathrm{ho}}$]", fontsize=11)
    axes[2].set_ylabel(r"Posición $x'$ [$a_{\mathrm{ho}}$]", fontsize=11)
    fig.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)

    _save_and_sync_figure(
        fig, save_dir, "01_obdm_mapa_calor_2d_comparativo.png", copy_overleaf=copy_overleaf
    )


def generar_cortes_obdm(
    x_vals,
    pimc_matrix,
    n_particles=N_PARTICLES,
    temp=TEMP,
    save_dir=OUTPUT_DIR,
    omega=OMEGA,
    copy_overleaf=True,
):
    """
    Genera cortes 1D de la matriz de densidad:
      Panel A: Diagonal local rho_1(x, x) (perfil de densidad espacial).
      Panel B: Corte anti-diagonal rho_1(x, -x) (coherencia cuántica simétrica).
      Panel C: Coherencia espacial relativa rho_1(x, 0) respecto al origen.
    """
    diag_pimc = np.diag(pimc_matrix)
    diag_bose = analytical_bose_obdm_1d(x_vals, x_vals, n_particles, temp, omega=omega)
    diag_dist = analytical_distinguishable_obdm_1d(x_vals, x_vals, n_particles, temp, omega=omega)

    # Corte anti-diagonal: x' = -x
    antidiag_pimc = np.diag(np.fliplr(pimc_matrix))
    antidiag_bose = analytical_bose_obdm_1d(x_vals, -x_vals, n_particles, temp, omega=omega)
    antidiag_dist = analytical_distinguishable_obdm_1d(x_vals, -x_vals, n_particles, temp, omega=omega)

    # Corte x' = 0 (interpolación lineal vectorizada sobre la malla ordenada x_vals)
    if 0.0 <= x_vals[0]:
        cut_center_pimc = pimc_matrix[:, 0]
    elif 0.0 >= x_vals[-1]:
        cut_center_pimc = pimc_matrix[:, -1]
    else:
        hi = int(np.searchsorted(x_vals, 0.0))
        lo = hi - 1
        weight = (0.0 - x_vals[lo]) / (x_vals[hi] - x_vals[lo])
        cut_center_pimc = (1.0 - weight) * pimc_matrix[:, lo] + weight * pimc_matrix[:, hi]

    zeros_x = np.zeros_like(x_vals)
    cut_center_bose = analytical_bose_obdm_1d(x_vals, zeros_x, n_particles, temp, omega=omega)
    cut_center_dist = analytical_distinguishable_obdm_1d(x_vals, zeros_x, n_particles, temp, omega=omega)

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

    _save_and_sync_figure(
        fig, save_dir, "02_obdm_cortes_coherencia_cuantica.png", copy_overleaf=copy_overleaf
    )


def generar_grafico_distancia_extremos(
    dist_path,
    n_particles=N_PARTICLES,
    temp=TEMP,
    save_dir=OUTPUT_DIR,
    tau=0.0,
    omega=OMEGA,
    copy_overleaf=True,
):
    """Grafica la distribución de distancia espacial entre los extremos Ira y Masha."""
    s_vals, p_s = load_worm_dist_data(dist_path)
    if s_vals is None:
        print("  [Aviso] No se encontró outworm_obdm_dist.dat, saltando gráfica de distancia.")
        return

    fig, ax = plt.subplots(figsize=(7.5, 4.8), constrained_layout=True)
    mask = s_vals <= 5.0
    s_plot = s_vals[mask]
    p_bose = analytical_bose_endpoint_dist(s_plot, n_particles, temp, omega=omega, tau=tau)
    p_dist = analytical_distinguishable_endpoint_dist(s_plot, temp, omega=omega, tau=tau)

    ax.plot(s_plot, p_bose, color="crimson", linewidth=2.2, label="Analítico Bose Canónico", zorder=2)
    ax.plot(s_plot, p_dist, color="black", linestyle="--", linewidth=1.8, label="Distinguibles", zorder=1)
    ax.plot(s_plot, p_s[mask], "o", color="teal", markersize=4.0, alpha=0.85,
            label=r"PIMC: $p(s)$ extremos worm (gap == 1)", zorder=3)
    ax.set_title(r"Distribución de Separación de Extremos del Worm $s = |\mathbf{r}_{\mathrm{ira}} - \mathbf{r}_{\mathrm{masha}}|$",
                 fontsize=11, fontweight="bold")
    ax.set_xlabel(r"Distancia de separación $s$ [$a_{\mathrm{ho}}$]", fontsize=10)
    ax.set_ylabel(r"Densidad de probabilidad $p(s)$", fontsize=10)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, linestyle=":", alpha=0.5)

    _save_and_sync_figure(
        fig, save_dir, "03_obdm_distribucion_distancia_extremos.png", copy_overleaf=copy_overleaf
    )


# ==============================================================================
# 4. Función Principal de Análisis
# ==============================================================================

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "data_dir",
        nargs="?",
        default=SCRIPT_DIR,
        help="Directorio que contiene outworm_obdm.dat (por defecto: raíz del proyecto)",
    )
    parser.add_argument("-T", "--temp", type=float, default=None, help="Temperatura T para la referencia analítica")
    parser.add_argument("-N", "--particles", type=int, default=None, help="Número de partículas N")
    parser.add_argument("--omega", type=float, default=None, help="Frecuencia armónica omega")
    parser.add_argument("-o", "--output-dir", default=OUTPUT_DIR, help="Directorio donde guardar las figuras")
    parser.add_argument("--no-overleaf", action="store_true", help="No copiar figuras a overleaf/imágenes")
    args = parser.parse_args(argv)

    data_dir = os.path.abspath(args.data_dir)
    obdm_file = os.path.join(data_dir, "outworm_obdm.dat")
    dist_file = os.path.join(data_dir, "outworm_obdm_dist.dat")

    if not os.path.isfile(obdm_file):
        print(f"[Aviso] No se encontró '{obdm_file}'. Ejecute primero una simulación PIMC con Worm.")
        return 1

    print(f"\n==================================================================")
    print(f"  ANÁLISIS DE LA MATRIZ DE DENSIDAD DE UNA PARTÍCULA (OBDM)")
    print(f"==================================================================")
    print(f"  Leyendo archivo: {obdm_file}")

    x_vals, pimc_norm, _, meta = load_worm_obdm_data(obdm_file)
    n_part, temp, omega, sources = resolve_obdm_parameters(
        data_dir, meta, temp=args.temp, n_particles=args.particles, omega=args.omega
    )
    samples = meta.get("samples", 0)
    step = meta.get("step", 0.0)

    print(f"  -> Parámetros        : N={n_part} ({sources['N']}), T={temp:g} ({sources['T']}), omega={omega:g} ({sources['omega']})")
    print(f"  -> Puntos de cuadrícula: {len(x_vals)}")
    print(f"  -> Muestras gap == 1  : {samples}")
    print(f"  -> Rango espacial     : [{x_vals[0]:.2f}, {x_vals[-1]:.2f}] (paso dx = {step:.4f})")
    print(f"  -> Integral traza raw : {meta.get('trace_integral', 0.0):.2f}")

    copy_overleaf = not args.no_overleaf
    # 1. Mapa de Calor 2D
    generar_mapa_calor_obdm(
        x_vals, pimc_norm, n_particles=n_part, temp=temp,
        save_dir=args.output_dir, omega=omega, copy_overleaf=copy_overleaf,
    )

    # 2. Cortes Unidimensionales
    generar_cortes_obdm(
        x_vals, pimc_norm, n_particles=n_part, temp=temp,
        save_dir=args.output_dir, omega=omega, copy_overleaf=copy_overleaf,
    )

    # 3. Distribución de Distancia de Extremos
    tau_gap1 = float(meta.get("tau", 0.0))
    generar_grafico_distancia_extremos(
        dist_file, n_particles=n_part, temp=temp,
        save_dir=args.output_dir, tau=tau_gap1, omega=omega, copy_overleaf=copy_overleaf,
    )

    print(f"==================================================================\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
