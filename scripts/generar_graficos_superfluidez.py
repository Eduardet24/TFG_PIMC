#!/usr/bin/env python3
r"""
generar_graficos_superfluidez.py
================================
Generador de Figuras Científicas de Verificación y Análisis de la
Fracción Superfluida (\rho_s / \rho) en el Algoritmo Worm (PIMC 2D).

Estudia y grafica todo lo relevante para el análisis de superfluidez en trampa armónica 2D:
  1. Dependencia con la temperatura T (Bosones con Swap vs Distinguibles sin Swap)
     comparada con las soluciones analíticas exactas (canónica de Bose y Boltzmann).
  2. Panel de residuos normalizados (Pull plot, regla 3-sigma).
  3. Correlación física entre \rho_s / \rho y la fracción de partículas en ciclos
     de permutación no triviales (L >= 2).
  4. Convergencia MCMC por bloques y comprobación de simetría (<A_z> ~ 0).
  5. Convergencia con el número de beads P frente al límite discreto exacto y continuo.
  6. Escalado con el número de partículas N (Bosones vs Distinguibles).
  7. Panel resumen 2x2 tipo publicación.

Organización de salida en Resultados/Verificacion_Superfluidez/:
  - 00_panel_completo_superfluidez.{png,pdf}
  - 01_Dependencia_Temperatura/
  - 02_Correlacion_Permutaciones_y_MCMC/
  - 03_Convergencia_Beads_y_Particulas/
  - datos/ (datos numéricos y resumen JSON/CSV)

Uso:
  python3 generar_graficos_superfluidez.py [--no-run]
"""

import os
import re
import sys
import json
import time
import shutil
import argparse
import subprocess
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent.parent
BASE_OUT = SCRIPT_DIR / "Resultados" / "Verificacion_Superfluidez"
DIR_TEMP = BASE_OUT / "01_Dependencia_Temperatura"
DIR_PERM = BASE_OUT / "02_Correlacion_Permutaciones_y_MCMC"
DIR_CONV = BASE_OUT / "03_Convergencia_Beads_y_Particulas"
DIR_DATA = BASE_OUT / "datos"
OVERLEAF_DIR = SCRIPT_DIR.parent / "overleaf" / "imágenes"

for d in (BASE_OUT, DIR_TEMP, DIR_PERM, DIR_CONV, DIR_DATA):
    d.mkdir(parents=True, exist_ok=True)
try:
    OVERLEAF_DIR.mkdir(parents=True, exist_ok=True)
except OSError:
    pass

# Estilo científico limpio
plt.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "legend.fontsize": 9.5,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linestyle": "--",
})

COLOR_BOSE = "#d90429"
COLOR_DIST = "#2b2d42"
COLOR_PIMC_SWAP = "#0077b6"
COLOR_PIMC_NOSWAP = "#e65100"
COLOR_PERM = "#2a9d8f"
COLOR_ACCENT = "#7209b7"

OMEGA = 1.0


# ==============================================================================
# 1. SOLUCIONES ANALÍTICAS EXACTAS (CONTINUO Y DISCRETO EN P)
# ==============================================================================
def exact_distinguishable_continuum(temp, omega=OMEGA):
    """Fracción superfluida exacta para partículas distinguibles en trampa armónica 2D."""
    beta = 1.0 / np.asarray(temp, dtype=float)
    u = beta * omega
    return 1.0 - u / np.sinh(u)


def _cycle_properties_discrete(k, n_beads, temp, omega=OMEGA):
    """Propiedades exactas de un ciclo de longitud k con P=n_beads en trampa armónica 2D."""
    beta = 1.0 / temp
    tau = beta / n_beads
    M = k * n_beads
    m = np.arange(M)
    step = 2.0 * np.pi * m / M
    lam = 4.0 * np.sin(0.5 * step) ** 2 / tau + tau * omega ** 2

    # Función de partición de un ciclo en 2D (factor tau^(-M) común omitido)
    # En 1D es prod (dimless)^(-1/2), por tanto en 2D es prod (dimless)^(-1)
    dimless = 4.0 * np.sin(0.5 * step) ** 2 + (tau * omega) ** 2
    log_z1_2d = -np.sum(np.log(dimless))

    icl_k = np.sum(2.0 / lam) / n_beads
    az2_k = np.sum(np.sin(step[1:]) ** 2 / lam[1:] ** 2)
    return log_z1_2d, icl_k, az2_k


def _cycle_properties_continuum(k, temp, omega=OMEGA):
    """Propiedades exactas en el límite continuo (P -> inf) de un ciclo de longitud k en 2D."""
    u = k * omega / temp
    # z_1(k*beta) = 1 / (2 sinh(u/2))^2 = exp(-u) / (1 - exp(-u))^2
    log_z1_2d = -u - 2.0 * np.log1p(-np.exp(-u))
    icl_k = (k / omega) / np.tanh(0.5 * u)
    beta_k = k / temp
    az2_k = (beta_k / (4.0 * omega)) / np.tanh(0.5 * u) - (beta_k ** 2) / (8.0 * np.sinh(0.5 * u) ** 2)
    return log_z1_2d, icl_k, az2_k


def exact_bose_superfluid(n_particles, temp, n_beads=None, omega=OMEGA):
    """
    Calcula exactamente <A_z^2>, <I_cl> y rho_s/rho = 4*T*<A_z^2>/<I_cl>
    para N bosones no interactuantes en una trampa armónica 2D canónica.
    Si n_beads is None, usa el límite continuo (P -> inf).
    """
    Z = np.zeros(n_particles + 1)
    Z[0] = 1.0
    icl_num = np.zeros(n_particles + 1)
    az2_num = np.zeros(n_particles + 1)

    props = {}
    log_scales = []
    for k in range(1, n_particles + 1):
        if n_beads is None:
            lz, icl_k, az2_k = _cycle_properties_continuum(k, temp, omega)
        else:
            lz, icl_k, az2_k = _cycle_properties_discrete(k, n_beads, temp, omega)
        props[k] = (lz, icl_k, az2_k)
        log_scales.append(lz / k)

    # Reescalado numérico para evitar overflow en productos de partición
    c_shift = max(log_scales)
    z1 = {k: np.exp(props[k][0] - k * c_shift) for k in range(1, n_particles + 1)}

    for count in range(1, n_particles + 1):
        for k in range(1, count + 1):
            _, icl_k, az2_k = props[k]
            w = z1[k] / count
            Z[count] += w * Z[count - k]
            icl_num[count] += w * (icl_k * Z[count - k] + icl_num[count - k])
            az2_num[count] += w * (az2_k * Z[count - k] + az2_num[count - k])

    icl_mean = icl_num[n_particles] / Z[n_particles]
    az2_mean = az2_num[n_particles] / Z[n_particles]
    rhos = (4.0 * temp / icl_mean) * az2_mean
    return rhos, az2_mean, icl_mean


# ==============================================================================
# 2. COMPILACIÓN Y EJECUCIÓN DE SIMULACIONES PIMC
# ==============================================================================
def get_or_compile_mc():
    exe = SCRIPT_DIR / "MC"
    sources = sorted(SCRIPT_DIR.glob("*.c"))
    newest_src = max(s.stat().st_mtime for s in sources)
    if not exe.is_file() or exe.stat().st_mtime < newest_src:
        cmd = ["gcc", "-O2", "-fopenmp", "-std=c99", *[str(s) for s in sources], "-lm", "-o", str(exe)]
        print(f"  [Compilando] {' '.join(cmd[:5])} ...")
        res = subprocess.run(cmd, cwd=SCRIPT_DIR, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"Error compilando MC:\n{res.stderr}")
    return exe


def run_pimc_case(exe, base_cfg_text, run_dir, n_particles, temp, n_beads,
                  blck_heating=25, blck=120, n_iter=1500, swap_enabled=True):
    run_dir.mkdir(parents=True, exist_ok=True)
    # Ajustar worm_C para mantener fracción Z saludable (~10-25%)
    worm_c = 0.03 / max(1, n_particles * n_beads * 0.25)

    cfg = base_cfg_text
    cfg = re.sub(r"N=\s*\d+", f"N= {n_particles}", cfg)
    cfg = re.sub(r"T=\s*[\d\.eE+-]+", f"T= {temp}", cfg)
    cfg = re.sub(r"Npop=\s*\d+", f"Npop= {n_beads}", cfg)
    cfg = re.sub(r"NwalkersMax=\s*\d+", f"NwalkersMax= {max(60, n_beads * 2)}", cfg)
    cfg = re.sub(r"worm_C=\s*[\d\.eE+-]+", f"worm_C= {worm_c:.6f}", cfg)
    cfg = re.sub(r"blck_heating=\s*\d+", f"blck_heating= {blck_heating}", cfg)
    cfg = re.sub(r"blck=\s*[\d\s]+", f"blck= {blck}\n", cfg)
    cfg = re.sub(r"Niter=\s*\d+", f"Niter= {n_iter}", cfg)
    cfg = re.sub(r"Nmeasure=\s*\d+", "Nmeasure= 1", cfg)

    cfg_file = run_dir / "MC.cfg"
    cfg_file.write_text(cfg, encoding="utf-8")

    env = os.environ.copy()
    env["WORM_DISABLE_SWAP"] = "0" if swap_enabled else "1"

    proc = subprocess.run([str(exe), "MC.cfg"], cwd=run_dir, env=env, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Error en simulación {run_dir}:\n{proc.stderr[-500:]}")

    sf_file = run_dir / "outworm_superfluid.dat"
    sf_blocks_file = run_dir / "outworm_superfluid_blocks.dat"
    perm_file = run_dir / "outworm_perm.dat"

    sf_vals = np.loadtxt(sf_file, comments="#")
    sf_blocks = np.loadtxt(sf_blocks_file, comments="#")
    perm_vals = np.loadtxt(perm_file, comments="#")
    if perm_vals.ndim == 1:
        perm_vals = perm_vals.reshape(1, -1)

    p1 = float(perm_vals[0, 2])
    p_non_trivial = float(np.sum(perm_vals[1:, 2])) if perm_vals.shape[0] > 1 else 0.0

    return {
        "N": n_particles,
        "T": temp,
        "P": n_beads,
        "swap": swap_enabled,
        "rho_s": float(sf_vals[0]),
        "rho_s_block_mean": float(sf_vals[1]),
        "rho_s_err": float(sf_vals[2]),
        "A2_mean": float(sf_vals[3]),
        "A_mean": float(sf_vals[4]),
        "Icl_mean": float(sf_vals[5]),
        "P1": p1,
        "P_nontrivial": p_non_trivial,
        "blocks_file": str(sf_blocks_file),
    }


def save_figure(fig, path_without_ext):
    png_path = Path(f"{path_without_ext}.png")
    pdf_path = Path(f"{path_without_ext}.pdf")
    fig.savefig(png_path, bbox_inches="tight", dpi=300)
    fig.savefig(pdf_path, bbox_inches="tight")
    if OVERLEAF_DIR.is_dir():
        try:
            shutil.copy2(png_path, OVERLEAF_DIR / png_path.name)
            shutil.copy2(pdf_path, OVERLEAF_DIR / pdf_path.name)
        except OSError:
            pass
    plt.close(fig)
    print(f"  [Figura] {png_path.relative_to(SCRIPT_DIR)}")


# ==============================================================================
# 3. GENERACIÓN DE GRÁFICOS
# ==============================================================================
def plot_temperature_dependence(data_temp):
    n_part = data_temp["N"]
    n_beads = data_temp["P"]
    temps = np.array([d["T"] for d in data_temp["swap"]])
    rhos_swap = np.array([d["rho_s_block_mean"] for d in data_temp["swap"]])
    err_swap = np.array([d["rho_s_err"] for d in data_temp["swap"]])
    rhos_noswap = np.array([d["rho_s_block_mean"] for d in data_temp["noswap"]])
    err_noswap = np.array([d["rho_s_err"] for d in data_temp["noswap"]])

    t_grid = np.linspace(min(temps) * 0.85, max(temps) * 1.05, 200)
    exact_bose_grid = np.array([exact_bose_superfluid(n_part, t, n_beads=n_beads)[0] for t in t_grid])
    exact_dist_grid = np.array([exact_bose_superfluid(1, t, n_beads=n_beads)[0] for t in t_grid])

    exact_bose_pts = np.array([exact_bose_superfluid(n_part, t, n_beads=n_beads)[0] for t in temps])
    exact_dist_pts = np.array([exact_bose_superfluid(1, t, n_beads=n_beads)[0] for t in temps])

    # --- Figura 1: rho_s/rho vs T (simple) ---
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    ax.plot(t_grid, exact_bose_grid, color=COLOR_BOSE, lw=2.2,
            label=f"Teoría exacta Bose (N={n_part}, P={n_beads})")
    ax.plot(t_grid, exact_dist_grid, color=COLOR_DIST, lw=2.0, ls="--",
            label=f"Teoría exacta Distinguibles (P={n_beads})")
    ax.errorbar(temps, rhos_swap, yerr=err_swap, fmt="o", color=COLOR_PIMC_SWAP,
                ms=6.5, capsize=3.5, elinewidth=1.4, label="Worm PIMC — Bosones (Swap)")
    ax.errorbar(temps, rhos_noswap, yerr=err_noswap, fmt="s", color=COLOR_PIMC_NOSWAP,
                ms=6.0, capsize=3.5, elinewidth=1.4, label="Worm PIMC — Distinguibles (Sin Swap)")

    ax.set_xlabel(r"Temperatura $T$ ($\hbar\omega / k_B$)")
    ax.set_ylabel(r"Fracción superfluida $\rho_s / \rho$")
    ax.set_title(f"Fracción superfluida vs Temperatura (Trampa armónica 2D, N={n_part})")
    ax.set_ylim(-0.02, 1.05)
    ax.legend(frameon=True, facecolor="white", framealpha=0.92)
    save_figure(fig, DIR_TEMP / "01_rhos_vs_T_bosones_vs_distinguibles")

    # --- Figura 2: rho_s/rho vs T con panel de residuos normalizados (pull plot) ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.6, 6.8), sharex=True,
                                   gridspec_kw={"height_ratios": [2.4, 1.0], "hspace": 0.08})
    ax1.plot(t_grid, exact_bose_grid, color=COLOR_BOSE, lw=2.2, label="Teoría exacta Bose")
    ax1.plot(t_grid, exact_dist_grid, color=COLOR_DIST, lw=2.0, ls="--", label="Teoría exacta Distinguibles")
    ax1.errorbar(temps, rhos_swap, yerr=err_swap, fmt="o", color=COLOR_PIMC_SWAP,
                 ms=6.5, capsize=3.5, label="PIMC Bosones (Swap)")
    ax1.errorbar(temps, rhos_noswap, yerr=err_noswap, fmt="s", color=COLOR_PIMC_NOSWAP,
                 ms=6.0, capsize=3.5, label="PIMC Distinguibles (Sin Swap)")
    ax1.set_ylabel(r"$\rho_s / \rho$")
    ax1.set_title(f"Verificación térmica de la fracción superfluida y residuos ($N={n_part}$, $P={n_beads}$)")
    ax1.set_ylim(-0.02, 1.05)
    ax1.legend(frameon=True, facecolor="white", framealpha=0.92)

    pull_swap = (rhos_swap - exact_bose_pts) / np.maximum(err_swap, 1e-12)
    pull_noswap = (rhos_noswap - exact_dist_pts) / np.maximum(err_noswap, 1e-12)
    ax2.axhspan(-3, 3, color="#e9ecef", alpha=0.7, label=r"Banda $\pm 3\sigma$")
    ax2.axhspan(-1, 1, color="#dee2e6", alpha=0.8, label=r"Banda $\pm 1\sigma$")
    ax2.axhline(0.0, color="black", lw=1.0, ls=":")
    ax2.plot(temps, pull_swap, "o-", color=COLOR_PIMC_SWAP, ms=5.5, lw=1.1, label="Pull Bosones")
    ax2.plot(temps, pull_noswap, "s--", color=COLOR_PIMC_NOSWAP, ms=5.0, lw=1.1, label="Pull Distinguibles")
    ax2.set_xlabel(r"Temperatura $T$")
    ax2.set_ylabel(r"Residuo $(y - y_{\mathrm{ex}})/\sigma$")
    ax2.set_ylim(-4.2, 4.2)
    ax2.legend(loc="upper right", ncol=2, fontsize=8.5, frameon=True)
    save_figure(fig, DIR_TEMP / "02_rhos_vs_T_con_residuos")


def plot_permutations_and_mcmc(data_temp, mcmc_case):
    temps = np.array([d["T"] for d in data_temp["swap"]])
    rhos_swap = np.array([d["rho_s_block_mean"] for d in data_temp["swap"]])
    err_swap = np.array([d["rho_s_err"] for d in data_temp["swap"]])
    rhos_noswap = np.array([d["rho_s_block_mean"] for d in data_temp["noswap"]])
    delta_rhos = rhos_swap - rhos_noswap
    p_nontriv = np.array([d["P_nontrivial"] for d in data_temp["swap"]])

    # --- Figura 3: Correlación entre exceso superfluido y ciclos de permutación ---
    fig, ax1 = plt.subplots(figsize=(7.6, 5.2))
    l1 = ax1.errorbar(temps, rhos_swap, yerr=err_swap, fmt="o-", color=COLOR_PIMC_SWAP,
                      ms=6.5, lw=1.8, capsize=3.5, label=r"$\rho_s / \rho$ (Bosones)")
    l2 = ax1.plot(temps, delta_rhos, "d--", color=COLOR_ACCENT, ms=6.0, lw=1.6,
                  label=r"Exceso bosónico $\Delta(\rho_s/\rho)_{\mathrm{Swap}}$")
    ax1.set_xlabel(r"Temperatura $T$")
    ax1.set_ylabel(r"Fracción superfluida $\rho_s / \rho$", color=COLOR_PIMC_SWAP)
    ax1.tick_params(axis="y", labelcolor=COLOR_PIMC_SWAP)
    ax1.set_ylim(-0.02, 1.05)

    ax2 = ax1.twinx()
    l3 = ax2.plot(temps, p_nontriv, "s-", color=COLOR_PERM, ms=6.0, lw=1.8,
                  label=r"Fracción en ciclos $L \geq 2$")
    ax2.set_ylabel(r"Fracción de partículas en intercambio $P(L \geq 2)$", color=COLOR_PERM)
    ax2.tick_params(axis="y", labelcolor=COLOR_PERM)
    ax2.set_ylim(-0.02, 1.05)
    ax2.grid(False)

    lines = [l1] + l2 + l3
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc="upper right", frameon=True, facecolor="white", framealpha=0.92)
    ax1.set_title(f"Origen bosónico de la superfluidez: ciclos de permutación (N={data_temp['N']})")
    save_figure(fig, DIR_PERM / "03_rhos_vs_ciclos_permutacion")

    # --- Figura 4: Convergencia MCMC por bloques y chequeo de simetría <A_z> ~ 0 ---
    blk_data = np.loadtxt(mcmc_case["blocks_file"], comments="#")
    blocks = blk_data[:, 0]
    samples = blk_data[:, 1]
    rho_blk = blk_data[:, 2]
    a_blk = blk_data[:, 4]

    cum_weights = np.cumsum(samples)
    cum_rho = np.cumsum(rho_blk * samples) / cum_weights
    cum_a = np.cumsum(a_blk * samples) / cum_weights
    exact_val, _, _ = exact_bose_superfluid(mcmc_case["N"], mcmc_case["T"], n_beads=mcmc_case["P"])

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(7.6, 6.4), sharex=True,
                                         gridspec_kw={"height_ratios": [1.6, 1.0], "hspace": 0.1})
    ax_top.plot(blocks, rho_blk, ".", color="#90e0ef", alpha=0.75, label="Estimación por bloque")
    ax_top.plot(blocks, cum_rho, "-", color=COLOR_PIMC_SWAP, lw=2.2, label="Media acumulada PIMC")
    ax_top.axhline(exact_val, color=COLOR_BOSE, lw=2.0, ls="--", label=f"Valor exacto ({exact_val:.4f})")
    ax_top.set_ylabel(r"$\rho_s / \rho$")
    ax_top.set_title(f"Convergencia MCMC y simetría temporal (N={mcmc_case['N']}, T={mcmc_case['T']}, P={mcmc_case['P']})")
    ax_top.legend(frameon=True, facecolor="white", framealpha=0.92)

    ax_bot.axhline(0.0, color="black", lw=1.2, ls="--", label=r"Simetría exacta $\langle A_z \rangle = 0$")
    ax_bot.plot(blocks, a_blk, ".", color="#adb5bd", alpha=0.7, label=r"$A_z$ por bloque")
    ax_bot.plot(blocks, cum_a, "-", color=COLOR_PERM, lw=2.0, label=r"Media acumulada $\langle A_z \rangle$")
    ax_bot.set_xlabel("Bloque de medición")
    ax_bot.set_ylabel(r"Área media $\langle A_z \rangle$")
    ax_bot.legend(frameon=True, facecolor="white", framealpha=0.92, ncol=2, fontsize=8.5)
    save_figure(fig, DIR_PERM / "04_convergencia_mcmc_y_simetria")


def plot_convergence_and_n_scaling(data_beads, data_n):
    # --- Figura 5: Convergencia con número de beads P ---
    p_vals = np.array([d["P"] for d in data_beads["runs"]])
    rhos_p = np.array([d["rho_s_block_mean"] for d in data_beads["runs"]])
    err_p = np.array([d["rho_s_err"] for d in data_beads["runs"]])

    n_part = data_beads["N"]
    temp = data_beads["T"]
    p_grid = np.arange(min(p_vals), max(p_vals) + 5, 2)
    exact_disc_p = np.array([exact_bose_superfluid(n_part, temp, n_beads=int(p))[0] for p in p_grid])
    exact_cont, _, _ = exact_bose_superfluid(n_part, temp, n_beads=None)

    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    ax.axhline(exact_cont, color=COLOR_BOSE, lw=2.0, ls=":", label=f"Límite continuo $P\\to\\infty$ ({exact_cont:.4f})")
    ax.plot(p_grid, exact_disc_p, color=COLOR_DIST, lw=2.0, ls="--", label="Teoría exacta discreta $P$")
    ax.errorbar(p_vals, rhos_p, yerr=err_p, fmt="o", color=COLOR_PIMC_SWAP,
                ms=7.0, capsize=4.0, elinewidth=1.5, label="Worm PIMC (Bosones)")
    ax.set_xlabel("Número de beads $P$")
    ax.set_ylabel(r"Fracción superfluida $\rho_s / \rho$")
    ax.set_title(f"Convergencia de Trotter con el número de beads (N={n_part}, T={temp})")
    ax.legend(frameon=True, facecolor="white", framealpha=0.92)
    save_figure(fig, DIR_CONV / "05_convergencia_beads_P")

    # --- Figura 6: Escalado con el número de partículas N ---
    n_vals = np.array([d["N"] for d in data_n["swap"]])
    rhos_n_swap = np.array([d["rho_s_block_mean"] for d in data_n["swap"]])
    err_n_swap = np.array([d["rho_s_err"] for d in data_n["swap"]])
    rhos_n_noswap = np.array([d["rho_s_block_mean"] for d in data_n["noswap"]])
    err_n_noswap = np.array([d["rho_s_err"] for d in data_n["noswap"]])

    n_grid = np.arange(1, max(n_vals) + 1)
    temp_n = data_n["T"]
    p_n = data_n["P"]
    exact_n_bose = np.array([exact_bose_superfluid(int(nn), temp_n, n_beads=p_n)[0] for nn in n_grid])
    exact_n_dist = np.full_like(n_grid, exact_bose_superfluid(1, temp_n, n_beads=p_n)[0], dtype=float)

    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    ax.plot(n_grid, exact_n_bose, "o-", color=COLOR_BOSE, lw=2.0, ms=4.5, label="Teoría exacta Bose")
    ax.plot(n_grid, exact_n_dist, "--", color=COLOR_DIST, lw=2.0, label="Teoría exacta Distinguibles")
    ax.errorbar(n_vals, rhos_n_swap, yerr=err_n_swap, fmt="o", color=COLOR_PIMC_SWAP,
                ms=7.0, capsize=4.0, label="Worm PIMC — Bosones (Swap)")
    ax.errorbar(n_vals, rhos_n_noswap, yerr=err_n_noswap, fmt="s", color=COLOR_PIMC_NOSWAP,
                ms=6.5, capsize=4.0, label="Worm PIMC — Distinguibles (Sin Swap)")
    ax.set_xlabel("Número de partículas $N$")
    ax.set_ylabel(r"Fracción superfluida $\rho_s / \rho$")
    ax.set_title(f"Dependencia con el número de partículas N (T={temp_n}, P={p_n})")
    ax.set_xticks(n_grid)
    ax.set_ylim(0.0, 1.02)
    ax.legend(frameon=True, facecolor="white", framealpha=0.92)
    save_figure(fig, DIR_CONV / "06_dependencia_N_particulas")


def plot_dashboard(data_temp, mcmc_case, data_beads, data_n):
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.8))
    ((ax1, ax2), (ax3, ax4)) = axes

    # Panel A: rho_s vs T
    n_part, n_beads = data_temp["N"], data_temp["P"]
    temps = np.array([d["T"] for d in data_temp["swap"]])
    rhos_s = np.array([d["rho_s_block_mean"] for d in data_temp["swap"]])
    err_s = np.array([d["rho_s_err"] for d in data_temp["swap"]])
    rhos_ns = np.array([d["rho_s_block_mean"] for d in data_temp["noswap"]])
    err_ns = np.array([d["rho_s_err"] for d in data_temp["noswap"]])
    t_grid = np.linspace(min(temps) * 0.85, max(temps) * 1.05, 160)
    ax1.plot(t_grid, [exact_bose_superfluid(n_part, t, n_beads)[0] for t in t_grid], color=COLOR_BOSE, lw=2.0, label="Exacta Bose")
    ax1.plot(t_grid, [exact_bose_superfluid(1, t, n_beads)[0] for t in t_grid], color=COLOR_DIST, lw=1.8, ls="--", label="Exacta Dist.")
    ax1.errorbar(temps, rhos_s, yerr=err_s, fmt="o", color=COLOR_PIMC_SWAP, ms=5.5, capsize=3, label="PIMC Swap")
    ax1.errorbar(temps, rhos_ns, yerr=err_ns, fmt="s", color=COLOR_PIMC_NOSWAP, ms=5.0, capsize=3, label="PIMC No-Swap")
    ax1.set_title(f"(a) Fracción superfluida vs T (N={n_part})")
    ax1.set_xlabel("Temperatura $T$")
    ax1.set_ylabel(r"$\rho_s / \rho$")
    ax1.set_ylim(-0.02, 1.05)
    ax1.legend(fontsize=8.5, frameon=True)

    # Panel B: Correlación con ciclos de permutación
    p_nontriv = np.array([d["P_nontrivial"] for d in data_temp["swap"]])
    ax2.errorbar(temps, rhos_s, yerr=err_s, fmt="o-", color=COLOR_PIMC_SWAP, ms=5.5, lw=1.5, capsize=3, label=r"$\rho_s/\rho$ (Bosones)")
    ax2.plot(temps, p_nontriv, "s-", color=COLOR_PERM, ms=5.5, lw=1.6, label=r"Fracción ciclos $L \geq 2$")
    ax2.plot(temps, rhos_s - rhos_ns, "d--", color=COLOR_ACCENT, ms=5.0, lw=1.4, label=r"$\Delta(\rho_s/\rho)_{\mathrm{Swap}}$")
    ax2.set_title("(b) Superfluidez y ciclos de intercambio")
    ax2.set_xlabel("Temperatura $T$")
    ax2.set_ylabel("Fracción")
    ax2.set_ylim(-0.02, 1.05)
    ax2.legend(fontsize=8.5, frameon=True)

    # Panel C: Convergencia en beads P
    p_vals = np.array([d["P"] for d in data_beads["runs"]])
    rhos_p = np.array([d["rho_s_block_mean"] for d in data_beads["runs"]])
    err_p = np.array([d["rho_s_err"] for d in data_beads["runs"]])
    p_grid = np.arange(min(p_vals), max(p_vals) + 4, 2)
    exact_cont, _, _ = exact_bose_superfluid(data_beads["N"], data_beads["T"], n_beads=None)
    ax3.axhline(exact_cont, color=COLOR_BOSE, lw=1.8, ls=":", label=f"Continuo ({exact_cont:.3f})")
    ax3.plot(p_grid, [exact_bose_superfluid(data_beads["N"], data_beads["T"], int(p))[0] for p in p_grid],
             color=COLOR_DIST, lw=1.8, ls="--", label="Discreto exacto")
    ax3.errorbar(p_vals, rhos_p, yerr=err_p, fmt="o", color=COLOR_PIMC_SWAP, ms=6.0, capsize=3, label="PIMC Bosones")
    ax3.set_title(f"(c) Convergencia con beads P (N={data_beads['N']}, T={data_beads['T']})")
    ax3.set_xlabel("Número de beads $P$")
    ax3.set_ylabel(r"$\rho_s / \rho$")
    ax3.legend(fontsize=8.5, frameon=True)

    # Panel D: Escalado con N partículas
    n_vals = np.array([d["N"] for d in data_n["swap"]])
    n_grid = np.arange(1, max(n_vals) + 1)
    ax4.plot(n_grid, [exact_bose_superfluid(int(nn), data_n["T"], data_n["P"])[0] for nn in n_grid],
             "o-", color=COLOR_BOSE, lw=1.8, ms=4, label="Exacta Bose")
    ax4.axhline(exact_bose_superfluid(1, data_n["T"], data_n["P"])[0], color=COLOR_DIST, lw=1.8, ls="--", label="Exacta Dist.")
    ax4.errorbar(n_vals, [d["rho_s_block_mean"] for d in data_n["swap"]],
                 yerr=[d["rho_s_err"] for d in data_n["swap"]], fmt="o", color=COLOR_PIMC_SWAP, ms=6, capsize=3, label="PIMC Swap")
    ax4.errorbar(n_vals, [d["rho_s_block_mean"] for d in data_n["noswap"]],
                 yerr=[d["rho_s_err"] for d in data_n["noswap"]], fmt="s", color=COLOR_PIMC_NOSWAP, ms=5.5, capsize=3, label="PIMC No-Swap")
    ax4.set_title(f"(d) Escalado con N partículas (T={data_n['T']}, P={data_n['P']})")
    ax4.set_xlabel("Número de partículas $N$")
    ax4.set_ylabel(r"$\rho_s / \rho$")
    ax4.set_xticks(n_grid)
    ax4.set_ylim(0.0, 1.02)
    ax4.legend(fontsize=8.5, frameon=True)

    fig.suptitle("Resumen de Verificación: Fracción Superfluida en Worm PIMC (Trampa Armónica 2D)",
                 fontsize=14, fontweight="bold", y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    save_figure(fig, BASE_OUT / "00_panel_completo_superfluidez")


# ==============================================================================
# 4. ORQUESTACIÓN PRINCIPAL
# ==============================================================================
def main():
    parser = argparse.ArgumentParser(description="Generar gráficos de fracción superfluida Worm PIMC")
    parser.add_argument("--no-run", action="store_true", help="Reutilizar datos JSON existentes sin re-ejecutar C")
    args = parser.parse_args()

    summary_json = DIR_DATA / "resumen_superfluidez.json"

    if args.no_run and summary_json.is_file():
        print("[1/2] Cargando datos existentes desde resumen_superfluidez.json...")
        dataset = json.loads(summary_json.read_text(encoding="utf-8"))
    else:
        print("=================================================================")
        print("  ESTUDIO COMPLETO DE FRACCIÓN SUPERFLUIDA (WORM PIMC 2D)")
        print("=================================================================")
        exe = get_or_compile_mc()
        base_cfg = (SCRIPT_DIR / "MC.cfg").read_text(encoding="utf-8", errors="ignore")
        t0 = time.time()

        # 1. Barrido en temperatura (N=4, P=30)
        n_temp = 4
        p_temp = 30
        temps = [0.25, 0.35, 0.5, 0.7, 1.0, 1.4, 2.0, 2.6]
        swap_runs, noswap_runs = [], []
        print(f"[1/3] Barrido en temperatura T ({len(temps)} puntos, N={n_temp}, P={p_temp})...")
        for t in temps:
            beads = 40 if t < 0.4 else p_temp
            r_s = run_pimc_case(exe, base_cfg, DIR_DATA / f"T_{t:g}_swap", n_temp, t, beads,
                                blck_heating=25, blck=100, n_iter=1500, swap_enabled=True)
            r_ns = run_pimc_case(exe, base_cfg, DIR_DATA / f"T_{t:g}_noswap", n_temp, t, beads,
                                 blck_heating=25, blck=100, n_iter=1500, swap_enabled=False)
            swap_runs.append(r_s)
            noswap_runs.append(r_ns)
            print(f"      T={t:4.2f} | Swap: {r_s['rho_s_block_mean']:.4f}+/-{r_s['rho_s_err']:.4f} | "
                  f"NoSwap: {r_ns['rho_s_block_mean']:.4f}+/-{r_ns['rho_s_err']:.4f} | P(L>=2)={r_s['P_nontrivial']:.3f}")

        # 2. Convergencia en beads P (N=4, T=0.5)
        p_list = [10, 20, 30, 50, 80]
        beads_runs = []
        print(f"[2/3] Convergencia en número de beads P ({p_list}, N=4, T=0.5)...")
        for p in p_list:
            r_p = run_pimc_case(exe, base_cfg, DIR_DATA / f"P_{p}", 4, 0.5, p,
                                blck_heating=25, blck=100, n_iter=1500, swap_enabled=True)
            beads_runs.append(r_p)
            print(f"      P={p:2d} | rho_s/rho = {r_p['rho_s_block_mean']:.4f}+/-{r_p['rho_s_err']:.4f}")

        # 3. Escalado con N partículas (T=0.5, P=30)
        n_list = [1, 2, 3, 5, 8]
        n_swap_runs, n_noswap_runs = [], []
        print(f"[3/3] Escalado con número de partículas N ({n_list}, T=0.5, P=30)...")
        for nn in n_list:
            r_ns_s = run_pimc_case(exe, base_cfg, DIR_DATA / f"N_{nn}_swap", nn, 0.5, 30,
                                   blck_heating=25, blck=100, n_iter=1500, swap_enabled=True)
            r_ns_ns = run_pimc_case(exe, base_cfg, DIR_DATA / f"N_{nn}_noswap", nn, 0.5, 30,
                                    blck_heating=25, blck=100, n_iter=1500, swap_enabled=False)
            n_swap_runs.append(r_ns_s)
            n_noswap_runs.append(r_ns_ns)
            print(f"      N={nn} | Swap: {r_ns_s['rho_s_block_mean']:.4f}+/-{r_ns_s['rho_s_err']:.4f} | "
                  f"NoSwap: {r_ns_ns['rho_s_block_mean']:.4f}+/-{r_ns_ns['rho_s_err']:.4f}")

        dataset = {
            "temp_sweep": {"N": n_temp, "P": p_temp, "swap": swap_runs, "noswap": noswap_runs},
            "beads_sweep": {"N": 4, "T": 0.5, "runs": beads_runs},
            "n_sweep": {"T": 0.5, "P": 30, "swap": n_swap_runs, "noswap": n_noswap_runs},
        }
        summary_json.write_text(json.dumps(dataset, indent=2), encoding="utf-8")
        print(f"\n  Simulaciones completadas en {time.time() - t0:.1f} s.")

    print("\nGenerando figuras científicas en Resultados/Verificacion_Superfluidez/...")
    mcmc_case = dataset["temp_sweep"]["swap"][2]  # T = 0.5
    plot_temperature_dependence(dataset["temp_sweep"])
    plot_permutations_and_mcmc(dataset["temp_sweep"], mcmc_case)
    plot_convergence_and_n_scaling(dataset["beads_sweep"], dataset["n_sweep"])
    plot_dashboard(dataset["temp_sweep"], mcmc_case, dataset["beads_sweep"], dataset["n_sweep"])
    print("\n¡Proceso completado!")


if __name__ == "__main__":
    main()
