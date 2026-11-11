#!/usr/bin/env python3
"""
generar_graficos_energia.py
===========================
Generador de Figuras Científicas de Verificación, Comprobación y Demostración
de la Energía en el Algoritmo Worm (Path Integral Monte Carlo - PIMC 2D).

Objetivos físicos y metodológicos:
  1. DEMOSTRACIÓN CUÁNTICA (Efecto Swap):
     Evidenciar el salto cuántico bosónico de energía entre partículas distinguibles
     (sin Swap, Boltzmann) y bosones idénticos (con Swap, Bose-Einstein).
  2. COMPROBACIÓN FÍSICA (Teorema del Virial Cuántico):
     Comprobar la equipartición virial <E_pot> = <E_kin> demostrando la validez
     del estimador virial centroidal por ciclos implementado en C (quantities.c).
  3. VERIFICACIÓN ESTADÍSTICA (Convergencia MCMC y Regla 3-Sigma):
     Verificar la ausencia de sesgo en el Sector Z mediante medias acumuladas
     ponderadas y panel de residuos normalizados (Pull plot |z| <= 3).
  4. VERIFICACIÓN TERMODINÁMICA (Curva E(T) y Fracción Condensada):
     Comprobar la energía en todo el rango térmico T in [0.2, 3.0], desde el estado
     fundamental (T -> 0, E/N -> 1.0) hasta el límite clásico (T >> 1, E/N -> 2T),
    correlacionada con la ocupación fundamental N0/N y los ciclos no triviales.
  5. DIVISIBILIDAD TOTAL:
     Genera figuras modulares individuales e independientes (formato simple y
     formato completo con residuos) y un Dashboard / Panel general 2x2.
     Cada figura se guarda automáticamente en PNG (300 DPI) y PDF vectorial en:
       - Resultados/Verificacion_Energia/
       - ../overleaf/imágenes/

Uso:
  python3 generar_graficos_energia.py [--particles 5] [--run]
"""

import os
import sys
import re
import shutil
import argparse
import subprocess
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# CONFIGURACIÓN DE RUTAS Y CONSTANTES
# ==============================================================================
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR
RESULTS_DIR = SCRIPT_DIR / "Resultados" / "Verificacion_Energia"
OVERLEAF_DIR = SCRIPT_DIR.parent / "overleaf" / "imágenes"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
OVERLEAF_DIR.mkdir(parents=True, exist_ok=True)

# Paleta cromática científica de alto contraste (ColorBrewer / PRX style)
COLOR_BOSE = "#d90429"         # Rojo carmesí intenso (Teoría Bose)
COLOR_DIST = "#2b2b2b"         # Gris marengo oscuro (Teoría Distinguibles)
COLOR_PIMC_SWAP = "#0077b6"    # Azul zafiro / Teal (PIMC Bosones)
COLOR_PIMC_NOSWAP = "#e65100"  # Ámbar terracota (PIMC Distinguibles)
COLOR_POT = "#1b7837"          # Verde esmeralda (Energía Potencial)
COLOR_KIN = "#d95f02"          # Naranja cálido (Energía Cinética Virial)
COLOR_VIRIAL = "#7b1fa2"       # Púrpura intenso (Residuo Virial)
COLOR_RESIDUES = "#3182bd"     # Azul para residuos normalizados pull

OMEGA = 1.0


# ==============================================================================
# MODELOS ANALÍTICOS EXACTOS (CANÓNICO BOSÓNICO Y DISTINGUIBLES)
# ==============================================================================
def canonical_bose_partition(n_particles, temp, omega=OMEGA):
    """Calcula las sumas de partición canónicas exactas para bosones en oscilador 2D."""
    if n_particles < 1 or temp <= 0 or omega <= 0:
        raise ValueError("N, T y omega deben ser positivos")

    beta = 1.0 / temp
    partition = np.zeros(n_particles + 1)
    excitation_energy_sum = np.zeros(n_particles + 1)
    excited_partition = np.zeros(n_particles + 1)
    partition[0] = 1.0
    excited_partition[0] = 1.0

    for count in range(1, n_particles + 1):
        for cycle in range(1, count + 1):
            q = np.exp(-cycle * beta * omega)
            one_cycle_partition = 1.0 / (1.0 - q)**2
            one_cycle_excitation = 2.0 * omega * q / (1.0 - q)
            excited_cycle_partition = one_cycle_partition - 1.0
            partition[count] += one_cycle_partition * partition[count - cycle] / count
            excitation_energy_sum[count] += one_cycle_partition * (
                cycle * one_cycle_excitation * partition[count - cycle]
                + excitation_energy_sum[count - cycle]
            ) / count
            excited_partition[count] += (
                excited_cycle_partition * excited_partition[count - cycle] / count
            )

    return partition, excitation_energy_sum, excited_partition


def analytical_bose_energy(n_particles, temp, omega=OMEGA):
    """Energía canónica exacta por partícula y fracción de condensado N0/N."""
    partition, excitation_energy_sum, excited_partition = canonical_bose_partition(
        n_particles, temp, omega
    )
    e_per_particle = omega + excitation_energy_sum[n_particles] / (n_particles * partition[n_particles])
    ground_state_occupancy = sum(
        (n_particles - count) * excited_partition[count]
        for count in range(n_particles)
    )
    ground_state_occupancy /= partition[n_particles]
    return e_per_particle, ground_state_occupancy / n_particles


def analytical_distinguishable_energy(temp, omega=OMEGA):
    """Energía por partícula para partículas distinguibles en oscilador armónico 2D."""
    beta = 1.0 / temp
    return omega / np.tanh(0.5 * beta * omega)


def analytical_bose_cycle_fractions(n_particles, temp, omega=OMEGA):
    """Fracción canónica de partículas en ciclos de longitud L."""
    partition, _, _ = canonical_bose_partition(n_particles, temp, omega)
    beta = 1.0 / temp
    return np.array([
        (1.0 / (1.0 - np.exp(-cycle * beta * omega))**2)
        * partition[n_particles - cycle] / (n_particles * partition[n_particles])
        for cycle in range(1, n_particles + 1)
    ])


# ==============================================================================
# PROCESAMIENTO ESTADÍSTICO DE BLOQUES MONTE CARLO
# ==============================================================================
def blocked_ratio_stats(values, counts):
    """Calcula la media ponderada y error estándar asintótico para estimador de razón."""
    values = np.asarray(values, dtype=float)
    counts = np.asarray(counts, dtype=float)
    if len(values) < 2 or np.any(counts < 0) or counts.sum() <= 0:
        raise ValueError("Se necesitan al menos dos bloques válidos para la estadística")
    mean = np.average(values, weights=counts)
    residuals = counts * (values - mean)
    error = np.std(residuals, ddof=1) / (np.sqrt(len(residuals)) * counts.mean())
    return mean, error


def save_figure(fig, filename_base):
    """Guarda en PNG (300 DPI) y PDF vectorial tanto en Resultados como en Overleaf."""
    for folder in [RESULTS_DIR, OVERLEAF_DIR]:
        png_path = folder / f"{filename_base}.png"
        pdf_path = folder / f"{filename_base}.pdf"
        fig.savefig(png_path, dpi=300, bbox_inches="tight")
        fig.savefig(pdf_path, bbox_inches="tight")
    print(f"  [GUARDADO] {filename_base}.png/.pdf -> Resultados/ y overleaf/imágenes/")
    plt.close(fig)


# ==============================================================================
# CARGA Y PREPARACIÓN DE DATOS DE SIMULACIÓN
# ==============================================================================
def load_energy_blocks(filepath, burn_in_frac=0.20):
    """Carga outworm_energy_blocks.dat y descarta el calentamiento inicial."""
    data = np.loadtxt(filepath, skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    # data cols: block, samples, Etot, Epot, Ekin
    useful = data[int(len(data) * burn_in_frac):]
    valid = useful[:, 1] > 0
    return useful[valid]


def load_cycle_blocks(filepath, burn_in_frac=0.20):
    """Carga bloques de ciclos, descartando calentamiento y bloques sin muestras."""
    data = np.loadtxt(filepath, skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    useful = data[int(len(data) * burn_in_frac):]
    return useful[useful[:, 1] > 0]


def ensure_simulation_data(n_particles=5, temp=0.5, force_run=False):
    """Asegura que los datos de simulación existen para las comparativas."""
    if n_particles == 2:
        sim_swap = SCRIPT_DIR / "Resultados" / "Verificacion_Densidad" / "sim_swap"
        sim_noswap = SCRIPT_DIR / "Resultados" / "Verificacion_Densidad" / "sim_noswap"
        if sim_swap.exists() and sim_noswap.exists() and not force_run:
            return sim_swap, sim_noswap

    base_n = SCRIPT_DIR / "Resultados" / "NParticulas2D" / f"N{n_particles}"
    sim_swap = base_n / f"sim_T{temp:g}"
    sim_noswap = base_n / f"sim_noswap_T{temp:g}"

    # Si falta la simulación sin swap para N=5, la ejecutamos rápidamente
    if force_run or not (sim_noswap / "outworm_energy_blocks.dat").exists():
        print(f"-> Generando simulación sin Swap para N={n_particles}, T={temp}...")
        sim_noswap.mkdir(parents=True, exist_ok=True)
        mc_bin = SCRIPT_DIR / "MC"
        if not mc_bin.exists():
            raise RuntimeError("Ejecutable MC no encontrado. Compilar primero.")
        shutil.copy2(mc_bin, sim_noswap / "MC")

        with open(SCRIPT_DIR / "MC.cfg", "r") as f:
            cfg = f.read()
        cfg = re.sub(r"N=\s*\d+", f"N= {n_particles}", cfg)
        cfg = re.sub(r"T=\s*[\d\.]+", f"T= {temp}", cfg)
        cfg = re.sub(r"Npop=\s*\d+", "Npop= 30", cfg)
        cfg = re.sub(r"NwalkersMax=\s*\d+", "NwalkersMax= 60", cfg)
        cfg = re.sub(r"blck_heating=\s*\d+", "blck_heating= 20", cfg)
        cfg = re.sub(r"blck=\s*[\d\s]+", "blck= 80", cfg)
        cfg = re.sub(r"Niter=\s*\d+", "Niter= 500", cfg)
        cfg = re.sub(r"Nmeasure=\s*\d+", "Nmeasure= 1", cfg)

        with open(sim_noswap / "MC.cfg", "w") as f:
            f.write(cfg)

        env = os.environ.copy()
        env["WORM_DISABLE_SWAP"] = "1"
        subprocess.run(["./MC", "MC.cfg"], cwd=sim_noswap, env=env, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print("  -> Simulación sin Swap completada con éxito.")

    return sim_swap, sim_noswap


# ==============================================================================
# 1. FIGURA 1: DEMOSTRACIÓN CUÁNTICA - EFECTO SWAP (BOSE VS DISTINGUIBLE)
# ==============================================================================
def plot_01_demostracion_swap(n_particles=5, temp=0.5):
    """
    Demostración visual del colapso energético cuántico inducido por las permutaciones
    (movimiento Swap) frente al gas clásico sin permutaciones.
    """
    sim_swap, sim_noswap = ensure_simulation_data(n_particles, temp)

    data_swap = load_energy_blocks(sim_swap / "outworm_energy_blocks.dat")
    data_noswap = load_energy_blocks(sim_noswap / "outworm_energy_blocks.dat")

    blocks_s = data_swap[:, 0]
    weights_s = data_swap[:, 1]
    etot_s = data_swap[:, 2]

    blocks_ns = data_noswap[:, 0]
    weights_ns = data_noswap[:, 1]
    etot_ns = data_noswap[:, 2]

    mean_s, err_s = blocked_ratio_stats(etot_s, weights_s)
    mean_ns, err_ns = blocked_ratio_stats(etot_ns, weights_ns)

    e_bose, n0_frac = analytical_bose_energy(n_particles, temp)
    e_dist = analytical_distinguishable_energy(temp)

    delta_quantum = e_dist - e_bose
    delta_pimc = mean_ns - mean_s
    agreement_pct = (
        100.0 * (1.0 - abs(delta_pimc - delta_quantum) / abs(delta_quantum))
        if delta_quantum != 0 else np.nan
    )

    # -------------------------------------------------------------------------
    # Versión Completa: Panel Temporal + Histogramas Marginales
    # -------------------------------------------------------------------------
    fig = plt.figure(figsize=(9.2, 5.8), dpi=300)
    gs = fig.add_gridspec(1, 2, width_ratios=[3.4, 1.0], wspace=0.06)

    ax_main = fig.add_subplot(gs[0, 0])
    ax_hist = fig.add_subplot(gs[0, 1], sharey=ax_main)

    # 1. Series por bloques en ax_main
    ax_main.plot(blocks_ns, etot_ns, color=COLOR_PIMC_NOSWAP, alpha=0.35, linewidth=1.1,
                 linestyle="--", label="PIMC Sin Swap (Bloques)")
    ax_main.plot(blocks_s, etot_s, color=COLOR_PIMC_SWAP, alpha=0.35, linewidth=1.1,
                 linestyle="-", label="PIMC Con Swap (Bloques)")

    # Medias acumuladas móviles
    cum_ns = np.cumsum(etot_ns * weights_ns) / np.cumsum(weights_ns)
    cum_s = np.cumsum(etot_s * weights_s) / np.cumsum(weights_s)
    ax_main.plot(blocks_ns, cum_ns, color=COLOR_PIMC_NOSWAP, linewidth=2.2,
                 label=rf"PIMC Sin Swap: $\langle E \rangle/N = {mean_ns:.4f} \pm {err_ns:.4f}$")
    ax_main.plot(blocks_s, cum_s, color=COLOR_PIMC_SWAP, linewidth=2.2,
                 label=rf"PIMC Con Swap: $\langle E \rangle/N = {mean_s:.4f} \pm {err_s:.4f}$")

    # Líneas analíticas de referencia
    ax_main.axhline(e_dist, color=COLOR_DIST, linestyle=":", linewidth=2.0,
                    label=rf"Analítico Distinguible ($E_{{\mathrm{{dist}}}}/N = {e_dist:.4f}$)")
    ax_main.axhline(e_bose, color=COLOR_BOSE, linestyle="--", linewidth=2.0,
                    label=rf"Analítico Bose Exacto ($E_{{\mathrm{{Bose}}}}/N = {e_bose:.4f}$)")

    # Flecha y anotación del salto cuántico
    arrow_x = blocks_s[len(blocks_s) // 2]
    ax_main.annotate(
        "", xy=(arrow_x, e_bose), xytext=(arrow_x, e_dist),
        arrowprops=dict(arrowstyle="<->", color="#5c0099", lw=2.2, shrinkA=0, shrinkB=0)
    )
    ax_main.text(
        arrow_x + 1.5, (e_dist + e_bose) / 2.0,
        rf"$\Delta E_{{\mathrm{{cu\acute{{a}}ntica}}}} = {delta_quantum:.4f}$" + "\n" +
        rf"Caída: ${(delta_quantum / e_dist) * 100:.1f}\%$",
        color="#5c0099", fontsize=9.8, fontweight="bold", va="center"
    )

    ax_main.set_xlabel("Bloque de Medición Monte Carlo", fontsize=11)
    ax_main.set_ylabel(r"Energía por partícula $\langle E \rangle / N$ [$\hbar\omega$]", fontsize=11)
    ax_main.set_title(
        rf"Demostración Cuántica: Efecto Swap en Oscilador 2D ($N={n_particles}$, $T={temp}$)",
        fontsize=12, fontweight="bold", pad=10
    )
    ax_main.grid(True, linestyle=":", alpha=0.55)
    ax_main.legend(loc="upper left", fontsize=8.8, framealpha=0.92)

    # 2. Histogramas marginales en ax_hist
    bins = np.linspace(min(etot_s.min(), e_bose - 0.1), max(etot_ns.max(), e_dist + 0.1), 28)
    ax_hist.hist(etot_ns, bins=bins, orientation="horizontal", color=COLOR_PIMC_NOSWAP,
                 alpha=0.45, density=True, edgecolor=COLOR_PIMC_NOSWAP, linewidth=1.0)
    ax_hist.hist(etot_s, bins=bins, orientation="horizontal", color=COLOR_PIMC_SWAP,
                 alpha=0.45, density=True, edgecolor=COLOR_PIMC_SWAP, linewidth=1.0)

    ax_hist.axhline(e_dist, color=COLOR_DIST, linestyle=":", linewidth=2.0)
    ax_hist.axhline(e_bose, color=COLOR_BOSE, linestyle="--", linewidth=2.0)
    ax_hist.set_xlabel("Densidad", fontsize=10)
    ax_hist.tick_params(axis="y", labelleft=False)
    ax_hist.grid(True, linestyle=":", alpha=0.55)

    # Cuadro explicativo
    box_text = (
        r"$\mathbf{Demostraci\acute{o}n\ Bose:}$" + "\n" +
        rf"$\bullet$ Sin Swap: {mean_ns:.4f} ($\approx$ Boltzmann)" + "\n" +
        rf"$\bullet$ Con Swap: {mean_s:.4f} ($\approx$ Bose)" + "\n" +
        rf"$\bullet$ $\Delta E_{{\mathrm{{PIMC}}}} = {delta_pimc:.4f}$" + "\n" +
        rf"$\bullet$ Acuerdo teórico: {agreement_pct:.1f}\%"
    )
    ax_main.annotate(
        box_text, xy=(0.03, 0.05), xycoords="axes fraction",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#ffffff", edgecolor="#adb5bd", alpha=0.95),
        fontsize=9.0
    )

    save_figure(fig, "01_energia_demostracion_swap_bose_vs_distinguible")

    # -------------------------------------------------------------------------
    # Versión Simple (Panel Único para incrustación compacta)
    # -------------------------------------------------------------------------
    fig_s, ax_s = plt.subplots(figsize=(8.0, 5.2), dpi=300)
    ax_s.plot(blocks_ns, etot_ns, color=COLOR_PIMC_NOSWAP, alpha=0.25, linewidth=1.0, linestyle="--")
    ax_s.plot(blocks_s, etot_s, color=COLOR_PIMC_SWAP, alpha=0.25, linewidth=1.0, linestyle="-")
    ax_s.plot(blocks_ns, cum_ns, color=COLOR_PIMC_NOSWAP, linewidth=2.3,
              label=rf"PIMC Sin Swap (Distinguibles): $\langle E \rangle/N = {mean_ns:.4f} \pm {err_ns:.4f}$")
    ax_s.plot(blocks_s, cum_s, color=COLOR_PIMC_SWAP, linewidth=2.3,
              label=rf"PIMC Con Swap (Bosones): $\langle E \rangle/N = {mean_s:.4f} \pm {err_s:.4f}$")
    ax_s.axhline(e_dist, color=COLOR_DIST, linestyle=":", linewidth=2.0, label=rf"Distinguible ({e_dist:.4f})")
    ax_s.axhline(e_bose, color=COLOR_BOSE, linestyle="--", linewidth=2.0, label=rf"Bose Exacto ({e_bose:.4f})")
    ax_s.set_xlabel("Bloque de Medición Monte Carlo", fontsize=11)
    ax_s.set_ylabel(r"Energía por partícula $\langle E \rangle / N$ [$\hbar\omega$]", fontsize=11)
    ax_s.set_title(rf"Demostración del Salto Cuántico Bosónico por Permutaciones ($N={n_particles}$, $T={temp}$)",
                   fontsize=12, fontweight="bold", pad=10)
    ax_s.grid(True, linestyle=":", alpha=0.55)
    ax_s.legend(loc="upper left", fontsize=9.0, framealpha=0.92)
    save_figure(fig_s, "01_energia_demostracion_swap_bose_vs_distinguible_simple")


# ==============================================================================
# 2. FIGURA 2: COMPROBACIÓN FÍSICA - TEOREMA DEL VIRIAL (<Epot> = <Ekin>)
# ==============================================================================
def plot_02_comprobacion_virial(n_particles=5, temp=0.5):
    """
    Comprueba la consistencia del estimador virial centroidal por ciclos con
    Epot - Ekin = 0 dentro de la incertidumbre estadística.
    """
    sim_swap, _ = ensure_simulation_data(n_particles, temp)
    data = load_energy_blocks(sim_swap / "outworm_energy_blocks.dat")

    blocks = data[:, 0]
    weights = data[:, 1]
    epot = data[:, 3]
    ekin = data[:, 4]
    diff = epot - ekin

    mean_pot, err_pot = blocked_ratio_stats(epot, weights)
    mean_kin, err_kin = blocked_ratio_stats(ekin, weights)
    mean_diff, err_diff = blocked_ratio_stats(diff, weights)

    z_virial = abs(mean_diff) / err_diff if err_diff > 0 else 0.0
    virial_ok = z_virial <= 3.0
    e_bose, _ = analytical_bose_energy(n_particles, temp)
    half_e = e_bose / 2.0

    # -------------------------------------------------------------------------
    # Versión Completa: Panel 1 (Epot, Ekin) + Panel 2 (Diferencia Virial)
    # -------------------------------------------------------------------------
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(8.2, 7.0), sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.1], "hspace": 0.08}, dpi=300
    )

    # Panel Superior: Epot vs Ekin
    ax_top.plot(blocks, epot, color=COLOR_POT, alpha=0.40, linewidth=1.1, label=r"$E_{{\mathrm{{pot}}}}$ (Bloques)")
    ax_top.plot(blocks, ekin, color=COLOR_KIN, alpha=0.40, linewidth=1.1, linestyle="--", label=r"$E_{{\mathrm{{kin}}}}$ Virial (Bloques)")

    cum_pot = np.cumsum(epot * weights) / np.cumsum(weights)
    cum_kin = np.cumsum(ekin * weights) / np.cumsum(weights)

    ax_top.plot(blocks, cum_pot, color=COLOR_POT, linewidth=2.2,
                label=rf"$\langle E_{{\mathrm{{pot}}}} \rangle = {mean_pot:.4f} \pm {err_pot:.4f}$")
    ax_top.plot(blocks, cum_kin, color=COLOR_KIN, linewidth=2.2, linestyle="--",
                label=rf"$\langle E_{{\mathrm{{kin}}}} \rangle = {mean_kin:.4f} \pm {err_kin:.4f}$")
    ax_top.axhline(half_e, color=COLOR_BOSE, linestyle=":", linewidth=1.8,
                   label=rf"Teoría Virial: $E_{{\mathrm{{Bose}}}} / 2 = {half_e:.4f}$")

    ax_top.set_ylabel(r"Energía por partícula [$\hbar\omega$]", fontsize=11)
    ax_top.set_title(
        rf"Comprobación del Teorema del Virial Cuántico ($\langle E_{{\mathrm{{pot}}}} \rangle = \langle E_{{\mathrm{{kin}}}} \rangle$)",
        fontsize=12, fontweight="bold", pad=10
    )
    ax_top.grid(True, linestyle=":", alpha=0.55)
    ax_top.legend(loc="upper right", fontsize=9.0, framealpha=0.92)

    # Panel Inferior: Residuo Virial Epot - Ekin
    ax_bot.axhspan(-3.0 * err_diff, 3.0 * err_diff, color="#d4edda", alpha=0.65,
                   label=rf"Banda de confianza $3\sigma$ ($\pm {3.0*err_diff:.4f}$)")
    ax_bot.axhspan(-1.0 * err_diff, 1.0 * err_diff, color="#c3e6cb", alpha=0.50,
                   label=rf"Banda de confianza $1\sigma$ ($\pm {err_diff:.4f}$)")
    ax_bot.axhline(0.0, color="black", linestyle="-", linewidth=1.0)

    cum_diff = cum_pot - cum_kin
    ax_bot.plot(blocks, diff, color=COLOR_VIRIAL, alpha=0.35, linewidth=1.0, linestyle="none", marker="o", markersize=3.0)
    ax_bot.plot(blocks, cum_diff, color=COLOR_VIRIAL, linewidth=2.2,
                label=rf"$\Delta_{{\mathrm{{vir}}}} = {mean_diff:+.4f} \pm {err_diff:.4f}$")

    ax_bot.set_xlabel("Bloque de Medición Monte Carlo", fontsize=11)
    ax_bot.set_ylabel(r"$E_{{\mathrm{{pot}}}} - E_{{\mathrm{{kin}}}}$", fontsize=10.5)
    ax_bot.grid(True, linestyle=":", alpha=0.55)
    ax_bot.legend(loc="lower right", fontsize=8.5, framealpha=0.92)

    # Cuadro de verificación del virial
    virial_box = (
        r"$\mathbf{Validaci\acute{o}n\ del\ Virial}:$" + "\n" +
        rf"$|\langle E_{{\mathrm{{pot}}}} \rangle - \langle E_{{\mathrm{{kin}}}} \rangle| = {abs(mean_diff):.4f}$" + "\n" +
        rf"Pull virial: $|z_{{\mathrm{{vir}}}}| = {z_virial:.2f}\sigma$" + "\n" +
        (r"Estado: $\mathbf{\checkmark\ compatible\ con\ 3\sigma}$" if virial_ok
         else r"Estado: $\mathbf{\times\ fuera\ de\ 3\sigma}$")
    )
    ax_top.annotate(
        virial_box, xy=(0.03, 0.06), xycoords="axes fraction",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#ffffff",
              edgecolor="#28a745" if virial_ok else "#dc3545", alpha=0.95),
        fontsize=9.0
    )

    save_figure(fig, "02_energia_comprobacion_teorema_virial")

    # -------------------------------------------------------------------------
    # Versión Simple: Solo panel de Epot y Ekin
    # -------------------------------------------------------------------------
    fig_s, ax_s = plt.subplots(figsize=(8.0, 5.2), dpi=300)
    ax_s.plot(blocks, cum_pot, color=COLOR_POT, linewidth=2.3,
              label=rf"$\langle E_{{\mathrm{{pot}}}} \rangle = {mean_pot:.4f} \pm {err_pot:.4f}$")
    ax_s.plot(blocks, cum_kin, color=COLOR_KIN, linewidth=2.3, linestyle="--",
              label=rf"$\langle E_{{\mathrm{{kin}}}} \rangle = {mean_kin:.4f} \pm {err_kin:.4f}$")
    ax_s.axhline(half_e, color=COLOR_BOSE, linestyle=":", linewidth=2.0,
                 label=rf"Teoría: $E_{{\mathrm{{Bose}}}} / 2 = {half_e:.4f}$")
    ax_s.set_xlabel("Bloque de Medición Monte Carlo", fontsize=11)
    ax_s.set_ylabel(r"Energía por partícula [$\hbar\omega$]", fontsize=11)
    ax_s.set_title(rf"Equipartición Virial: $\langle E_{{\mathrm{{pot}}}} \rangle$ vs $\langle E_{{\mathrm{{kin}}}} \rangle$ ($N={n_particles}$, $T={temp}$)",
                   fontsize=12, fontweight="bold", pad=10)
    ax_s.grid(True, linestyle=":", alpha=0.55)
    ax_s.legend(loc="upper right", fontsize=9.2, framealpha=0.92)
    save_figure(fig_s, "02_energia_comprobacion_teorema_virial_simple")


# ==============================================================================
# 3. FIGURA 3: VERIFICACIÓN ESTADÍSTICA - CONVERGENCIA MCMC Y PULL PLOT
# ==============================================================================
def plot_03_verificacion_convergencia_mcmc(n_particles=5, temp=0.5):
    """
    Verificación rigurosa de convergencia estocástica en Sector Z:
      - Media acumulada ponderada con banda de error 1-sigma y 3-sigma.
      - Panel inferior con pull estandarizado z(b) = (E_cum - E_exact) / sigma_cum.
    """
    sim_swap, _ = ensure_simulation_data(n_particles, temp)
    data = load_energy_blocks(sim_swap / "outworm_energy_blocks.dat")

    blocks = data[:, 0]
    weights = data[:, 1]
    etot = data[:, 2]

    e_exact, n0_frac = analytical_bose_energy(n_particles, temp)
    mean_tot, err_tot = blocked_ratio_stats(etot, weights)

    # Calcular medias y errores acumulados bloque a bloque
    n_blocks = len(etot)
    cum_means = np.zeros(n_blocks)
    cum_errors = np.zeros(n_blocks)

    for k in range(1, n_blocks + 1):
        sub_v = etot[:k]
        sub_w = weights[:k]
        m = np.average(sub_v, weights=sub_w)
        cum_means[k - 1] = m
        if k >= 3:
            res = sub_w * (sub_v - m)
            cum_errors[k - 1] = np.std(res, ddof=1) / (np.sqrt(k) * sub_w.mean())
        else:
            cum_errors[k - 1] = err_tot

    z_scores = (cum_means - e_exact) / np.maximum(cum_errors, 1e-6)
    diff_pct = abs(mean_tot - e_exact) / e_exact * 100.0
    discrepancy_ok = diff_pct < 1.0
    pull_ok = abs(z_scores[-1]) <= 3.0

    # -------------------------------------------------------------------------
    # Versión Completa: Panel 1 (Convergencia) + Panel 2 (Pull Normalizado)
    # -------------------------------------------------------------------------
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(8.2, 7.2), sharex=True,
        gridspec_kw={"height_ratios": [2.3, 1.0], "hspace": 0.08}, dpi=300
    )

    # Panel Superior: Convergencia de E
    ax_top.plot(blocks, etot, color="#6c757d", alpha=0.35, linewidth=1.0, label="E(bloque) PIMC")

    # Bandas de confianza sobre la media acumulada
    ax_top.fill_between(blocks, cum_means - 3.0 * cum_errors, cum_means + 3.0 * cum_errors,
                        color="#a8dadc", alpha=0.35, label=r"Banda de confianza acumulada $3\sigma$ (99.7%)")
    ax_top.fill_between(blocks, cum_means - 1.0 * cum_errors, cum_means + 1.0 * cum_errors,
                        color="#457b9d", alpha=0.30, label=r"Banda de confianza acumulada $1\sigma$ (68.3%)")

    ax_top.plot(blocks, cum_means, color=COLOR_PIMC_SWAP, linewidth=2.3,
                label=rf"Media acumulada: $\langle E \rangle / N = {mean_tot:.4f} \pm {err_tot:.4f}$")
    ax_top.axhline(e_exact, color=COLOR_BOSE, linestyle="--", linewidth=2.0,
                   label=rf"Analítico Bose Exacto ($E_{{\mathrm{{Bose}}}}/N = {e_exact:.4f}$)")

    ax_top.set_ylabel(r"Energía por partícula $\langle E \rangle / N$ [$\hbar\omega$]", fontsize=11)
    ax_top.set_title(
        rf"Verificación de Convergencia MCMC y Regla 3-Sigma ($N={n_particles}$, $T={temp}$)",
        fontsize=12, fontweight="bold", pad=10
    )
    ax_top.grid(True, linestyle=":", alpha=0.55)
    ax_top.legend(loc="upper right", fontsize=8.8, framealpha=0.92)

    # Panel Inferior: Pull Plot
    ax_bot.axhspan(-3.0, 3.0, color="#d4edda", alpha=0.65, label=r"Zona permitida $\pm 3\sigma$ (99.73%)")
    ax_bot.axhspan(-1.0, 1.0, color="#c3e6cb", alpha=0.50, label=r"Zona $\pm 1\sigma$ (68.3%)")
    ax_bot.axhline(0.0, color="black", linestyle="-", linewidth=0.9)
    ax_bot.axhline(+3.0, color="#dc3545", linestyle=":", linewidth=1.1)
    ax_bot.axhline(-3.0, color="#dc3545", linestyle=":", linewidth=1.1)

    ax_bot.plot(blocks[2:], z_scores[2:], color=COLOR_PIMC_SWAP, linewidth=1.8, marker="o", markersize=3.2)
    ax_bot.set_ylim(-3.8, 3.8)
    ax_bot.set_xlabel("Bloque de Medición Monte Carlo", fontsize=11)
    ax_bot.set_ylabel(r"Pull $z = \frac{\bar{E} - E_{\mathrm{exact}}}{\sigma_{\bar{E}}}$", fontsize=10.5)
    ax_bot.grid(True, linestyle=":", alpha=0.55)
    ax_bot.legend(loc="lower right", fontsize=8.2, framealpha=0.92)

    # Cuadro explicativo
    stat_box = (
        r"$\mathbf{Validaci\acute{o}n\ Num\acute{e}rica}:$" + "\n" +
        rf"$\bullet$ Discrepancia: {diff_pct:.2f}\%" + "\n" +
        rf"$\bullet$ Pull final: $z = {z_scores[-1]:+.2f}\sigma$" + "\n" +
        (r"$\bullet$ Discrepancia $<1\%$: $\mathbf{\checkmark}$" if discrepancy_ok
         else r"$\bullet$ Discrepancia $<1\%$: $\mathbf{\times}$") + "\n" +
        (r"$\bullet$ Pull $|z|\leq3$: $\mathbf{\checkmark}$" if pull_ok
         else r"$\bullet$ Pull $|z|\leq3$: $\mathbf{\times}$")
    )
    ax_top.annotate(
        stat_box, xy=(0.03, 0.06), xycoords="axes fraction",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#ffffff",
              edgecolor="#28a745" if discrepancy_ok and pull_ok else "#dc3545", alpha=0.95),
        fontsize=9.0
    )

    save_figure(fig, "03_energia_verificacion_convergencia_mcmc")

    # -------------------------------------------------------------------------
    # Versión Simple: Solo panel de convergencia
    # -------------------------------------------------------------------------
    fig_s, ax_s = plt.subplots(figsize=(8.0, 5.2), dpi=300)
    ax_s.plot(blocks, etot, color="#6c757d", alpha=0.30, linewidth=1.0)
    ax_s.fill_between(blocks, cum_means - 3.0 * cum_errors, cum_means + 3.0 * cum_errors,
                      color="#a8dadc", alpha=0.35, label=r"Banda de confianza $3\sigma$")
    ax_s.plot(blocks, cum_means, color=COLOR_PIMC_SWAP, linewidth=2.3,
              label=rf"Media acumulada: {mean_tot:.4f} $\pm$ {err_tot:.4f}")
    ax_s.axhline(e_exact, color=COLOR_BOSE, linestyle="--", linewidth=2.0,
                 label=rf"Analítico Bose ({e_exact:.4f})")
    ax_s.set_xlabel("Bloque de Medición Monte Carlo", fontsize=11)
    ax_s.set_ylabel(r"Energía por partícula $\langle E \rangle / N$ [$\hbar\omega$]", fontsize=11)
    ax_s.set_title(rf"Convergencia de Energía en Sector Z ($N={n_particles}$, $T={temp}$)",
                   fontsize=12, fontweight="bold", pad=10)
    ax_s.grid(True, linestyle=":", alpha=0.55)
    ax_s.legend(loc="upper right", fontsize=9.2, framealpha=0.92)
    save_figure(fig_s, "03_energia_verificacion_convergencia_mcmc_simple")


# ==============================================================================
# 4. FIGURA 4: VERIFICACIÓN TERMODINÁMICA - CURVA E(T) Y CONDENSADO N0/N
# ==============================================================================
def plot_04_termodinamica_temperatura(n_particles=5):
    """
    Curva completa E(T) vs Temperatura con:
      - Comparativa Bose Analítico vs Distinguible Analítico vs PIMC Worm.
    - Panel inferior con ocupación del estado fundamental N0/N y ciclos no triviales L >= 2.
    """
    temps = [0.2, 0.5, 1.0, 1.5, 2.0, 3.0]
    base_dir = SCRIPT_DIR / "Resultados" / "NParticulas2D" / f"N{n_particles}"
    if not base_dir.exists():
        base_dir = SCRIPT_DIR / "Resultados" / "NParticulas2D" / "N5"
        n_particles = 5

    pimc_e = []
    pimc_err = []
    pimc_permuted_cycles = []
    pimc_permuted_err = []

    for t in temps:
        sim_dir = base_dir / f"sim_T{t:g}"
        f_e = sim_dir / "outworm_energy_blocks.dat"
        data_e = load_energy_blocks(f_e)
        m, e = blocked_ratio_stats(data_e[:, 2], data_e[:, 1])
        pimc_e.append(m)
        pimc_err.append(e)

        f_c = sim_dir / "outworm_cycle_blocks.dat"
        if f_c.exists():
            u_c = load_cycle_blocks(f_c)
            w_c = u_c[:, 1]
            macro_frac = np.sum(u_c[:, 3:], axis=1)  # L >= 2
            m_c, e_c = blocked_ratio_stats(macro_frac, w_c)
            pimc_permuted_cycles.append(m_c)
            pimc_permuted_err.append(e_c)
        else:
            pimc_permuted_cycles.append(0.0)
            pimc_permuted_err.append(0.0)

    # Curvas analíticas continuas
    t_fine = np.linspace(0.15, 3.2, 200)
    e_dist_fine = [analytical_distinguishable_energy(t) for t in t_fine]
    bose_res = [analytical_bose_energy(n_particles, t) for t in t_fine]
    e_bose_fine = [r[0] for r in bose_res]
    n0_fine = [r[1] for r in bose_res]

    # -------------------------------------------------------------------------
    # Versión Completa: Panel 1 (E(T)) + Panel 2 (Ocupación fundamental y ciclos)
    # -------------------------------------------------------------------------
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(8.4, 7.4), sharex=True,
        gridspec_kw={"height_ratios": [2.3, 1.1], "hspace": 0.08}, dpi=300
    )

    # Panel Superior: E(T)
    ax_top.plot(t_fine, e_dist_fine, color=COLOR_DIST, linestyle="-", linewidth=2.2,
                label="Analítico Distinguibles (Boltzmannons)")
    ax_top.plot(t_fine, e_bose_fine, color=COLOR_BOSE, linestyle="--", linewidth=2.2,
                label=rf"Analítico Bosones Canónico ($N={n_particles}$)")

    ax_top.errorbar(temps, pimc_e, yerr=pimc_err, fmt="o", color=COLOR_PIMC_SWAP,
                    ecolor="black", elinewidth=1.2, capsize=3.5, markersize=5.8,
                    label=rf"PIMC Worm Algorithm ($N={n_particles}$)", zorder=5)

    # Anotaciones físicas del límite fundamental y clásico
    ax_top.annotate(
        r"Límite Fundamental ($T \to 0$):" + "\n" + r"$E/N \to 1.0\ \hbar\omega$",
        xy=(0.2, 1.0), xytext=(0.43, 0.08), textcoords="axes fraction",
        arrowprops=dict(arrowstyle="->", color=COLOR_BOSE, lw=1.3),
        fontsize=8.8, bbox=dict(boxstyle="round,pad=0.4", facecolor="#fff", edgecolor="#ccc", alpha=0.9)
    )
    ax_top.annotate(
        r"Límite Clásico ($T \gg 1$):" + "\n" + r"$E/N \to 2 k_B T$",
        xy=(3.0, analytical_distinguishable_energy(3.0)), xytext=(2.1, 4.9),
        arrowprops=dict(arrowstyle="->", color=COLOR_DIST, lw=1.3),
        fontsize=8.8, bbox=dict(boxstyle="round,pad=0.4", facecolor="#fff", edgecolor="#ccc", alpha=0.9)
    )

    ax_top.set_ylabel(r"Energía por partícula $\langle E \rangle / N$ [$\hbar\omega$]", fontsize=11)
    ax_top.set_title(
        rf"Termodinámica Energética: Curva $E(T)$ en Oscilador Armónico 2D ($N={n_particles}$)",
        fontsize=12, fontweight="bold", pad=10
    )
    ax_top.grid(True, linestyle=":", alpha=0.55)
    ax_top.legend(loc="upper left", fontsize=9.0, framealpha=0.92)

    # Panel Inferior: Ocupación fundamental y ciclos no triviales
    ax_bot.plot(t_fine, n0_fine, color=COLOR_BOSE, linewidth=2.0,
                label=r"Ocupación del Estado Fundamental $N_0/N$ (Bose)")
    ax_bot.errorbar(temps, pimc_permuted_cycles, yerr=pimc_permuted_err, fmt="s", color="#5c0099",
                    ecolor="black", elinewidth=1.1, capsize=3.0, markersize=5.0,
                    label=r"Fracción PIMC en Ciclos No Triviales ($L \geq 2$)", zorder=4)

    ax_bot.set_ylim(-0.05, 1.05)
    ax_bot.set_xlabel(r"Temperatura reducida $T$ [$k_B / \hbar\omega$]", fontsize=11)
    ax_bot.set_ylabel("Fracción de partículas", fontsize=10.5)
    ax_bot.grid(True, linestyle=":", alpha=0.55)
    ax_bot.legend(loc="upper right", fontsize=8.8, framealpha=0.92)

    save_figure(fig, "04_energia_termodinamica_vs_temperatura")

    # -------------------------------------------------------------------------
    # Versión Simple: Solo curva E(T)
    # -------------------------------------------------------------------------
    fig_s, ax_s = plt.subplots(figsize=(8.0, 5.2), dpi=300)
    ax_s.plot(t_fine, e_dist_fine, color=COLOR_DIST, linestyle="-", linewidth=2.2, label="Distinguibles (Boltzmann)")
    ax_s.plot(t_fine, e_bose_fine, color=COLOR_BOSE, linestyle="--", linewidth=2.2, label=rf"Bosones Canónico ($N={n_particles}$)")
    ax_s.errorbar(temps, pimc_e, yerr=pimc_err, fmt="o", color=COLOR_PIMC_SWAP, ecolor="black",
                  elinewidth=1.2, capsize=3.5, markersize=6.0, label="PIMC Worm", zorder=5)
    ax_s.set_xlabel(r"Temperatura reducida $T$ [$k_B / \hbar\omega$]", fontsize=11)
    ax_s.set_ylabel(r"Energía por partícula $\langle E \rangle / N$ [$\hbar\omega$]", fontsize=11)
    ax_s.set_title(rf"Energía vs Temperatura en Oscilador Armónico 2D ($N={n_particles}$)",
                   fontsize=12, fontweight="bold", pad=10)
    ax_s.grid(True, linestyle=":", alpha=0.55)
    ax_s.legend(loc="upper left", fontsize=9.2, framealpha=0.92)
    save_figure(fig_s, "04_energia_termodinamica_vs_temperatura_simple")


# ==============================================================================
# 5. FIGURA 5: PANEL COMPLETO MULTICOMPONENTE (DASHBOARD 2x2)
# ==============================================================================
def plot_05_panel_completo_dashboard(n_particles=5, temp=0.5):
    """
    Panel panorámico 2x2 que reúne las 4 demostraciones y comprobaciones clave:
      (a) Efecto cuántico Swap (Caída bosónica)
      (b) Teorema del Virial (<Epot> = <Ekin>)
      (c) Convergencia MCMC en Sector Z
      (d) Curva termodinámica E(T)
    """
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.5), dpi=300)
    ax_a, ax_b = axes[0, 0], axes[0, 1]
    ax_c, ax_d = axes[1, 0], axes[1, 1]

    sim_swap, sim_noswap = ensure_simulation_data(n_particles, temp)
    data_s = load_energy_blocks(sim_swap / "outworm_energy_blocks.dat")
    data_ns = load_energy_blocks(sim_noswap / "outworm_energy_blocks.dat")

    # (a) Efecto Swap
    b_s, w_s, e_s = data_s[:, 0], data_s[:, 1], data_s[:, 2]
    b_ns, w_ns, e_ns = data_ns[:, 0], data_ns[:, 1], data_ns[:, 2]
    cum_s = np.cumsum(e_s * w_s) / np.cumsum(w_s)
    cum_ns = np.cumsum(e_ns * w_ns) / np.cumsum(w_ns)
    e_bose, _ = analytical_bose_energy(n_particles, temp)
    e_dist = analytical_distinguishable_energy(temp)

    ax_a.plot(b_ns, cum_ns, color=COLOR_PIMC_NOSWAP, linewidth=2.0, label=r"PIMC Sin Swap")
    ax_a.plot(b_s, cum_s, color=COLOR_PIMC_SWAP, linewidth=2.0, label=r"PIMC Con Swap")
    ax_a.axhline(e_dist, color=COLOR_DIST, linestyle=":", linewidth=1.8, label=rf"Distinguible ({e_dist:.3f})")
    ax_a.axhline(e_bose, color=COLOR_BOSE, linestyle="--", linewidth=1.8, label=rf"Bose ({e_bose:.3f})")
    ax_a.set_title(r"(a) Demostración del Salto Cuántico (Efecto Swap)", fontsize=11, fontweight="bold")
    ax_a.set_xlabel("Bloques MC", fontsize=10)
    ax_a.set_ylabel(r"$\langle E \rangle / N$ [$\hbar\omega$]", fontsize=10)
    ax_a.grid(True, linestyle=":", alpha=0.55)
    ax_a.legend(loc="upper left", fontsize=8.2, framealpha=0.9)

    # (b) Teorema del Virial
    epot, ekin = data_s[:, 3], data_s[:, 4]
    cum_pot = np.cumsum(epot * w_s) / np.cumsum(w_s)
    cum_kin = np.cumsum(ekin * w_s) / np.cumsum(w_s)
    half_e = e_bose / 2.0
    ax_b.plot(b_s, cum_pot, color=COLOR_POT, linewidth=2.0, label=rf"$\langle E_{{\mathrm{{pot}}}} \rangle$")
    ax_b.plot(b_s, cum_kin, color=COLOR_KIN, linewidth=2.0, linestyle="--", label=rf"$\langle E_{{\mathrm{{kin}}}} \rangle$ Virial")
    ax_b.axhline(half_e, color=COLOR_BOSE, linestyle=":", linewidth=1.8, label=rf"Virial Exacto ({half_e:.3f})")
    ax_b.set_title(r"(b) Comprobación del Teorema del Virial", fontsize=11, fontweight="bold")
    ax_b.set_xlabel("Bloques MC", fontsize=10)
    ax_b.set_ylabel(r"Componentes [$\hbar\omega$]", fontsize=10)
    ax_b.grid(True, linestyle=":", alpha=0.55)
    ax_b.legend(loc="upper right", fontsize=8.2, framealpha=0.9)

    # (c) Convergencia MCMC
    mean_s, err_s = blocked_ratio_stats(e_s, w_s)
    ax_c.plot(b_s, e_s, color="#888", alpha=0.25, linewidth=0.8)
    ax_c.plot(b_s, cum_s, color=COLOR_PIMC_SWAP, linewidth=2.0, label=rf"Media PIMC: {mean_s:.4f} $\pm$ {err_s:.4f}")
    ax_c.axhline(e_bose, color=COLOR_BOSE, linestyle="--", linewidth=1.8, label=rf"Bose ({e_bose:.4f})")
    ax_c.set_title(r"(c) Verificación de Convergencia MCMC (Sector Z)", fontsize=11, fontweight="bold")
    ax_c.set_xlabel("Bloques MC", fontsize=10)
    ax_c.set_ylabel(r"$\langle E \rangle / N$ [$\hbar\omega$]", fontsize=10)
    ax_c.grid(True, linestyle=":", alpha=0.55)
    ax_c.legend(loc="upper right", fontsize=8.2, framealpha=0.9)

    # (d) Termodinámica E(T)
    temps = [0.2, 0.5, 1.0, 1.5, 2.0, 3.0]
    base_dir = SCRIPT_DIR / "Resultados" / "NParticulas2D" / f"N{n_particles}"
    if not base_dir.exists():
        base_dir = SCRIPT_DIR / "Resultados" / "NParticulas2D" / "N5"
        n_particles = 5
    p_e, p_err = [], []
    for t in temps:
        d = load_energy_blocks(base_dir / f"sim_T{t:g}" / "outworm_energy_blocks.dat")
        m, e = blocked_ratio_stats(d[:, 2], d[:, 1])
        p_e.append(m)
        p_err.append(e)

    t_f = np.linspace(0.15, 3.2, 150)
    ax_d.plot(t_f, [analytical_distinguishable_energy(t) for t in t_f], color=COLOR_DIST, linewidth=1.8, label="Distinguible")
    ax_d.plot(t_f, [analytical_bose_energy(n_particles, t)[0] for t in t_f], color=COLOR_BOSE, linestyle="--", linewidth=1.8, label="Bose Canónico")
    ax_d.errorbar(temps, p_e, yerr=p_err, fmt="o", color=COLOR_PIMC_SWAP, ecolor="black", capsize=3, markersize=5, label="PIMC Worm")
    ax_d.set_title(r"(d) Verificación Térmica $E(T)$ vs Temperatura", fontsize=11, fontweight="bold")
    ax_d.set_xlabel(r"Temperatura $T$ [$k_B/\hbar\omega$]", fontsize=10)
    ax_d.set_ylabel(r"$\langle E \rangle / N$ [$\hbar\omega$]", fontsize=10)
    ax_d.grid(True, linestyle=":", alpha=0.55)
    ax_d.legend(loc="upper left", fontsize=8.2, framealpha=0.9)

    fig.tight_layout(pad=2.2)
    save_figure(fig, "05_energia_panel_completo_dashboard")


# ==============================================================================
# IMPRESIÓN DE INFORME Y TABLAS CUANTITATIVAS
# ==============================================================================
def print_verification_report(n_particles=5):
    """Imprime por consola la tabla completa de verificación energética."""
    temps = [0.2, 0.5, 1.0, 1.5, 2.0, 3.0]
    base_dir = SCRIPT_DIR / "Resultados" / "NParticulas2D" / f"N{n_particles}"
    if not base_dir.exists():
        base_dir = SCRIPT_DIR / "Resultados" / "NParticulas2D" / "N5"
        n_particles = 5

    print("\n" + "=" * 122)
    print(f"  INFORME CUANTITATIVO DE VERIFICACIÓN Y COMPROBACIÓN ENERGÉTICA (N = {n_particles})")
    print("=" * 122)
    print(f"{'T':>5} | {'E/N PIMC':>17} | {'E/N Bose':>10} | {'E/N Dist':>10} | {'Discrep.':>9} | {'<1%':>6} | {'Pull z':>12} | {'Epot - Ekin':>21} | {'|zvir|':>5} | {'3sigma':>6}")
    print("-" * 122)

    latex_rows = []

    for t in temps:
        d = load_energy_blocks(base_dir / f"sim_T{t:g}" / "outworm_energy_blocks.dat")
        w = d[:, 1]
        etot, err_tot = blocked_ratio_stats(d[:, 2], w)
        diff, err_diff = blocked_ratio_stats(d[:, 3] - d[:, 4], w)

        e_bose, n0 = analytical_bose_energy(n_particles, t)
        e_dist = analytical_distinguishable_energy(t)
        disc_pct = abs(etot - e_bose) / e_bose * 100.0
        z_score = (etot - e_bose) / err_tot if err_tot > 0 else 0.0
        z_vir = abs(diff) / err_diff if err_diff > 0 else 0.0
        discrepancy_ok = "OK" if disc_pct < 1.0 else "FUERA"
        vir_ok = "OK" if z_vir <= 3.0 else "FALLO"

        print(f"{t:5.2f} | {etot:8.4f} +/- {err_tot:.4f} | {e_bose:10.4f} | {e_dist:10.4f} | {disc_pct:8.2f}% | {discrepancy_ok:>6} | {z_score:+7.2f} sigma | {diff:+8.4f} +/- {err_diff:.4f} | {z_vir:5.2f} | {vir_ok:>6}")
        latex_discrepancy = r"\checkmark" if disc_pct < 1.0 else r"\times"
        latex_virial = r"\checkmark" if z_vir <= 3.0 else r"\times"
        latex_rows.append(f"  {t:.1f} & ${etot:.4f} \\pm {err_tot:.4f}$ & ${e_bose:.4f}$ & ${e_dist:.4f}$ & ${disc_pct:.2f}\\%$ & ${diff:+.4f} \\pm {err_diff:.4f}$ & ${latex_discrepancy}$ & ${latex_virial}$ \\\\")

    print("=" * 122)
    print("\n[TABLA LATEX LISTA PARA OVERLEAF / MEMORIA TFG]:")
    print("\\begin{table}[htbp]")
    print("  \\centering")
    print(f"  \\caption{{Verificación de observables de energía y consistencia virial en oscilador armónico 2D ($N={n_particles}$).}}")
    print("  \\begin{tabular}{cccccccc}")
    print("    \\toprule")
    print("    $T$ & $\\langle E \\rangle / N$ (PIMC) & $E_{\\mathrm{Bose}}/N$ & $E_{\\mathrm{dist}}/N$ & Error vs Bose & $\\Delta_{\\mathrm{vir}}$ & $<1\\%$ & $|z_{\\mathrm{vir}}|\\leq3$ \\\\")
    print("    \\midrule")
    for r in latex_rows:
        print(r)
    print("    \\bottomrule")
    print("  \\end{tabular}")
    print("\\end{table}\n")


# ==============================================================================
# PUNTO DE ENTRADA PRINCIPAL
# ==============================================================================
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--particles", type=int, default=5,
                        help="Número de partículas (por defecto: 5)")
    parser.add_argument("--temp", type=float, default=0.5,
                        help="Temperatura de referencia para pruebas de convergencia (por defecto: 0.5)")
    parser.add_argument("--run", action="store_true",
                        help="Forzar reejecución de simulaciones C si no existen")
    args = parser.parse_args()

    print("=" * 72)
    print("  GENERADOR DE FIGURAS Y COMPROBACIÓN DE ENERGÍA (PIMC WORM 2D)")
    print("=" * 72)
    print(f"  Partículas N : {args.particles}")
    print(f"  Temperatura  : {args.temp}")
    print(f"  Destinos     : {RESULTS_DIR}")
    print(f"                 {OVERLEAF_DIR}")
    print("=" * 72)

    # 1. Figura 1: Demostración cuántica del Swap
    print("\n[1/5] Generando Figura 1: Demostración Cuántica (Efecto Swap)...")
    plot_01_demostracion_swap(args.particles, args.temp)

    # 2. Figura 2: Comprobación del Teorema del Virial
    print("\n[2/5] Generando Figura 2: Comprobación del Teorema del Virial...")
    plot_02_comprobacion_virial(args.particles, args.temp)

    # 3. Figura 3: Verificación de Convergencia MCMC y Pull 3-Sigma
    print("\n[3/5] Generando Figura 3: Verificación de Convergencia MCMC (3-Sigma)...")
    plot_03_verificacion_convergencia_mcmc(args.particles, args.temp)

    # 4. Figura 4: Verificación Termodinámica E(T) y Fracción de Condensado
    print("\n[4/5] Generando Figura 4: Termodinámica E(T) y Condensación...")
    plot_04_termodinamica_temperatura(args.particles)

    # 5. Figura 5: Dashboard Panorámico Completo
    print("\n[5/5] Generando Figura 5: Panel Completo Dashboard 2x2...")
    plot_05_panel_completo_dashboard(args.particles, args.temp)

    # Informe por terminal
    print_verification_report(args.particles)
    print(">>> ¡Proceso completado con éxito! Todas las figuras han sido actualizadas.")


if __name__ == "__main__":
    main()
