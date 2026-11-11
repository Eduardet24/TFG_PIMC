#!/usr/bin/env python3
"""
Generador de Gráficos de Densidad Radial PIMC en Oscilador Armónico 2D:
  1. Gráfico SIN SWAP: Partículas distinguibles (sin intercambio, L=1).
     Muestra cómo los puntos PIMC caen exactamente sobre la curva analítica distinguible.
  2. Gráfico CON SWAP: Bosones idénticos (Worm Algorithm con permutaciones).
     Muestra cómo los puntos PIMC caen sobre Bose canónico, destacando claramente
     la línea distinguible para una comparación visual inmediata y nítida.
"""

import os
import sys
import shutil
import subprocess
import numpy as np
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from Script_Comparacion_N_Particulas_2D import (
    analytical_distinguishable_density,
    analytical_bose_density,
    analytical_distinguishable_energy,
    analytical_bose_energy,
)

OUTPUT_DIR = os.path.join(SCRIPT_DIR, "Resultados", "Verificacion_Densidad")
OVERLEAF_IMG_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "overleaf", "imágenes")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(OVERLEAF_IMG_DIR, exist_ok=True)

N_PARTICLES = 2
TEMP = 0.5
OMEGA = 1.0
BEADS = 30
WORM_C = 0.1 / (0.40 * N_PARTICLES * BEADS)
BLCK_HEATING = 30
BLCK = 150
NITER = 500


def run_pimc_simulation(disable_swap, run_dir):
    """Ejecuta una corrida PIMC aislando directorio, con o sin swap."""
    os.makedirs(run_dir, exist_ok=True)
    mc_bin = os.path.join(SCRIPT_DIR, "MC")
    if not os.path.exists(mc_bin):
        raise RuntimeError("No se encontró el ejecutable MC. Compilar primero.")
    
    shutil.copy2(mc_bin, os.path.join(run_dir, "MC"))
    cfg_src = os.path.join(SCRIPT_DIR, "MC.cfg")
    with open(cfg_src, "r") as f:
        cfg = f.read()

    import re
    cfg = re.sub(r"N=\s*\d+", f"N= {N_PARTICLES}", cfg)
    cfg = re.sub(r"T=\s*[\d\.]+", f"T= {TEMP}", cfg)
    cfg = re.sub(r"Npop=\s*\d+", f"Npop= {BEADS}", cfg)
    cfg = re.sub(r"NwalkersMax=\s*\d+", f"NwalkersMax= {BEADS * 2}", cfg)
    cfg = re.sub(r"worm_C=\s*[\d\.]+", f"worm_C= {WORM_C}", cfg)
    cfg = re.sub(r"blck_heating=\s*\d+", f"blck_heating= {BLCK_HEATING}", cfg)
    cfg = re.sub(r"blck=\s*[\d\s]+", f"blck= {BLCK}", cfg)
    cfg = re.sub(r"Niter=\s*\d+", f"Niter= {NITER}", cfg)
    cfg = re.sub(r"Nmeasure=\s*\d+", f"Nmeasure= 1", cfg)

    cfg_dst = os.path.join(run_dir, "MC.cfg")
    with open(cfg_dst, "w") as f:
        f.write(cfg)

    env = os.environ.copy()
    if disable_swap:
        env["WORM_DISABLE_SWAP"] = "1"
    else:
        env.pop("WORM_DISABLE_SWAP", None)

    tag = "SIN SWAP (Distinguibles)" if disable_swap else "CON SWAP (Bose)"
    print(f"-> Ejecutando simulación PIMC {tag} (N={N_PARTICLES}, T={TEMP}, P={BEADS})...")

    proc = subprocess.run(["./MC", "MC.cfg"], cwd=run_dir, env=env,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Error en simulación {tag}:\n{proc.stdout[-2000:]}")

    # Extraer perfil de outrd.dat
    rd_path = os.path.join(run_dir, "outrd.dat")
    data_rd = np.loadtxt(rd_path)
    splits = np.where(np.diff(data_rd[:, 0]) < 0)[0] + 1
    blocks = np.split(data_rd, splits)
    valid_blocks = [b for b in blocks if len(b) > 0 and b[:, 1].sum() > 0]
    
    # Descartar calentamiento
    sel_blocks = valid_blocks[len(valid_blocks) // 4:]
    r_grid = sel_blocks[-1][:, 0]

    profiles = []
    for b in sel_blocks:
        p = b[:, 1]
        try:
            norm = np.trapezoid(2.0 * np.pi * r_grid * p, r_grid)
        except AttributeError:
            norm = np.trapz(2.0 * np.pi * r_grid * p, r_grid)
        if norm > 0:
            profiles.append(p / norm)

    profiles = np.array(profiles)
    dens_mean = profiles.mean(axis=0)
    dens_err = profiles.std(axis=0, ddof=1) / np.sqrt(len(profiles))

    # Extraer energía media
    e_path = os.path.join(run_dir, "oute.dat")
    data_e = np.loadtxt(e_path)
    etot_mean = np.mean(data_e[len(data_e)//4:, 0])

    print(f"   Completada {tag}: <E>/N = {etot_mean:.4f}, bloques={len(profiles)}")
    return r_grid, dens_mean, dens_err, etot_mean


def main():
    dir_noswap = os.path.join(OUTPUT_DIR, "sim_noswap")
    dir_swap = os.path.join(OUTPUT_DIR, "sim_swap")

    r_noswap, d_noswap, err_noswap, e_noswap = run_pimc_simulation(disable_swap=True, run_dir=dir_noswap)
    r_swap, d_swap, err_swap, e_swap = run_pimc_simulation(disable_swap=False, run_dir=dir_swap)

    r_fine = np.linspace(0.001, 3.5, 250)
    n_dist_fine = analytical_distinguishable_density(r_fine, TEMP, OMEGA)
    n_bose_fine = analytical_bose_density(r_fine, N_PARTICLES, TEMP, OMEGA)
    e_dist_exact = analytical_distinguishable_energy(TEMP, OMEGA)
    e_bose_exact, n0_frac = analytical_bose_energy(N_PARTICLES, TEMP, OMEGA)

    print("\n" + "="*70)
    print("COMPARACIÓN FÍSICA A T = 0.5, N = 2:")
    print(f"  Analítico Distinguible: E/N = {e_dist_exact:.4f}, n(0) = {n_dist_fine[0]:.4f}")
    print(f"  Analítico Bose        : E/N = {e_bose_exact:.4f}, n(0) = {n_bose_fine[0]:.4f}")
    print(f"  PIMC sin Swap (Sim)   : E/N = {e_noswap:.4f}, n(0) = {d_noswap[0]:.4f}")
    print(f"  PIMC con Swap (Sim)   : E/N = {e_swap:.4f}, n(0) = {d_swap[0]:.4f}")
    print("="*70 + "\n")

    # =========================================================================
    # GRÁFICA 1: SIN SWAP (Partículas Distinguibles)
    # =========================================================================
    fig1, ax1 = plt.subplots(figsize=(8.0, 5.5), dpi=300)

    # Curva distinguible analítica (donde DEBEN caer los puntos)
    ax1.plot(r_fine, n_dist_fine, color="#1f77b4", linestyle="-", linewidth=2.2,
             label=r"Analítico: Distinguibles (Gaussiana térmica)", zorder=2)
    
    # Curva Bose mostrada como referencia para evidenciar la ausencia de intercambio
    ax1.plot(r_fine, n_bose_fine, color="#d62728", linestyle="--", linewidth=1.8, alpha=0.85,
             label=rf"Analítico: Bosones idénticos ($N={N_PARTICLES}$, referencia)", zorder=1)

    # Datos PIMC sin swap
    mask1 = r_noswap <= 3.5
    ax1.errorbar(r_noswap[mask1], d_noswap[mask1], yerr=err_noswap[mask1],
                 fmt="o", color="#2ca02c", ecolor="black", elinewidth=1.1,
                 capsize=2.5, markersize=4.5, alpha=0.9,
                 label=rf"PIMC sin Swap ($N={N_PARTICLES}$, $P={BEADS}$)", zorder=3)

    ax1.set_xlim(0, 3.5)
    ax1.set_ylim(bottom=0)
    ax1.set_xlabel(r"Distancia radial $r = \sqrt{x^2 + y^2}$ [$a_{\mathrm{ho}}$]", fontsize=11)
    ax1.set_ylabel(r"Densidad radial por partícula $n(r)$", fontsize=11)
    ax1.set_title(rf"Perfil de Densidad Radial $n(r)$ SIN Swap ($N={N_PARTICLES}$, $T={TEMP}$)" + "\n" +
                  r"(Ausencia de intercambios $\rightarrow$ Estadística de Maxwell-Boltzmann)",
                  fontsize=12, fontweight="bold", pad=10)
    ax1.grid(True, linestyle=":", alpha=0.55)
    ax1.legend(loc="upper right", fontsize=9.5, framealpha=0.92)

    # Anotación explicativa en la gráfica
    ax1.annotate(r"Sin Swap: $P(L=1) = 100\%$" + "\n" +
                 rf"$\langle E \rangle / N = {e_noswap:.3f} \approx E_{{\mathrm{{dist}}}} = {e_dist_exact:.3f}$",
                 xy=(1.8, 0.15), xytext=(1.8, 0.17),
                 bbox=dict(boxstyle="round,pad=0.5", facecolor="#f8f9fa", edgecolor="#6c757d", alpha=0.9),
                 fontsize=9.5)

    fig1.tight_layout()
    path_noswap = os.path.join(OUTPUT_DIR, "verificacion_densidad_radial_sin_swap.png")
    fig1.savefig(path_noswap, dpi=300)
    plt.close(fig1)
    print(f"-> Gráfica guardada: {path_noswap}")
    shutil.copy2(path_noswap, os.path.join(OVERLEAF_IMG_DIR, "verificacion_fase4_densidad_radial_sin_swap.png"))

    # =========================================================================
    # GRÁFICA 2: CON SWAP (Línea de Distinguible MÁS CLARA y nítida)
    # =========================================================================
    fig2, ax2 = plt.subplots(figsize=(8.0, 5.5), dpi=300)

    # Curva Bose canónico analítica (donde DEBEN caer los puntos con swap)
    ax2.plot(r_fine, n_bose_fine, color="#d62728", linestyle="-", linewidth=2.3,
             label=rf"Analítico: Bose Canónico ($N={N_PARTICLES}$, exacto)", zorder=2)

    # Línea de Distinguible: Se hace más clara, nítida y visible con estilo discontinuo de alto contraste
    ax2.plot(r_fine, n_dist_fine, color="#2b2b2b", linestyle="--", linewidth=2.0, alpha=0.9,
             label=r"Analítico: Distinguible (referencia Boltzmann)", zorder=1)

    # Datos PIMC con swap
    mask2 = r_swap <= 3.5
    ax2.errorbar(r_swap[mask2], d_swap[mask2], yerr=err_swap[mask2],
                 fmt="o", color="#9467bd", ecolor="black", elinewidth=1.1,
                 capsize=2.5, markersize=4.5, alpha=0.9,
                 label=rf"PIMC con Worm Swap ($N={N_PARTICLES}$, $P={BEADS}$)", zorder=3)

    ax2.set_xlim(0, 3.5)
    ax2.set_ylim(bottom=0)
    ax2.set_xlabel(r"Distancia radial $r = \sqrt{x^2 + y^2}$ [$a_{\mathrm{ho}}$]", fontsize=11)
    ax2.set_ylabel(r"Densidad radial por partícula $n(r)$", fontsize=11)
    ax2.set_title(rf"Perfil de Densidad Radial $n(r)$ CON Swap ($N={N_PARTICLES}$, $T={TEMP}$)" + "\n" +
                  r"(Permutaciones bosónicas $\rightarrow$ Condensación y pico central)",
                  fontsize=12, fontweight="bold", pad=10)
    ax2.grid(True, linestyle=":", alpha=0.55)
    ax2.legend(loc="upper right", fontsize=9.5, framealpha=0.92)

    # Anotación explicativa
    ax2.annotate(r"Con Swap: ciclos $L \geq 2$" + "\n" +
                 rf"$\langle E \rangle / N = {e_swap:.3f} \approx E_{{\mathrm{{bose}}}} = {e_bose_exact:.3f}$" + "\n" +
                 rf"Pico central $n(0)_{{\mathrm{{Bose}}}} > n(0)_{{\mathrm{{dist}}}}$",
                 xy=(1.8, 0.15), xytext=(1.8, 0.17),
                 bbox=dict(boxstyle="round,pad=0.5", facecolor="#f8f9fa", edgecolor="#6c757d", alpha=0.9),
                 fontsize=9.5)

    fig2.tight_layout()
    path_swap = os.path.join(OUTPUT_DIR, "verificacion_fase4_densidad_radial.png")
    fig2.savefig(path_swap, dpi=300)
    plt.close(fig2)
    print(f"-> Gráfica guardada: {path_swap}")
    shutil.copy2(path_swap, os.path.join(OVERLEAF_IMG_DIR, "verificacion_fase4_densidad_radial.png"))
    shutil.copy2(path_swap, os.path.join(SCRIPT_DIR, "Resultados", "Verificacion_Fase4_Worm", "verificacion_fase4_densidad_radial.png"))

    # =========================================================================
    # GRÁFICA 3: COMPARATIVA DIRECTA EN UN SOLO PANEL (SIN SWAP vs CON SWAP)
    # =========================================================================
    fig3, ax3 = plt.subplots(figsize=(8.5, 5.8), dpi=300)

    # Curvas analíticas
    ax3.plot(r_fine, n_bose_fine, color="#d62728", linestyle="-", linewidth=2.2,
             label=rf"Analítico: Bose Canónico ($N={N_PARTICLES}$)", zorder=1)
    ax3.plot(r_fine, n_dist_fine, color="#343a40", linestyle="--", linewidth=2.0,
             label="Analítico: Distinguible (Boltzmann)", zorder=1)

    # Datos PIMC
    ax3.errorbar(r_swap[mask2], d_swap[mask2], yerr=err_swap[mask2],
                 fmt="s", color="#9467bd", ecolor="black", elinewidth=1.1,
                 capsize=2.5, markersize=4.5, alpha=0.85,
                 label=r"PIMC CON Swap $\rightarrow$ Bosones", zorder=3)

    ax3.errorbar(r_noswap[mask1], d_noswap[mask1], yerr=err_noswap[mask1],
                 fmt="o", color="#2ca02c", ecolor="black", elinewidth=1.1,
                 capsize=2.5, markersize=4.5, alpha=0.85,
                 label=r"PIMC SIN Swap $\rightarrow$ Distinguibles", zorder=2)

    ax3.set_xlim(0, 3.5)
    ax3.set_ylim(bottom=0)
    ax3.set_xlabel(r"Distancia radial $r = \sqrt{x^2 + y^2}$ [$a_{\mathrm{ho}}$]", fontsize=11)
    ax3.set_ylabel(r"Densidad radial por partícula $n(r)$", fontsize=11)
    ax3.set_title(rf"Efecto del Movimiento Swap en el Perfil de Densidad 2D ($N={N_PARTICLES}$, $T={TEMP}$)",
                  fontsize=12, fontweight="bold", pad=10)
    ax3.grid(True, linestyle=":", alpha=0.55)
    ax3.legend(loc="upper right", fontsize=10, framealpha=0.92)

    fig3.tight_layout()
    path_comp = os.path.join(OUTPUT_DIR, "comparacion_densidad_swap_vs_noswap.png")
    fig3.savefig(path_comp, dpi=300)
    plt.close(fig3)
    print(f"-> Gráfica guardada: {path_comp}")
    shutil.copy2(path_comp, os.path.join(OVERLEAF_IMG_DIR, "comparacion_densidad_swap_vs_noswap.png"))

    print("\n¡Todas las gráficas de densidad se han generado exitosamente!")


if __name__ == "__main__":
    main()
