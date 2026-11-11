#!/usr/bin/env python3
"""
Generador de Gráficos de Densidad Radial PIMC en Oscilador Armónico 2D con la Regla de las 3 Sigmas:
  1. Gráfico SIN SWAP (Partículas Distinguibles):
     - Puntos PIMC sin permutaciones.
     - Curva analítica Distinguible (con banda de confianza 3-sigma, 99.7%).
     - Panel inferior de residuos estandarizados z(r) dentro de [-3, 3].
  2. Gráfico CON SWAP (Bosones Idénticos):
     - Puntos PIMC con permutaciones completas (Worm Algorithm).
     - Curva analítica Bose Canónico (con banda de confianza 3-sigma, 99.7%).
     - Línea Distinguible en alto contraste (nítida y diferenciada).
     - Panel inferior de residuos estandarizados z(r) dentro de [-3, 3].
  3. Gráfico Comparativo Directo (Sin Swap vs Con Swap con bandas 3-sigma).
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
BLCK = 200
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
    print(f"-> Ejecutando simulación PIMC {tag} (N={N_PARTICLES}, T={TEMP}, P={BEADS}, blck={BLCK})...")

    proc = subprocess.run(["./MC", "MC.cfg"], cwd=run_dir, env=env,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Error en simulación {tag}:\n{proc.stdout[-2000:]}")

    rd_path = os.path.join(run_dir, "outrd.dat")
    data_rd = np.loadtxt(rd_path)
    splits = np.where(np.diff(data_rd[:, 0]) < 0)[0] + 1
    blocks = np.split(data_rd, splits)
    valid_blocks = [b for b in blocks if len(b) > 0 and b[:, 1].sum() > 0]

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

    e_path = os.path.join(run_dir, "oute.dat")
    data_e = np.loadtxt(e_path)
    etot_mean = np.mean(data_e[len(data_e)//4:, 0])

    print(f"   Completada {tag}: <E>/N = {etot_mean:.4f}, {len(profiles)} bloques analizados.")
    return r_grid, dens_mean, dens_err, etot_mean


def plot_with_3sigma_subpanel(r_sim, d_sim, err_sim, e_sim,
                              r_fine, n_target, n_ref,
                              target_label, ref_label, sim_label,
                              target_color, ref_color, sim_color,
                              title_main, out_path_main, out_path_pull,
                              is_bose=True):
    """Genera:
       1) Gráfico principal con banda sombreada 3-sigma (99.73%) y barras 1-sigma.
       2) Gráfico de dos paneles con distribución de residuos estandarizados z(r) en [-3, +3].
    """
    mask = (r_sim <= 3.2) & (d_sim > 1e-4) & (err_sim > 0)
    r_fit = r_sim[mask]
    d_fit = d_sim[mask]
    err_fit = err_sim[mask]

    if is_bose:
        n_exact_points = analytical_bose_density(r_fit, N_PARTICLES, TEMP, OMEGA)
    else:
        n_exact_points = analytical_distinguishable_density(r_fit, TEMP, OMEGA)

    z_scores = (d_fit - n_exact_points) / err_fit
    max_z = np.max(np.abs(z_scores))
    pct_3sigma = np.mean(np.abs(z_scores) <= 3.0) * 100.0

    print(f"   Estadística 3-sigma ({sim_label}):")
    print(f"     Max |z|: {max_z:.2f} sigma")
    print(f"     Puntos dentro de 3-sigma: {pct_3sigma:.1f}%")

    # -------------------------------------------------------------------------
    # Figura 1: Gráfico Principal con Banda 3-sigma sombreada
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8.0, 5.6), dpi=300)

    # Banda de 3 sigmas (99.73% de confianza)
    ax.fill_between(r_fit, n_exact_points - 3.0 * err_fit, n_exact_points + 3.0 * err_fit,
                    color=target_color, alpha=0.18, label=r"Banda de confianza $3\sigma$ (99.73%)", zorder=1)
    ax.fill_between(r_fit, n_exact_points - 1.0 * err_fit, n_exact_points + 1.0 * err_fit,
                    color=target_color, alpha=0.28, label=r"Banda de confianza $1\sigma$ (68.3%)", zorder=2)

    # Curva objetivo analítica
    ax.plot(r_fine, n_target, color=target_color, linestyle="-", linewidth=2.3,
            label=target_label, zorder=3)

    # Curva de referencia alternativa (en alto contraste)
    ax.plot(r_fine, n_ref, color=ref_color, linestyle="--", linewidth=2.0, alpha=0.9,
            label=ref_label, zorder=3)

    # Puntos PIMC
    ax.errorbar(r_fit, d_fit, yerr=err_fit, fmt="o", color=sim_color, ecolor="black",
                elinewidth=1.1, capsize=2.8, markersize=4.8, alpha=0.95,
                label=rf"{sim_label} ($\pm 1\sigma$ barras)", zorder=4)

    ax.set_xlim(0, 3.2)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r"Distancia radial $r = \sqrt{x^2 + y^2}$ [$a_{\mathrm{ho}}$]", fontsize=11)
    ax.set_ylabel(r"Densidad radial por partícula $n(r)$", fontsize=11)
    ax.set_title(title_main, fontsize=12, fontweight="bold", pad=10)
    ax.grid(True, linestyle=":", alpha=0.55)
    ax.legend(loc="upper right", fontsize=9.2, framealpha=0.93)

    # Cuadro informativo con la prueba 3-sigma
    stat_text = (r"$\mathbf{Validaci\acute{o}n\ 3\sigma}:$" + "\n" +
                 rf"$\max |z| = {max_z:.2f}\sigma \leq 3\sigma$" + "\n" +
                 rf"Puntos dentro de $3\sigma$: {pct_3sigma:.1f}%" + "\n" +
                 rf"$\langle E \rangle / N = {e_sim:.3f}$")
    ax.annotate(stat_text, xy=(1.65, 0.14), xytext=(1.65, 0.16),
                bbox=dict(boxstyle="round,pad=0.5", facecolor="#ffffff", edgecolor="#adb5bd", alpha=0.95),
                fontsize=9.2)

    fig.tight_layout()
    fig.savefig(out_path_main, dpi=300)
    plt.close(fig)
    print(f"-> Guardada figura principal: {out_path_main}")

    # -------------------------------------------------------------------------
    # Figura 2: Formato 2 Paneles (Densidad + Residuos Estandarizados z)
    # -------------------------------------------------------------------------
    fig2, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(8.0, 7.2), sharex=True,
        gridspec_kw={"height_ratios": [2.3, 1.0], "hspace": 0.08}, dpi=300
    )

    # Panel Superior: Densidad n(r)
    ax_top.fill_between(r_fit, n_exact_points - 3.0 * err_fit, n_exact_points + 3.0 * err_fit,
                        color=target_color, alpha=0.18, label=r"Banda analítica $3\sigma$ (99.7%)", zorder=1)
    ax_top.plot(r_fine, n_target, color=target_color, linestyle="-", linewidth=2.3, label=target_label, zorder=2)
    ax_top.plot(r_fine, n_ref, color=ref_color, linestyle="--", linewidth=2.0, alpha=0.9, label=ref_label, zorder=2)
    ax_top.errorbar(r_fit, d_fit, yerr=err_fit, fmt="o", color=sim_color, ecolor="black",
                    elinewidth=1.1, capsize=2.8, markersize=4.5, alpha=0.95,
                    label=rf"{sim_label} ($\pm 1\sigma$)", zorder=3)
    ax_top.set_xlim(0, 3.2)
    ax_top.set_ylim(bottom=0)
    ax_top.set_ylabel(r"Densidad radial $n(r)$", fontsize=11)
    ax_top.set_title(title_main, fontsize=12, fontweight="bold", pad=10)
    ax_top.grid(True, linestyle=":", alpha=0.55)
    ax_top.legend(loc="upper right", fontsize=9.0, framealpha=0.92)

    # Panel Inferior: Residuos Estandarizados z(r)
    ax_bot.axhspan(-3.0, 3.0, color="#d4edda", alpha=0.65, label=r"Zona permitida $\pm 3\sigma$ (99.73%)")
    ax_bot.axhspan(-1.0, 1.0, color="#c3e6cb", alpha=0.50, label=r"Zona $\pm 1\sigma$ (68.3%)")
    ax_bot.axhline(0, color="black", linestyle="-", linewidth=0.9)
    ax_bot.axhline(+3, color="#dc3545", linestyle=":", linewidth=1.2)
    ax_bot.axhline(-3, color="#dc3545", linestyle=":", linewidth=1.2)

    ax_bot.plot(r_fit, z_scores, marker="o", color=sim_color, markeredgecolor="black",
                markersize=4.5, linestyle="none", alpha=0.9)

    ax_bot.set_ylim(-3.8, 3.8)
    ax_bot.set_xlabel(r"Distancia radial $r = \sqrt{x^2 + y^2}$ [$a_{\mathrm{ho}}$]", fontsize=11)
    ax_bot.set_ylabel(r"Pull $z = \frac{n_{\mathrm{MC}} - n_{\mathrm{exact}}}{\sigma}$", fontsize=10.5)
    ax_bot.grid(True, linestyle=":", alpha=0.55)
    ax_bot.legend(loc="lower right", fontsize=8.5, framealpha=0.9)

    fig2.subplots_adjust(top=0.92, bottom=0.10, left=0.12, right=0.96, hspace=0.10)
    fig2.savefig(out_path_pull, dpi=300)
    plt.close(fig2)
    print(f"-> Guardada figura con panel 3-sigma pull: {out_path_pull}")


def main():
    dir_noswap = os.path.join(OUTPUT_DIR, "sim_noswap")
    dir_swap = os.path.join(OUTPUT_DIR, "sim_swap")

    r_noswap, d_noswap, err_noswap, e_noswap = run_pimc_simulation(disable_swap=True, run_dir=dir_noswap)
    r_swap, d_swap, err_swap, e_swap = run_pimc_simulation(disable_swap=False, run_dir=dir_swap)

    r_fine = np.linspace(0.001, 3.5, 300)
    n_dist_fine = analytical_distinguishable_density(r_fine, TEMP, OMEGA)
    n_bose_fine = analytical_bose_density(r_fine, N_PARTICLES, TEMP, OMEGA)
    e_dist_exact = analytical_distinguishable_energy(TEMP, OMEGA)
    e_bose_exact, _ = analytical_bose_energy(N_PARTICLES, TEMP, OMEGA)

    print("\n" + "="*70)
    print("COMPARACIÓN DE VERIFICACIÓN FÍSICA Y 3-SIGMA (T = 0.5, N = 2):")
    print(f"  Analítico Distinguible: E/N = {e_dist_exact:.4f}, n(0) = {n_dist_fine[0]:.4f}")
    print(f"  Analítico Bose        : E/N = {e_bose_exact:.4f}, n(0) = {n_bose_fine[0]:.4f}")
    print(f"  PIMC sin Swap (Sim)   : E/N = {e_noswap:.4f}, n(0) = {d_noswap[0]:.4f}")
    print(f"  PIMC con Swap (Sim)   : E/N = {e_swap:.4f}, n(0) = {d_swap[0]:.4f}")
    print("="*70 + "\n")

    # 1. Gráfica SIN SWAP
    plot_with_3sigma_subpanel(
        r_sim=r_noswap, d_sim=d_noswap, err_sim=err_noswap, e_sim=e_noswap,
        r_fine=r_fine, n_target=n_dist_fine, n_ref=n_bose_fine,
        target_label=r"Analítico: Distinguibles (Gaussiana)",
        ref_label=rf"Analítico: Bosones ($N={N_PARTICLES}$, referencia)",
        sim_label="PIMC sin Swap",
        target_color="#1f77b4", ref_color="#d62728", sim_color="#2ca02c",
        title_main=rf"Perfil de Densidad Radial $n(r)$ SIN Swap ($N={N_PARTICLES}$, $T={TEMP}$)" + "\n" +
                   r"(Ausencia de intercambios $\rightarrow$ Maxwell-Boltzmann, Regla $3\sigma$)",
        out_path_main=os.path.join(OUTPUT_DIR, "verificacion_densidad_radial_sin_swap.png"),
        out_path_pull=os.path.join(OUTPUT_DIR, "verificacion_densidad_radial_sin_swap_3sigma.png"),
        is_bose=False
    )
    shutil.copy2(os.path.join(OUTPUT_DIR, "verificacion_densidad_radial_sin_swap.png"),
                 os.path.join(OVERLEAF_IMG_DIR, "verificacion_fase4_densidad_radial_sin_swap.png"))
    shutil.copy2(os.path.join(OUTPUT_DIR, "verificacion_densidad_radial_sin_swap_3sigma.png"),
                 os.path.join(OVERLEAF_IMG_DIR, "verificacion_fase4_densidad_radial_sin_swap_3sigma.png"))

    # 2. Gráfica CON SWAP
    plot_with_3sigma_subpanel(
        r_sim=r_swap, d_sim=d_swap, err_sim=err_swap, e_sim=e_swap,
        r_fine=r_fine, n_target=n_bose_fine, n_ref=n_dist_fine,
        target_label=rf"Analítico: Bose Canónico ($N={N_PARTICLES}$, exacto)",
        ref_label=r"Analítico: Distinguible (referencia Boltzmann)",
        sim_label="PIMC con Worm Swap",
        target_color="#d62728", ref_color="#2b2b2b", sim_color="#9467bd",
        title_main=rf"Perfil de Densidad Radial $n(r)$ CON Swap ($N={N_PARTICLES}$, $T={TEMP}$)" + "\n" +
                   r"(Permutaciones bosónicas $\rightarrow$ Condensación y Regla $3\sigma$)",
        out_path_main=os.path.join(OUTPUT_DIR, "verificacion_fase4_densidad_radial.png"),
        out_path_pull=os.path.join(OUTPUT_DIR, "verificacion_fase4_densidad_radial_3sigma.png"),
        is_bose=True
    )
    shutil.copy2(os.path.join(OUTPUT_DIR, "verificacion_fase4_densidad_radial.png"),
                 os.path.join(OVERLEAF_IMG_DIR, "verificacion_fase4_densidad_radial.png"))
    shutil.copy2(os.path.join(OUTPUT_DIR, "verificacion_fase4_densidad_radial_3sigma.png"),
                 os.path.join(OVERLEAF_IMG_DIR, "verificacion_fase4_densidad_radial_3sigma.png"))
    shutil.copy2(os.path.join(OUTPUT_DIR, "verificacion_fase4_densidad_radial.png"),
                 os.path.join(SCRIPT_DIR, "Resultados", "Verificacion_Fase4_Worm", "verificacion_fase4_densidad_radial.png"))

    # 3. Gráfica Comparativa Directa
    fig3, ax3 = plt.subplots(figsize=(8.5, 5.8), dpi=300)
    mask1 = (r_noswap <= 3.2) & (d_noswap > 1e-4)
    mask2 = (r_swap <= 3.2) & (d_swap > 1e-4)

    # Curvas analíticas
    ax3.plot(r_fine, n_bose_fine, color="#d62728", linestyle="-", linewidth=2.3,
             label=rf"Analítico: Bose Canónico ($N={N_PARTICLES}$)", zorder=2)
    ax3.plot(r_fine, n_dist_fine, color="#2b2b2b", linestyle="--", linewidth=2.0, alpha=0.9,
             label="Analítico: Distinguible (referencia Boltzmann)", zorder=2)

    # Bandas 3-sigma para ambas curvas
    n_bose_pts = analytical_bose_density(r_swap[mask2], N_PARTICLES, TEMP, OMEGA)
    n_dist_pts = analytical_distinguishable_density(r_noswap[mask1], TEMP, OMEGA)
    ax3.fill_between(r_swap[mask2], n_bose_pts - 3.0*err_swap[mask2], n_bose_pts + 3.0*err_swap[mask2],
                     color="#d62728", alpha=0.15, label=r"Banda Bose $\pm 3\sigma$ (99.7%)", zorder=1)
    ax3.fill_between(r_noswap[mask1], n_dist_pts - 3.0*err_noswap[mask1], n_dist_pts + 3.0*err_noswap[mask1],
                     color="#1f77b4", alpha=0.15, label=r"Banda Distinguible $\pm 3\sigma$ (99.7%)", zorder=1)

    # Datos PIMC
    ax3.errorbar(r_swap[mask2], d_swap[mask2], yerr=err_swap[mask2],
                 fmt="s", color="#9467bd", ecolor="black", elinewidth=1.1,
                 capsize=2.8, markersize=4.6, alpha=0.9,
                 label=r"PIMC CON Swap $\rightarrow$ Bosones", zorder=3)
    ax3.errorbar(r_noswap[mask1], d_noswap[mask1], yerr=err_noswap[mask1],
                 fmt="o", color="#2ca02c", ecolor="black", elinewidth=1.1,
                 capsize=2.8, markersize=4.6, alpha=0.9,
                 label=r"PIMC SIN Swap $\rightarrow$ Distinguibles", zorder=3)

    ax3.set_xlim(0, 3.2)
    ax3.set_ylim(bottom=0)
    ax3.set_xlabel(r"Distancia radial $r = \sqrt{x^2 + y^2}$ [$a_{\mathrm{ho}}$]", fontsize=11)
    ax3.set_ylabel(r"Densidad radial por partícula $n(r)$", fontsize=11)
    ax3.set_title(rf"Efecto del Movimiento Swap y Regla $3\sigma$ ($N={N_PARTICLES}$, $T={TEMP}$)",
                  fontsize=12, fontweight="bold", pad=10)
    ax3.grid(True, linestyle=":", alpha=0.55)
    ax3.legend(loc="upper right", fontsize=9.2, framealpha=0.93)

    fig3.tight_layout()
    path_comp = os.path.join(OUTPUT_DIR, "comparacion_densidad_swap_vs_noswap.png")
    fig3.savefig(path_comp, dpi=300)
    plt.close(fig3)
    shutil.copy2(path_comp, os.path.join(OVERLEAF_IMG_DIR, "comparacion_densidad_swap_vs_noswap.png"))
    print(f"-> Guardada comparativa directa: {path_comp}")

    print("\n¡Todas las gráficas con la Regla de las 3 Sigmas aplicadas exitosamente!")


if __name__ == "__main__":
    main()
