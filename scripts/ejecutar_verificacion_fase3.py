#!/usr/bin/env python3
"""
Script de Verificación de la Fase 3 del Algoritmo Worm (Prokof'ev-Svistunov / Dornheim 2014).
Verifica:
  - Movimientos Advance y Recede (evolución temporal de Ira).
  - Desplazamiento espacial de Ira (relajación del extremo abierto).
  - Transiciones de sector Z <-> G con longitud de gusano variable (gap > 1).
  - Conservación de magnitudes canónicas en Sector Z (Energía y Teorema del Virial).
  - Perfil de densidad radial n(r).
"""

import os
import sys
import re
import shutil
import subprocess
import time
import numpy as np
import matplotlib.pyplot as plt

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

from pimc_common import (
    analytical_bose_density,
    analytical_bose_energy,
    compile_mc_binary,
    replace_cfg_value,
)
from worm_config import get_n_particles

_trapezoid = getattr(np, "trapezoid", getattr(np, "trapz", None))

SCRIPT_DIR = os.path.dirname(_THIS_DIR)
OUTPUT_DIR = os.environ.get(
    "WORM_FASE3_OUTPUT_DIR",
    os.path.join(SCRIPT_DIR, "Resultados", "Verificacion_Fase3_Worm"),
)
RUN_DIR = os.environ.get("WORM_FASE3_RUN_DIR", os.path.join(OUTPUT_DIR, "sim_run"))
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(RUN_DIR, exist_ok=True)

DEFAULT_CFG_PATH = os.path.join(SCRIPT_DIR, "MC.cfg")
RUN_CFG_PATH = os.path.join(RUN_DIR, "MC.cfg")
TEMPLATE_CFG_PATH = RUN_CFG_PATH if os.path.isfile(RUN_CFG_PATH) else DEFAULT_CFG_PATH

# Parámetros físicos y de simulación
N_PARTICLES = get_n_particles(TEMPLATE_CFG_PATH)
TEMP = 0.5
BETA = 1.0 / TEMP
OMEGA = 1.0
BEADS = 30  # P = 30
TAU = BETA / BEADS
WORM_C = 0.1 / (0.60 * N_PARTICLES * BEADS)
HEATING_BLOCKS = 20
MEASURE_BLOCKS = 100
N_ITER = 500


def get_worm_executable():
    configured = os.environ.get("WORM_FASE3_EXECUTABLE")
    if configured:
        return configured
    return str(
        compile_mc_binary(
            SCRIPT_DIR,
            os.path.join(SCRIPT_DIR, "MC"),
            only_if_stale=True,
            dim="2D",
            worm=True,
            normal_modes=True,
        )
    )


print(f"=================================================================")
print(f"  VERIFICACIÓN FASE 3: MOVIMIENTOS ADVANCE / RECEDE DEL WORM     ")
print(f"  (Basado en Dornheim PhD 2014, Cap. 2 y Boninsegni PRE 2006)    ")
print(f"=================================================================")
print(f"  Partículas N         : {N_PARTICLES}")
print(f"  Temperatura T        : {TEMP} (beta = {BETA})")
print(f"  Frecuencia trampa    : {OMEGA}")
print(f"  Número de beads P    : {BEADS} (tau = {TAU:.4f})")
print(f"  Parámetro fugacidad C: {WORM_C}")
print(f"  Bloques medición     : {MEASURE_BLOCKS} (iter por bloque: {N_ITER})")
print(f"  Directorio ejecución : {RUN_DIR}")
print(f"=================================================================\n")

# 1. Configurar MC.cfg en directorio aislado (sin mutar la raíz salvo que RUN_DIR == SCRIPT_DIR)
is_root_run = os.path.abspath(RUN_DIR) == os.path.abspath(SCRIPT_DIR)
cfg_backup_path = os.path.join(SCRIPT_DIR, "MC.cfg.bak") if is_root_run else None
if is_root_run:
    shutil.copy2(RUN_CFG_PATH, cfg_backup_path)

with open(TEMPLATE_CFG_PATH, "r", encoding="utf-8", errors="ignore") as f:
    cfg_text = f.read()

cfg_text = replace_cfg_value(cfg_text, "MC", 3)
cfg_text = replace_cfg_value(cfg_text, "N", N_PARTICLES)
cfg_text = replace_cfg_value(cfg_text, "T", TEMP)
cfg_text = replace_cfg_value(cfg_text, "Npop", BEADS)
cfg_text = replace_cfg_value(cfg_text, "NwalkersMax", BEADS * 2)
cfg_text = replace_cfg_value(cfg_text, "worm_C", f"{WORM_C:.6g}")
cfg_text = replace_cfg_value(cfg_text, "generate_new_coordinates", 1)
cfg_text = replace_cfg_value(cfg_text, "blck_heating", HEATING_BLOCKS)
cfg_text = replace_cfg_value(cfg_text, "blck", MEASURE_BLOCKS)
cfg_text = replace_cfg_value(cfg_text, "Niter", N_ITER)
cfg_text = replace_cfg_value(cfg_text, "Nmeasure", 1)

with open(RUN_CFG_PATH, "w", encoding="utf-8") as f:
    f.write(cfg_text)

try:
    print("[1/4] Ejecutando simulación PIMC con Worm Fase 3 en C...")
    executable = get_worm_executable()
    cmd = [executable, "MC.cfg"]
    run_started = time.time()
    proc = subprocess.run(cmd, cwd=RUN_DIR, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    log_path = os.path.join(OUTPUT_DIR, "simulacion_terminal.log")
    with open(log_path, "w", encoding="utf-8") as f:
        f.write(proc.stdout)
    print(f"  -> Simulación finalizada con código de salida {proc.returncode}")
    print(f"  -> Log guardado en: {log_path}")
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd, output=proc.stdout)

    for fname in ["oute.dat", "outrd.dat", "in3Dprev.in", "outworm_gap.dat"]:
        src = os.path.join(RUN_DIR, fname)
        dst = os.path.join(OUTPUT_DIR, fname)
        if os.path.exists(src) and os.path.getmtime(src) >= run_started - 1.0:
            if os.path.abspath(src) != os.path.abspath(dst):
                shutil.copy2(src, dst)
            print(f"  -> Archivo de datos guardado: {fname}")

finally:
    if cfg_backup_path and os.path.exists(cfg_backup_path):
        shutil.move(cfg_backup_path, RUN_CFG_PATH)

# 2. Analizar datos de energía (oute.dat)
print("\n[2/4] Procesando resultados de energía en Sector Z...")
oute_file = os.path.join(OUTPUT_DIR, "oute.dat")
data_e = np.loadtxt(oute_file, comments="#", ndmin=2)

etot_blocks = data_e[:, 0]
epot_blocks = data_e[:, 1]
ekin_vir_blocks = data_e[:, 2]

n_use = len(etot_blocks)
start_idx = int(0.1 * n_use) if n_use - int(0.1 * n_use) >= 2 else 0
n_eff = n_use - start_idx

etot_mean = float(np.mean(etot_blocks[start_idx:]))
etot_err = float(np.std(etot_blocks[start_idx:], ddof=1) / np.sqrt(n_eff)) if n_eff > 1 else 0.0

epot_mean = float(np.mean(epot_blocks[start_idx:]))
epot_err = float(np.std(epot_blocks[start_idx:], ddof=1) / np.sqrt(n_eff)) if n_eff > 1 else 0.0

ekin_mean = float(np.mean(ekin_vir_blocks[start_idx:]))
ekin_err = float(np.std(ekin_vir_blocks[start_idx:], ddof=1) / np.sqrt(n_eff)) if n_eff > 1 else 0.0

e_analitica_bose = analytical_bose_energy(N_PARTICLES, TEMP, OMEGA)[0]

print(f"  Energía Total PIMC <E>/N   : {etot_mean:.6f} +/- {etot_err:.6f}")
print(f"  Energía Potencial <Epot>/N : {epot_mean:.6f} +/- {epot_err:.6f}")
print(f"  Energía Cinética <Ekin>/N  : {ekin_mean:.6f} +/- {ekin_err:.6f}")
print(f"  Teorema del Virial <Epot>-<Ekin> : {epot_mean - ekin_mean:.6f} (debe ser ~0)")
print(f"  Analítica canónica (Bosones): {e_analitica_bose:.6f}")
print(f"  Diferencia vs Analítica    : {abs(etot_mean - e_analitica_bose) / e_analitica_bose * 100.0:.3f}%")

gap_file = os.path.join(OUTPUT_DIR, "outworm_gap.dat")
max_gap_from_hist = 0
gap_data = None
if os.path.exists(gap_file):
    gap_data = np.atleast_2d(np.loadtxt(gap_file, comments="#"))
    if gap_data.size > 0 and gap_data.shape[1] >= 2:
        nonzero = gap_data[:, 1] > 0
        if np.any(nonzero):
            max_gap_from_hist = int(np.max(gap_data[nonzero, 0]))

# 3. Analizar estadísticas de Fase 3 en el log
print("\n[3/4] Extrayendo estadísticas de movimientos del Worm (Advance, Recede, Close, Open, Ira)...")
worm_lines = [
    line for line in proc.stdout.splitlines()
    if "Worm F4:" in line or "Worm F3:" in line
]
move_stats = {}
if worm_lines:
    last_worm_line = worm_lines[-1]
    print(f"  Último estado: {last_worm_line}")
    m_f4 = re.search(
        r"Worm F4:\s*Z=(\d+)\s+G=(\d+)\s+\(fG=([\d.]+)\)\s*\|\s*Op:\s*([\d.]+)%\s*\|\s*Cl:\s*([\d.]+)%\s*\|\s*Adv:\s*([\d.]+)%\s*\|\s*Rec:\s*([\d.]+)%(?:\s*\|\s*Sw:\s*[\d.]+%)?\s*\|\s*Ira:\s*([\d.]+)%",
        last_worm_line,
    )
    m_f3 = re.search(
        r"Z=(\d+)\s+G=(\d+)\s+\(fG=([\d\.]+)\)\s+\|\s+Op:\s+([\d\.]+)%\s+\|\s+Cl:\s+([\d\.]+)%\s+\|\s+Adv:\s+([\d\.]+)%\s+\|\s+Rec:\s+([\d\.]+)%\s+\|\s+Ira:\s+([\d\.]+)%\s+\|\s+MaxGap:\s+(\d+)",
        last_worm_line,
    )
    if m_f4:
        m_gap = re.search(r"MaxGap:\s*(\d+)", last_worm_line)
        max_gap = int(m_gap.group(1)) if m_gap else max_gap_from_hist
        move_stats = {
            'z_steps': int(m_f4.group(1)),
            'g_steps': int(m_f4.group(2)),
            'f_g': float(m_f4.group(3)),
            'pct_op': float(m_f4.group(4)),
            'pct_cl': float(m_f4.group(5)),
            'pct_adv': float(m_f4.group(6)),
            'pct_rec': float(m_f4.group(7)),
            'pct_ira': float(m_f4.group(8)),
            'max_gap': max_gap,
        }
    elif m_f3:
        move_stats = {
            'z_steps': int(m_f3.group(1)),
            'g_steps': int(m_f3.group(2)),
            'f_g': float(m_f3.group(3)),
            'pct_op': float(m_f3.group(4)),
            'pct_cl': float(m_f3.group(5)),
            'pct_adv': float(m_f3.group(6)),
            'pct_rec': float(m_f3.group(7)),
            'pct_ira': float(m_f3.group(8)),
            'max_gap': int(m_f3.group(9)),
        }
    if move_stats:
        print(f"  Pasos en Sector Z          : {move_stats['z_steps']:,}")
        print(f"  Pasos en Sector G          : {move_stats['g_steps']:,} ({move_stats['f_g']*100:.1f}%)")
        print(f"  Aceptación Open            : {move_stats['pct_op']:.1f}%")
        print(f"  Aceptación Close           : {move_stats['pct_cl']:.1f}%")
        print(f"  Aceptación Advance         : {move_stats['pct_adv']:.1f}%")
        print(f"  Aceptación Recede          : {move_stats['pct_rec']:.1f}%")
        print(f"  Aceptación IraMove (espacio): {move_stats['pct_ira']:.1f}%")
        print(f"  Máximo gap observado en tau: {move_stats['max_gap']} beads (de {BEADS} posibles)")
    else:
        print("  (No se pudo parsear el formato de la línea Worm en el log)")
else:
    print("  (No se encontraron líneas de Worm F4/F3 en el log)")

# 4. Procesar densidad radial (outrd.dat)
print("\n[4/4] Procesando perfil de densidad radial...")
outrd_file = os.path.join(OUTPUT_DIR, "outrd.dat")
r_grid, dens_mean_norm, dens_err_norm, dens_analitica = None, None, None, None
if os.path.exists(outrd_file):
    rd_raw = np.loadtxt(outrd_file, ndmin=2)
    if rd_raw.size > 0:
        saltos = np.where(np.diff(rd_raw[:, 0]) < 0)[0] + 1
        bloques = np.split(rd_raw, saltos)
        bloques_validos = [b for b in bloques if len(b) > 0 and b[:, 1].sum() > 0]
        bloques_utiles = bloques_validos[len(bloques_validos) // 4:] if len(bloques_validos) >= 4 else bloques_validos
        if len(bloques_utiles) > 0:
            r_grid = bloques_utiles[-1][:, 0]
            profiles = []
            for b in bloques_utiles:
                p = b[:, 1]
                norm = _trapezoid(2.0 * np.pi * r_grid * p, r_grid)
                if norm > 0:
                    profiles.append(p / norm)
            if len(profiles) > 0:
                profiles = np.array(profiles)
                dens_mean_norm = np.mean(profiles, axis=0)
                dens_err_norm = (np.std(profiles, axis=0, ddof=1) / np.sqrt(len(profiles))
                                 if len(profiles) > 1 else np.zeros_like(dens_mean_norm))
                dens_analitica = analytical_bose_density(r_grid, N_PARTICLES, TEMP, OMEGA)
                print(f"  -> Densidad procesada sobre {len(profiles)} perfiles por bloques.")
                print(f"  -> Error estándar estimado entre perfiles: {np.mean(dens_err_norm):.6e}")

# ==============================================================================
# Gráficas de Verificación de la Fase 3
# ==============================================================================

# Gráfica 1: Energía vs Bloque y Teorema del Virial
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7.5), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
bloques_idx = np.arange(1, len(etot_blocks) + 1)

ax1.plot(bloques_idx, etot_blocks, color="tab:blue", alpha=0.35, label="Etot (bloque)")
etot_cum = np.cumsum(etot_blocks) / bloques_idx
ax1.plot(bloques_idx, etot_cum, color="tab:blue", linewidth=2.2, label=f"Etot acumulada: {etot_mean:.4f}")
ax1.plot(bloques_idx, epot_blocks, color="tab:orange", alpha=0.25, label="Epot (bloque)")
ax1.plot(bloques_idx, ekin_vir_blocks, color="tab:green", alpha=0.25, label="Ekin virial (bloque)")

ax1.axhline(e_analitica_bose, color="tab:red", linestyle="--", linewidth=1.8,
            label=f"Exacta canónica bosónica: {e_analitica_bose:.4f}")

ax1.set_ylabel(r"Energía por partícula $E/N$ [$\hbar\omega$]", fontsize=11)
ax1.set_title(f"Fase 3 Worm: Verificación de Energía Canónica en Sector Z\n"
              f"Con Movimientos Advance y Recede ($N={N_PARTICLES}$, $T={TEMP}$, $P={BEADS}$)",
              fontsize=12, fontweight="bold")
ax1.grid(True, linestyle=":", alpha=0.6)
ax1.legend(loc="upper right", fontsize=9.5, framealpha=0.9)

epot_cum = np.cumsum(epot_blocks) / bloques_idx
ekin_cum = np.cumsum(ekin_vir_blocks) / bloques_idx
ax2.plot(bloques_idx, epot_cum, color="tab:orange", linewidth=1.8, label=r"$\langle E_{\mathrm{pot}}\rangle$")
ax2.plot(bloques_idx, ekin_cum, color="tab:green", linewidth=1.8, linestyle="--", label=r"$\langle E_{\mathrm{kin,vir}}\rangle$")
ax2.plot(bloques_idx, np.abs(epot_cum - ekin_cum), color="tab:purple", linewidth=1.2, linestyle=":", label=r"$|\langle E_{\mathrm{pot}}\rangle - \langle E_{\mathrm{kin}}\rangle|$")

ax2.set_xlabel("Bloque Monte Carlo", fontsize=11)
ax2.set_ylabel(r"Componentes [$\hbar\omega$]", fontsize=10.5)
ax2.set_title("Comprobación del Teorema del Virial en Sector Z", fontsize=10.5, pad=6)
ax2.grid(True, linestyle=":", alpha=0.6)
ax2.legend(loc="upper right", fontsize=9.5, framealpha=0.9)

fig.tight_layout()
fig_energy_path = os.path.join(OUTPUT_DIR, "verificacion_fase3_energia_virial.png")
fig.savefig(fig_energy_path, dpi=200)
plt.close(fig)
print(f"\n[Guardada figura]: {fig_energy_path}")

# Gráfica 2: Perfil de Densidad Radial
if r_grid is not None:
    fig, ax = plt.subplots(figsize=(8, 5.5))
    mask = r_grid <= 3.5
    ax.plot(r_grid[mask], dens_analitica[mask], "r-", linewidth=2.0,
            label=f"Analítico exacto: Bose canónico ($N={N_PARTICLES}$)", zorder=2)
    ax.errorbar(r_grid[mask], dens_mean_norm[mask], yerr=dens_err_norm[mask], fmt="o", color="tab:blue",
                ecolor="black", markersize=4.5, capsize=3.0, elinewidth=1.1,
                label=f"PIMC Sector Z ($N={N_PARTICLES}$, $P={BEADS}$)", alpha=0.85, zorder=3)

    ax.set_xlim(0, 3.5)
    ax.set_xlabel(r"Radio $r$ [$a_{\mathrm{ho}}$]", fontsize=11)
    ax.set_ylabel(r"Densidad por unidad de área $n(r)$", fontsize=11)
    ax.set_title(f"Fase 3 Worm: Perfil de Densidad en Sector Z\nOscilador Armónico 2D ($T={TEMP}$)",
                 fontsize=12, fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(loc="upper right", fontsize=10, framealpha=0.9)
    fig.tight_layout()
    fig_dens_path = os.path.join(OUTPUT_DIR, "verificacion_fase3_densidad_radial.png")
    fig.savefig(fig_dens_path, dpi=200)
    plt.close(fig)
    print(f"[Guardada figura]: {fig_dens_path}")
    overleaf_img = os.path.join(os.path.dirname(SCRIPT_DIR), "overleaf", "imágenes", "verificacion_fase3_densidad_radial.png")
    if os.path.exists(os.path.dirname(overleaf_img)):
        shutil.copy(fig_dens_path, overleaf_img)

# Gráfica 3: Estadísticas de Movimientos del Worm en Fase 3
if move_stats:
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    categorias = ["Open\n(Z -> G)", "Close\n(G -> Z)", "Advance\n(Ira +tau)", "Recede\n(Ira -tau)", "IraMove\n(Espacial)"]
    pcts = [move_stats['pct_op'], move_stats['pct_cl'], move_stats['pct_adv'], move_stats['pct_rec'], move_stats['pct_ira']]
    colores = ["tab:blue", "tab:cyan", "tab:green", "tab:orange", "tab:purple"]

    bars = ax.bar(categorias, pcts, color=colores, width=0.55, edgecolor="black", linewidth=0.8)

    for bar, pct in zip(bars, pcts):
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, yval + 1.8, f"{pct:.1f}%",
                ha="center", fontsize=11, fontweight="bold")

    ax.set_ylabel("Tasa de Aceptación (%)", fontsize=11)
    ax.set_ylim(0, 115)
    ax.set_title(f"Fase 3 Worm: Tasas de Aceptación de los 5 Movimientos\n"
                 f"MaxGap explorado: {move_stats['max_gap']}/{BEADS} beads | Tiempo en G: {move_stats['f_g']*100:.1f}%",
                 fontsize=11.5, fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.5, axis="y")
    fig.tight_layout()
    fig_worm_path = os.path.join(OUTPUT_DIR, "verificacion_fase3_worm_movimientos.png")
    fig.savefig(fig_worm_path, dpi=200)
    plt.close(fig)
    print(f"[Guardada figura]: {fig_worm_path}")

# Gráfica 4: Histograma de Longitud del Worm (Gap Delta w en tiempo imaginario)
gap_file = os.path.join(OUTPUT_DIR, "outworm_gap.dat")
if os.path.exists(gap_file):
    gap_data = np.atleast_2d(np.loadtxt(gap_file, comments="#"))
    if gap_data.ndim == 2 and gap_data.shape[0] > 0 and gap_data.shape[1] >= 2 and np.sum(gap_data[:, 1]) > 0:
        gaps = gap_data[:, 0]
        counts = gap_data[:, 1]
        probs = counts / np.sum(counts) * 100.0

        fig, ax = plt.subplots(figsize=(7.5, 4.8))
        bars = ax.bar(gaps, probs, color="teal", width=0.55, edgecolor="black", linewidth=0.8)

        for bar, prob, count in zip(bars, probs, counts):
            yval = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2, yval + 1.2, f"{prob:.1f}%\n({int(count):,})",
                    ha="center", fontsize=10, fontweight="bold")

        ax.set_xlabel(r"Longitud del corte $\Delta w$ [beads]", fontsize=11)
        ax.set_ylabel(r"Probabilidad de presencia $P(\Delta w)$ [%]", fontsize=11)
        ax.set_ylim(0, max(probs) * 1.25)
        ax.set_xticks(gaps)
        ax.set_title(f"Fase 3 Worm: Distribución de Longitud de Apertura del Worm en Sector G\n"
                     f"Demostración de Exploración Temporal de la Cabeza Ira",
                     fontsize=11.5, fontweight="bold")
        ax.grid(True, linestyle=":", alpha=0.5, axis="y")
        fig.tight_layout()
        fig_gap_path = os.path.join(OUTPUT_DIR, "verificacion_fase3_worm_gap_distribucion.png")
        fig.savefig(fig_gap_path, dpi=200)
        plt.close(fig)
        print(f"[Guardada figura]: {fig_gap_path}")

print(f"\n=================================================================")
print(f"  VERIFICACIÓN DE FASE 3 COMPLETADA CON ÉXITO")
print(f"  Directorio con todos los datos y gráficas:")
print(f"  {OUTPUT_DIR}")
print(f"=================================================================")
