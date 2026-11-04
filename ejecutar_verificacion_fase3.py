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
import numpy as np
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RUN_DIR = os.environ.get("WORM_FASE3_RUN_DIR", SCRIPT_DIR)
OUTPUT_DIR = os.environ.get(
    "WORM_FASE3_OUTPUT_DIR",
    os.path.join(SCRIPT_DIR, "Resultados", "Verificacion_Fase3_Worm"),
)
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Parámetros físicos y de simulación
N_PARTICLES = 5
TEMP = 0.5
BETA = 1.0 / TEMP
OMEGA = 1.0
BEADS = 30  # P = 30
TAU = BETA / BEADS
WORM_C = 0.1 / (0.60 * N_PARTICLES * BEADS)
HEATING_BLOCKS = 20
MEASURE_BLOCKS = 100
N_ITER = 500

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
print(f"=================================================================\n")

# 1. Configurar MC.cfg temporal para la prueba
cfg_path = os.path.join(RUN_DIR, "MC.cfg")
cfg_backup_path = os.path.join(SCRIPT_DIR, "MC.cfg.bak")
shutil.copy(cfg_path, cfg_backup_path)

with open(cfg_path, "r") as f:
    cfg_text = f.read()

cfg_text = re.sub(r"N=\s*\d+", f"N= {N_PARTICLES}", cfg_text)
cfg_text = re.sub(r"T=\s*[\d\.]+", f"T= {TEMP}", cfg_text)
cfg_text = re.sub(r"Npop=\s*\d+", f"Npop= {BEADS}", cfg_text)
cfg_text = re.sub(r"NwalkersMax=\s*\d+", f"NwalkersMax= {BEADS * 2}", cfg_text)
cfg_text = re.sub(r"worm_C=\s*[\d\.]+", f"worm_C= {WORM_C}", cfg_text)
cfg_text = re.sub(r"blck_heating=\s*\d+", f"blck_heating= {HEATING_BLOCKS}", cfg_text)
cfg_text = re.sub(r"blck=\s*[\d\s]+", f"blck= {MEASURE_BLOCKS}", cfg_text)
cfg_text = re.sub(r"Niter=\s*\d+", f"Niter= {N_ITER}", cfg_text)
cfg_text = re.sub(r"Nmeasure=\s*\d+", f"Nmeasure= 1", cfg_text)

with open(cfg_path, "w") as f:
    f.write(cfg_text)

try:
    print("[1/4] Ejecutando simulación PIMC con Worm Fase 3 en C...")
    cmd = [os.environ.get("WORM_FASE3_EXECUTABLE", os.path.join(RUN_DIR, "MC")), "MC.cfg"]
    proc = subprocess.run(cmd, cwd=RUN_DIR, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    log_path = os.path.join(OUTPUT_DIR, "simulacion_terminal.log")
    with open(log_path, "w") as f:
        f.write(proc.stdout)
    print(f"  -> Simulación finalizada con código de salida {proc.returncode}")
    print(f"  -> Log guardado en: {log_path}")
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd, output=proc.stdout)

    for fname in ["oute.dat", "outrd.dat", "in3Dprev.in", "outworm_gap.dat"]:
        src = os.path.join(RUN_DIR, fname)
        if os.path.exists(src):
            shutil.copy(src, os.path.join(OUTPUT_DIR, fname))
            print(f"  -> Archivo de datos guardado: {fname}")

finally:
    if os.path.exists(cfg_backup_path):
        shutil.move(cfg_backup_path, cfg_path)

# 2. Analizar datos de energía (oute.dat)
print("\n[2/4] Procesando resultados de energía en Sector Z...")
oute_file = os.path.join(OUTPUT_DIR, "oute.dat")
data_e = np.loadtxt(oute_file)

etot_blocks = data_e[:, 0]
epot_blocks = data_e[:, 1]
ekin_vir_blocks = data_e[:, 2]

n_use = len(etot_blocks)
start_idx = int(0.1 * n_use)

etot_mean = np.mean(etot_blocks[start_idx:])
etot_err = np.std(etot_blocks[start_idx:], ddof=1) / np.sqrt(n_use - start_idx)

epot_mean = np.mean(epot_blocks[start_idx:])
epot_err = np.std(epot_blocks[start_idx:], ddof=1) / np.sqrt(n_use - start_idx)

ekin_mean = np.mean(ekin_vir_blocks[start_idx:])
ekin_err = np.std(ekin_vir_blocks[start_idx:], ddof=1) / np.sqrt(n_use - start_idx)

e_analitica_continuo = OMEGA / np.tanh(0.5 * BETA * OMEGA)
omega_eff = (2.0 / TAU) * np.arcsinh(0.5 * OMEGA * TAU)
e_analitica_discreta = omega_eff / np.tanh(0.5 * BETA * omega_eff)

print(f"  Energía Total PIMC <E>/N   : {etot_mean:.6f} +/- {etot_err:.6f}")
print(f"  Energía Potencial <Epot>/N : {epot_mean:.6f} +/- {epot_err:.6f}")
print(f"  Energía Cinética <Ekin>/N  : {ekin_mean:.6f} +/- {ekin_err:.6f}")
print(f"  Teorema del Virial <Epot>-<Ekin> : {epot_mean - ekin_mean:.6f} (debe ser ~0)")
print(f"  Analítica (Discreta P={BEADS}): {e_analitica_discreta:.6f}")
print(f"  Diferencia vs Analítica    : {abs(etot_mean - e_analitica_discreta) / e_analitica_discreta * 100.0:.3f}%")

# 3. Analizar estadísticas de Fase 3 en el log
print("\n[3/4] Extrayendo estadísticas de movimientos del Worm (Advance, Recede, Close, Open, Ira)...")
worm_lines = [line for line in proc.stdout.splitlines() if "Worm F3:" in line]
move_stats = {}
if worm_lines:
    last_worm_line = worm_lines[-1]
    print(f"  Último estado: {last_worm_line}")
    # Worm F3: Z=... G=... (fG=...) | Op: ...% | Cl: ...% | Adv: ...% | Rec: ...% | Ira: ...% | MaxGap: ...
    m = re.search(r"Z=(\d+)\s+G=(\d+)\s+\(fG=([\d\.]+)\)\s+\|\s+Op:\s+([\d\.]+)%\s+\|\s+Cl:\s+([\d\.]+)%\s+\|\s+Adv:\s+([\d\.]+)%\s+\|\s+Rec:\s+([\d\.]+)%\s+\|\s+Ira:\s+([\d\.]+)%\s+\|\s+MaxGap:\s+(\d+)", last_worm_line)
    if m:
        move_stats = {
            'z_steps': int(m.group(1)),
            'g_steps': int(m.group(2)),
            'f_g': float(m.group(3)),
            'pct_op': float(m.group(4)),
            'pct_cl': float(m.group(5)),
            'pct_adv': float(m.group(6)),
            'pct_rec': float(m.group(7)),
            'pct_ira': float(m.group(8)),
            'max_gap': int(m.group(9)),
        }
        print(f"  Pasos en Sector Z          : {move_stats['z_steps']:,}")
        print(f"  Pasos en Sector G          : {move_stats['g_steps']:,} ({move_stats['f_g']*100:.1f}%)")
        print(f"  Aceptación Open            : {move_stats['pct_op']:.1f}%")
        print(f"  Aceptación Close           : {move_stats['pct_cl']:.1f}%")
        print(f"  Aceptación Advance         : {move_stats['pct_adv']:.1f}%")
        print(f"  Aceptación Recede          : {move_stats['pct_rec']:.1f}%")
        print(f"  Aceptación IraMove (espacio): {move_stats['pct_ira']:.1f}%")
        print(f"  Máximo gap observado en tau: {move_stats['max_gap']} beads (de {BEADS} posibles)")

# 4. Procesar densidad radial (outrd.dat)
print("\n[4/4] Procesando perfil de densidad radial...")
outrd_file = os.path.join(OUTPUT_DIR, "outrd.dat")
r_grid, dens_mean_norm, dens_err_norm, dens_analitica = None, None, None, None
if os.path.exists(outrd_file):
    rd_raw = np.loadtxt(outrd_file)
    saltos = np.where(np.diff(rd_raw[:, 0]) < 0)[0] + 1
    bloques = np.split(rd_raw, saltos)
    bloques_utiles = bloques[len(bloques) // 2:]
    r_grid = bloques_utiles[-1][:, 0]
    dens_profiles = np.stack([b[:, 1] for b in bloques_utiles])
    dens_mean = np.mean(dens_profiles, axis=0)
    dens_err = np.std(dens_profiles, axis=0, ddof=1) / np.sqrt(len(bloques_utiles))

    sigma2 = (0.5 / OMEGA) / np.tanh(0.5 * BETA * OMEGA)
    dens_analitica = (1.0 / (2.0 * np.pi * sigma2)) * np.exp(-r_grid**2 / (2.0 * sigma2))
    norm_sim = np.trapz(2.0 * np.pi * r_grid * dens_mean, r_grid)
    if norm_sim > 0:
        dens_mean_norm = dens_mean / norm_sim
        dens_err_norm = dens_err / norm_sim
    else:
        dens_mean_norm = dens_mean
        dens_err_norm = dens_err

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

ax1.axhline(e_analitica_discreta, color="tab:red", linestyle="--", linewidth=1.8,
            label=f"Exacta discreta P={BEADS}: {e_analitica_discreta:.4f}")
ax1.axhline(e_analitica_continuo, color="black", linestyle=":", linewidth=1.5,
            label=f"Exacta continuo: {e_analitica_continuo:.4f}")

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
    ax.plot(r_grid, dens_analitica, "k-", linewidth=2.0, label="Analítico Exacto: Gaussiana Térmica")
    ax.errorbar(r_grid, dens_mean_norm, yerr=dens_err_norm, fmt="o", color="tab:blue",
                markersize=4.5, capsize=2.5, elinewidth=1.0, label=f"PIMC Sector Z ($N={N_PARTICLES}$, $P={BEADS}$)")

    ax.set_xlabel(r"Radio $r$ [$a_{\mathrm{ho}}$]", fontsize=11)
    ax.set_ylabel(r"Densidad radial $n(r)$", fontsize=11)
    ax.set_title(f"Fase 3 Worm: Perfil de Densidad en Sector Z\nOscilador Armónico 2D ($T={TEMP}$)",
                 fontsize=12, fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(loc="upper right", fontsize=10, framealpha=0.9)
    fig.tight_layout()
    fig_dens_path = os.path.join(OUTPUT_DIR, "verificacion_fase3_densidad_radial.png")
    fig.savefig(fig_dens_path, dpi=200)
    plt.close(fig)
    print(f"[Guardada figura]: {fig_dens_path}")

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
    gap_data = np.loadtxt(gap_file)
    if gap_data.ndim == 2 and len(gap_data) > 0:
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
