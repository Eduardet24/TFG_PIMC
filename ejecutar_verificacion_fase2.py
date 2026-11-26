#!/usr/bin/env python3
"""
Script de Verificación de la Fase 2 del Algoritmo Worm (Prokof'ev-Svistunov / Dornheim 2014).
Ejecuta la simulación PIMC con transiciones Open/Close en oscilador armónico 2D,
extrae los datos reales de energía, densidad radial y estadísticas del worm,
y los compara contra las soluciones analíticas exactas de la bibliografía.
"""

import os
import sys
import re
import shutil
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from Script_Comparacion_N_Particulas_2D import analytical_bose_density, analytical_bose_energy
from worm_config import get_n_particles

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "Resultados", "Verificacion_Fase2_Worm")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Parámetros físicos y de simulación
N_PARTICLES = get_n_particles(os.path.join(SCRIPT_DIR, "MC.cfg"))
TEMP = 0.5
BETA = 1.0 / TEMP
OMEGA = 1.0
BEADS = 30  # P = 30
TAU = BETA / BEADS
WORM_C = 1.0
HEATING_BLOCKS = 20
MEASURE_BLOCKS = 100
N_ITER = 500

print(f"=================================================================")
print(f"  VERIFICACIÓN FASE 2: ALGORITMO WORM EN OSCILADOR ARMÓNICO 2D   ")
print(f"  (Basado en Dornheim PhD 2014, Cap. 2 y Ceperley RMP 1995)      ")
print(f"=================================================================")
print(f"  Partículas N         : {N_PARTICLES}")
print(f"  Temperatura T        : {TEMP} (beta = {BETA})")
print(f"  Frecuencia trampa    : {OMEGA}")
print(f"  Número de beads P    : {BEADS} (tau = {TAU:.4f})")
print(f"  Parámetro fugacidad C: {WORM_C}")
print(f"  Bloques medición     : {MEASURE_BLOCKS} (iter por bloque: {N_ITER})")
print(f"=================================================================\n")

# 1. Configurar MC.cfg temporal para la prueba
cfg_path = os.path.join(SCRIPT_DIR, "MC.cfg")
cfg_backup_path = os.path.join(SCRIPT_DIR, "MC.cfg.bak")
shutil.copy(cfg_path, cfg_backup_path)

with open(cfg_path, "r") as f:
    cfg_text = f.read()

# Modificar parámetros clave
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
    print("[1/4] Ejecutando simulación PIMC en C...")
    cmd = ["./MC", "MC.cfg"]
    proc = subprocess.run(cmd, cwd=SCRIPT_DIR, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    log_path = os.path.join(OUTPUT_DIR, "simulacion_terminal.log")
    with open(log_path, "w") as f:
        f.write(proc.stdout)
    print(f"  -> Simulación finalizada con código de salida {proc.returncode}")
    print(f"  -> Log guardado en: {log_path}")

    # Copiar archivos generados al directorio de resultados
    for fname in ["oute.dat", "outrd.dat", "in3Dprev.in"]:
        src = os.path.join(SCRIPT_DIR, fname)
        if os.path.exists(src):
            shutil.copy(src, os.path.join(OUTPUT_DIR, fname))
            print(f"  -> Archivo de datos guardado: {fname}")

finally:
    # Restaurar MC.cfg original
    if os.path.exists(cfg_backup_path):
        shutil.move(cfg_backup_path, cfg_path)

# 2. Analizar datos de energía (oute.dat)
print("\n[2/4] Procesando resultados de energía...")
oute_file = os.path.join(OUTPUT_DIR, "oute.dat")
data_e = np.loadtxt(oute_file)

# oute.dat: Etot, Epot, Ekin (virial), Ekin_prim
etot_blocks = data_e[:, 0]
epot_blocks = data_e[:, 1]
ekin_vir_blocks = data_e[:, 2]

# Descartar primeros 10% de bloques como equilibrado residual si fuera necesario
n_use = len(etot_blocks)
start_idx = int(0.1 * n_use)

etot_mean = np.mean(etot_blocks[start_idx:])
etot_err = np.std(etot_blocks[start_idx:], ddof=1) / np.sqrt(n_use - start_idx)

epot_mean = np.mean(epot_blocks[start_idx:])
epot_err = np.std(epot_blocks[start_idx:], ddof=1) / np.sqrt(n_use - start_idx)

ekin_mean = np.mean(ekin_vir_blocks[start_idx:])
ekin_err = np.std(ekin_vir_blocks[start_idx:], ddof=1) / np.sqrt(n_use - start_idx)

e_analitica_bose = analytical_bose_energy(N_PARTICLES, TEMP, OMEGA)[0]

print(f"  Energía Total PIMC <E>/N   : {etot_mean:.6f} +/- {etot_err:.6f}")
print(f"  Energía Potencial <Epot>/N : {epot_mean:.6f} +/- {epot_err:.6f}")
print(f"  Energía Cinética <Ekin>/N  : {ekin_mean:.6f} +/- {ekin_err:.6f}")
print(f"  Teorema del Virial <Epot>-<Ekin> : {epot_mean - ekin_mean:.6f} (debe ser ~0)")
print(f"  Analítica canónica (Bosones): {e_analitica_bose:.6f}")
diferencia_pct = abs(etot_mean - e_analitica_bose) / e_analitica_bose * 100.0
print(f"  Diferencia vs Analítica    : {diferencia_pct:.3f}%")

# 3. Analizar estadísticas del Worm en el log
print("\n[3/4] Extrayendo estadísticas de transición del Worm (Z <-> G)...")
worm_lines = [line for line in proc.stdout.splitlines() if "Worm F2:" in line]
if worm_lines:
    last_worm_line = worm_lines[-1]
    print(f"  Último estado: {last_worm_line}")
    # Extraer Z, G, frac_G, acc_open, acc_close
    match = re.search(r"Z=(\d+)\s+G=(\d+)\s+\(frac_G=([\d\.]+)\).*Open:\s+(\d+)/(\d+)\s+\(([\d\.]+)%\).*Close:\s+(\d+)/(\d+)\s+\(([\d\.]+)%\)", last_worm_line)
    if match:
        z_steps = int(match.group(1))
        g_steps = int(match.group(2))
        frac_g = float(match.group(3))
        acc_open_num = int(match.group(4))
        try_open_num = int(match.group(5))
        pct_open = float(match.group(6))
        acc_close_num = int(match.group(7))
        try_close_num = int(match.group(8))
        pct_close = float(match.group(9))

        print(f"  Pasos en Sector Z (cerrado): {z_steps:,}")
        print(f"  Pasos en Sector G (abierto): {g_steps:,}")
        print(f"  Fracción de tiempo en G    : {frac_g*100:.1f}%")
        print(f"  Aceptación Open (Z -> G)   : {acc_open_num}/{try_open_num} ({pct_open:.1f}%)")
        print(f"  Aceptación Close (G -> Z)  : {acc_close_num}/{try_close_num} ({pct_close:.1f}%)")
else:
    print("  (No se encontraron líneas de Worm F2 en el log)")

# 4. Analizar densidad radial (outrd.dat)
print("\n[4/4] Procesando perfil de densidad radial...")
outrd_file = os.path.join(OUTPUT_DIR, "outrd.dat")
r_grid, dens_mean_norm, dens_err_norm, dens_analitica = None, None, None, None
if os.path.exists(outrd_file):
    rd_raw = np.loadtxt(outrd_file)
    saltos = np.where(np.diff(rd_raw[:, 0]) < 0)[0] + 1
    bloques = np.split(rd_raw, saltos)
    bloques_validos = [b for b in bloques if len(b) > 0 and b[:, 1].sum() > 0]
    bloques_utiles = bloques_validos[len(bloques_validos) // 4:] if len(bloques_validos) >= 4 else bloques_validos
    if len(bloques_utiles) > 0:
        r_grid = bloques_utiles[-1][:, 0]
        profiles = []
        for b in bloques_utiles:
            p = b[:, 1]
            norm = np.trapz(2.0 * np.pi * r_grid * p, r_grid)
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
# Generar Figuras de Verificación
# ==============================================================================
# Gráfica 1: Energía vs Bloque y Teorema del Virial
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7.5), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
bloques_idx = np.arange(1, len(etot_blocks) + 1)

ax1.plot(bloques_idx, etot_blocks, color="tab:blue", alpha=0.45, label="Etot (bloque)")
# Media acumulada
etot_cum = np.cumsum(etot_blocks) / bloques_idx
ax1.plot(bloques_idx, etot_cum, color="tab:blue", linewidth=2.2, label=f"Etot acumulada: {etot_mean:.4f}")
ax1.plot(bloques_idx, epot_blocks, color="tab:orange", alpha=0.3, label="Epot (bloque)")
ax1.plot(bloques_idx, ekin_vir_blocks, color="tab:green", alpha=0.3, label="Ekin virial (bloque)")

ax1.axhline(e_analitica_bose, color="tab:red", linestyle="--", linewidth=1.8,
            label=f"Exacta canónica bosónica: {e_analitica_bose:.4f}")

ax1.set_ylabel(r"Energía por partícula $E/N$ [$\hbar\omega$]", fontsize=11)
ax1.set_title(f"Fase 2 Worm: Verificación de Energía Canónica (Sector Z)\n"
              f"Oscilador Armónico 2D ($N={N_PARTICLES}$, $T={TEMP}$, $P={BEADS}$)",
              fontsize=12, fontweight="bold")
ax1.grid(True, linestyle=":", alpha=0.6)
ax1.legend(loc="upper right", fontsize=9.5, framealpha=0.9)

# Panel 2: Comprobación del Teorema del Virial: Epot vs Ekin
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
fig_energy_path = os.path.join(OUTPUT_DIR, "verificacion_energia_virial.png")
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
    ax.set_title(f"Fase 2 Worm: Perfil de Densidad en Sector Z\nOscilador Armónico 2D ($T={TEMP}$)",
                 fontsize=12, fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(loc="upper right", fontsize=10, framealpha=0.9)
    fig.tight_layout()
    fig_dens_path = os.path.join(OUTPUT_DIR, "verificacion_densidad_radial.png")
    fig.savefig(fig_dens_path, dpi=200)
    plt.close(fig)
    print(f"[Guardada figura]: {fig_dens_path}")
    overleaf_img = os.path.join(os.path.dirname(SCRIPT_DIR), "overleaf", "imágenes", "verificacion_densidad_radial.png")
    if os.path.exists(os.path.dirname(overleaf_img)):
        shutil.copy(fig_dens_path, overleaf_img)

# Gráfica 3: Estadísticas del Worm
if worm_lines:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    categorias = ["Open (Z -> G)", "Close (G -> Z)"]
    aceptadas = [acc_open_num, acc_close_num]
    totales = [try_open_num, try_close_num]
    pcts = [pct_open, pct_close]

    x = np.arange(len(categorias))
    width = 0.35

    ax.bar(x - width/2, totales, width, label="Intentos", color="tab:gray", alpha=0.6)
    ax.bar(x + width/2, aceptadas, width, label="Aceptados", color="tab:blue")

    for i in range(len(categorias)):
        ax.text(x[i] + width/2, aceptadas[i] + max(totales)*0.02, f"{pcts[i]:.1f}%",
                ha="center", fontsize=11, fontweight="bold", color="tab:blue")

    ax.set_ylabel("Número de movimientos", fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels(categorias, fontsize=11)
    ax.set_title(f"Dinámica Worm: Sector Z ({z_steps:,} pasos) vs Sector G ({g_steps:,} pasos, {frac_g*100:.1f}%)",
                 fontsize=11.5, fontweight="bold")
    ax.legend(loc="upper right", fontsize=10)
    ax.grid(True, linestyle=":", alpha=0.5, axis="y")
    fig.tight_layout()
    fig_worm_path = os.path.join(OUTPUT_DIR, "verificacion_worm_estadisticas.png")
    fig.savefig(fig_worm_path, dpi=200)
    plt.close(fig)
    print(f"[Guardada figura]: {fig_worm_path}")

print(f"\n=================================================================")
print(f"  VERIFICACIÓN COMPLETADA CON ÉXITO")
print(f"  Directorio con todos los datos y gráficas:")
print(f"  {OUTPUT_DIR}")
print(f"=================================================================")

