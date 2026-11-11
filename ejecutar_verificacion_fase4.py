#!/usr/bin/env python3
"""
Script de Verificación de la Fase 4 del Algoritmo Worm (Prokof'ev-Svistunov / Dornheim 2014).
Verifica:
  - Movimiento Swap (intercambio local de enlaces y generación de permutaciones bosónicas).
    - Movimientos Advance y Recede con N>1 y gap variable del worm.
  - Distribución de longitudes de ciclo de permutación P(L) en Sector Z.
  - Caída de la energía desde el valor distinguible (E_dist/N = 1.313) hacia el valor bosónico (E_Bose/N = 1.093 a T=0.5).
  - Conservación del Teorema del Virial <Epot> = <Ekin> para polímeros permutados.
  - Perfil de densidad radial n(r) con condensación de Bose-Einstein.
"""

import os
import sys
import re
import shutil
import subprocess
import time
import glob
import numpy as np
import matplotlib.pyplot as plt
from worm_config import get_n_particles

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RUN_DIR = os.environ.get("WORM_FASE4_RUN_DIR", SCRIPT_DIR)
OUTPUT_DIR = os.environ.get(
    "WORM_FASE4_OUTPUT_DIR",
    os.path.join(SCRIPT_DIR, "Resultados", "Verificacion_Fase4_Worm"),
)
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Parámetros físicos y de simulación
N_PARTICLES = get_n_particles(os.path.join(RUN_DIR, "MC.cfg"), "WORM_FASE4_PARTICLES")
if N_PARTICLES < 2 and "WORM_PARTICLES" not in os.environ and "WORM_FASE4_PARTICLES" not in os.environ:
    N_PARTICLES = 2
TEMP = 0.5
BETA = 1.0 / TEMP
OMEGA = 1.0
BEADS = 30  # P = 30
TAU = BETA / BEADS
WORM_C = 0.1 / (0.40 * N_PARTICLES * BEADS)
HEATING_BLOCKS = 30
MEASURE_BLOCKS = 120
N_ITER = 500


def get_worm_executable():
    configured_executable = os.environ.get("WORM_FASE4_EXECUTABLE")
    if configured_executable:
        return configured_executable

    sources = sorted(glob.glob(os.path.join(SCRIPT_DIR, "*.c")))
    if not sources:
        raise RuntimeError(f"No se encontraron fuentes C en {SCRIPT_DIR}")

    executable = os.path.join(SCRIPT_DIR, "MC_snapshots")
    newest_source = max(os.path.getmtime(source) for source in sources)
    if not os.path.isfile(executable) or os.path.getmtime(executable) < newest_source:
        compiler = os.environ.get("CC", "gcc")
        command = [compiler, "-O2", "-fopenmp", "-std=c99", *sources,
                   "-lm", "-o", executable]
        print(f"  -> Compilando ejecutable Worm con snapshots: {compiler}")
        result = subprocess.run(
            command, cwd=SCRIPT_DIR, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True
        )
        if result.returncode != 0:
            raise RuntimeError(f"No se pudo compilar el ejecutable Worm:\n{result.stdout[-6000:]}")
        print(f"  -> Ejecutable actualizado: {executable}")

    return executable

print(f"=================================================================")
print(f"  VERIFICACIÓN FASE 4: MOVIMIENTO SWAP Y PERMUTACIONES BOSÓNICAS  ")
print(f"  (Basado en Dornheim PhD 2014, Cap. 2 y Boninsegni PRE 2006)    ")
print(f"=================================================================")
print(f"  Partículas N         : {N_PARTICLES}")
print(f"  Temperatura T        : {TEMP} (beta = {BETA})")
print(f"  Frecuencia trampa    : {OMEGA}")
print(f"  Número de beads P    : {BEADS} (tau = {TAU:.4f})")
print(f"  Parámetro fugacidad C: {WORM_C:.6f}")
print(f"  Bloques medición     : {MEASURE_BLOCKS} (iter por bloque: {N_ITER})")
print(f"=================================================================\n")

# Modelos analíticos exactos (de Script_Comparacion_N_Particulas_2D.py)
def canonical_bose_partition(n_particles, temp, omega=OMEGA):
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


def analytical_bose_cycle_fractions(n_particles, temp, omega=OMEGA):
    partition, _, _ = canonical_bose_partition(n_particles, temp, omega)
    beta = 1.0 / temp
    return np.array([
        (1.0 / (1.0 - np.exp(-cycle * beta * omega))**2)
        * partition[n_particles - cycle] / (n_particles * partition[n_particles])
        for cycle in range(1, n_particles + 1)
    ])


def blocked_ratio_stats(values, counts):
    values = np.asarray(values, dtype=float)
    counts = np.asarray(counts, dtype=float)
    if len(values) < 2 or np.any(counts < 0) or counts.sum() <= 0:
        raise ValueError("Se necesitan al menos dos bloques para estimar la incertidumbre")
    mean = np.average(values, weights=counts)
    residuals = counts * (values - mean)
    error = np.std(residuals, ddof=1) / (np.sqrt(len(residuals)) * counts.mean())
    return mean, error

def analytical_bose_density(r_grid, n_particles, temp, omega=OMEGA):
    partition, _, _ = canonical_bose_partition(n_particles, temp, omega)
    beta = 1.0 / temp
    normalization = n_particles * partition[n_particles]

    density = np.zeros_like(r_grid)
    for cycle in range(1, n_particles + 1):
        q = np.exp(-cycle * beta * omega)
        cycle_partition = 1.0 / (1.0 - q)**2
        weight = cycle_partition * partition[n_particles - cycle] / normalization
        sigma_k2 = (0.5 / omega) / np.tanh(0.5 * cycle * beta * omega)
        gauss = (1.0 / (2.0 * np.pi * sigma_k2)) * np.exp(-r_grid**2 / (2.0 * sigma_k2))
        density += weight * gauss
    return density

def analytical_distinguishable_energy(temp, omega=OMEGA):
    beta = 1.0 / temp
    return omega / np.tanh(0.5 * beta * omega)

def analytical_distinguishable_density(r_grid, temp, omega=OMEGA):
    beta = 1.0 / temp
    sigma2 = (0.5 / omega) / np.tanh(0.5 * beta * omega)
    return (1.0 / (2.0 * np.pi * sigma2)) * np.exp(-r_grid**2 / (2.0 * sigma2))

# 1. Configurar MC.cfg temporal para la prueba
cfg_path = os.path.join(RUN_DIR, "MC.cfg")
cfg_backup_path = os.path.join(SCRIPT_DIR, "MC.cfg.bak4")
shutil.copy(cfg_path, cfg_backup_path)

with open(cfg_path, "r") as f:
    cfg_text = f.read()

cfg_text = re.sub(r"N=\s*\d+", f"N= {N_PARTICLES}", cfg_text)
cfg_text = re.sub(r"T=\s*[\d\.]+", f"T= {TEMP}", cfg_text)
cfg_text = re.sub(r"Npop=\s*\d+", f"Npop= {BEADS}", cfg_text)
cfg_text = re.sub(r"NwalkersMax=\s*\d+", f"NwalkersMax= {BEADS * 2}", cfg_text)
cfg_text = re.sub(r"worm_C=\s*[\d\.]+", f"worm_C= {WORM_C:.6f}", cfg_text)
cfg_text = re.sub(r"blck_heating=\s*\d+", f"blck_heating= {HEATING_BLOCKS}", cfg_text)
cfg_text = re.sub(r"blck=\s*[\d\s]+", f"blck= {MEASURE_BLOCKS}", cfg_text)
cfg_text = re.sub(r"Niter=\s*\d+", f"Niter= {N_ITER}", cfg_text)
cfg_text = re.sub(r"Nmeasure=\s*\d+", f"Nmeasure= 1", cfg_text)

with open(cfg_path, "w") as f:
    f.write(cfg_text)

try:
    print("[1/4] Ejecutando simulación PIMC con Worm Fase 4 (Swap + Permutaciones) en C...")
    executable = get_worm_executable()
    print(f"  -> Ejecutable: {executable}")
    cmd = [executable, "MC.cfg"]
    run_started = time.time()
    proc = subprocess.run(cmd, cwd=RUN_DIR, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    log_path = os.path.join(OUTPUT_DIR, "simulacion_terminal.log")
    with open(log_path, "w") as f:
        f.write(proc.stdout)
    print(f"  -> Simulación finalizada con código de salida {proc.returncode}")
    print(f"  -> Log guardado en: {log_path}")
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd, output=proc.stdout)

    required_outputs = ["oute.dat", "outworm_perm.dat", "outworm_energy_blocks.dat",
                        "outworm_cycle_blocks.dat", "outworm_snapshots.dat"]
    for fname in required_outputs:
        src = os.path.join(RUN_DIR, fname)
        if not os.path.isfile(src) or os.path.getmtime(src) < run_started:
            raise RuntimeError(f"El archivo {fname} falta o no se regeneró en esta corrida")

    for fname in ["oute.dat", "outrd.dat", "in3Dprev.in", "outworm_gap.dat", "outworm_perm.dat",
                  "outworm_energy_blocks.dat", "outworm_cycle_blocks.dat", "outworm_snapshots.dat"]:
        src = os.path.join(RUN_DIR, fname)
        if os.path.exists(src):
            shutil.copy(src, os.path.join(OUTPUT_DIR, fname))
            print(f"  -> Archivo de datos guardado: {fname}")

finally:
    if os.path.exists(cfg_backup_path):
        shutil.move(cfg_backup_path, cfg_path)

# 2. Analizar datos de energía (oute.dat)
print("\n[2/4] Procesando resultados de energía por bloques MC en Sector Z...")
energy_blocks_file = os.path.join(OUTPUT_DIR, "outworm_energy_blocks.dat")
data_e = np.loadtxt(energy_blocks_file, skiprows=1)
if data_e.ndim == 1:
    data_e = data_e.reshape(1, -1)

etot_blocks = data_e[:, 2]
epot_blocks = data_e[:, 3]
ekin_vir_blocks = data_e[:, 4]

n_use = len(data_e)
useful_start = int(0.15 * n_use)
useful = data_e[useful_start:]
energy_counts = useful[:, 1]
etot_mean, etot_err = blocked_ratio_stats(useful[:, 2], energy_counts)
epot_mean, epot_err = blocked_ratio_stats(useful[:, 3], energy_counts)
ekin_mean, ekin_err = blocked_ratio_stats(useful[:, 4], energy_counts)

e_dist = analytical_distinguishable_energy(TEMP)
e_bose, n0_frac = analytical_bose_energy(N_PARTICLES, TEMP)

print(f"  Energía Total PIMC <E>/N   : {etot_mean:.6f} +/- {etot_err:.6f}")
print(f"  Energía Potencial <Epot>/N : {epot_mean:.6f} +/- {epot_err:.6f}")
print(f"  Energía Cinética <Ekin>/N  : {ekin_mean:.6f} +/- {ekin_err:.6f}")
print(f"  Virial <Epot>-<Ekin>       : {epot_mean - ekin_mean:.6f} (debe ser ~0)")
print(f"  -------------------------------------------------------------")
print(f"  Referencia Distinguible    : {e_dist:.6f}")
print(f"  Referencia Bose (Exacta)   : {e_bose:.6f}  (N0/N = {n0_frac*100:.1f}%)")
print(f"  Discrepancia vs Bose       : {abs(etot_mean - e_bose) / e_bose * 100.0:.3f}%")
print(f"  Caída vs Distinguible      : {(e_dist - etot_mean) / (e_dist - e_bose) * 100.0:.1f}% del salto cuántico")

# 3. Analizar estadísticas de Fase 4 en el log y ciclos de permutación
print("\n[3/4] Analizando estadísticas del Worm y ciclos de permutación...")
worm_lines = [line for line in proc.stdout.splitlines() if "Worm F4:" in line]
if not worm_lines:
    raise RuntimeError("La simulación no emitió estadísticas Worm F4")
print(f"  Último estado: {worm_lines[-1]}")
swap_match = re.search(r"\bSw:\s*([\d.]+)%", worm_lines[-1])
if swap_match is None or float(swap_match.group(1)) <= 0.0:
    raise RuntimeError("El movimiento Swap no registró ninguna aceptación")
swap_rate = float(swap_match.group(1))
if N_PARTICLES > 1:
    advance_match = re.search(r"\bAdv:\s*([\d.]+)%", worm_lines[-1])
    recede_match = re.search(r"\bRec:\s*([\d.]+)%", worm_lines[-1])
    if (advance_match is None or recede_match is None
            or float(advance_match.group(1)) <= 0.0
            or float(recede_match.group(1)) <= 0.0):
        raise RuntimeError("Advance y Recede deben registrar aceptaciones para N>1")
    print(f"  Aceptación Advance/Recede: {advance_match.group(1)}% / {recede_match.group(1)}%")

perm_file = os.path.join(OUTPUT_DIR, "outworm_perm.dat")
perm_data = None
if not os.path.exists(perm_file):
    raise FileNotFoundError("Falta outworm_perm.dat; no se puede validar la estructura de ciclos")
perm_data = np.loadtxt(perm_file, skiprows=1)
if perm_data.ndim == 1 and len(perm_data) > 0:
    perm_data = perm_data.reshape(1, -1)
if (perm_data.size == 0 or perm_data.shape[1] < 3
        or not np.isclose(np.sum(perm_data[:, 2]), 1.0, atol=1e-3)):
    raise RuntimeError("La distribución de ciclos no está normalizada o está vacía")

cycle_blocks = np.loadtxt(os.path.join(OUTPUT_DIR, "outworm_cycle_blocks.dat"), skiprows=1)
if cycle_blocks.ndim == 1:
    cycle_blocks = cycle_blocks.reshape(1, -1)
cycle_blocks = cycle_blocks[len(cycle_blocks) // 2:]
cycle_counts = cycle_blocks[:, 1]
cycle_means, cycle_errors = [], []
for cycle in range(N_PARTICLES):
    mean = np.average(cycle_blocks[:, cycle + 2], weights=cycle_counts)
    residuals = cycle_counts * (cycle_blocks[:, cycle + 2] - mean)
    error = np.std(residuals, ddof=1) / (np.sqrt(len(residuals)) * cycle_counts.mean())
    cycle_means.append(mean)
    cycle_errors.append(error)
cycle_expected = analytical_bose_cycle_fractions(N_PARTICLES, TEMP)
cycle_z_scores = [abs(mean - expected) / error
                  for mean, expected, error in zip(cycle_means, cycle_expected, cycle_errors)
                  if expected >= 0.01 and error > 0.0]
max_cycle_z = max(cycle_z_scores, default=0.0)
print(f"  Máxima discrepancia de ciclos: {max_cycle_z:.2f} errores estándar")
if max_cycle_z > 3.0:
    raise RuntimeError("La distribución de ciclos se aparta de la referencia canónica (>3 sigma)")
print(f"\n  Distribución de Ciclos de Permutación P(L):")
print(f"  {'Longitud L':>10} | {'Conteos Part':>14} | {'Fracción P(L)':>14}")
print(f"  " + "-" * 44)
for row in perm_data:
    L_val = int(row[0])
    cnt = int(row[1])
    frac = row[2]
    print(f"  {L_val:>10} | {cnt:>14,d} | {frac:>13.2%}")

# 4. Generar cuatro gráficas independientes de verificación
print("\n[4/4] Generando cuatro gráficas independientes de verificación...")
plot_counts = data_e[:, 1]
valid_energy = plot_counts > 0
plot_blocks = data_e[valid_energy, 0]
plot_counts = plot_counts[valid_energy]
plot_etot = etot_blocks[valid_energy]
plot_epot = epot_blocks[valid_energy]
plot_ekin = ekin_vir_blocks[valid_energy]
cum_e = np.cumsum(plot_etot * plot_counts) / np.cumsum(plot_counts)

fig, ax = plt.subplots(figsize=(8, 5.5))
ax.plot(plot_blocks, plot_etot, color="black", alpha=0.35, label="E(bloque) PIMC")
ax.plot(plot_blocks, cum_e, color="blue", linewidth=2.0, label="Media acumulada ponderada")
ax.axhline(e_bose, color="red", linestyle="--", linewidth=1.8, label=f"Bose Analítico ({e_bose:.4f})")
ax.axhline(e_dist, color="gray", linestyle=":", linewidth=1.5, label=f"Distinguible ({e_dist:.4f})")
ax.set_xlabel("Bloque de Medición", fontsize=11)
ax.set_ylabel(r"$\langle E \rangle / N$", fontsize=11)
ax.set_title(f"Convergencia de Energía (T={TEMP}, N={N_PARTICLES}, P={BEADS})", fontsize=12, fontweight="bold")
ax.legend(fontsize=9, loc="best")
ax.grid(True, alpha=0.3)
fig.tight_layout()
plot_path = os.path.join(OUTPUT_DIR, "verificacion_fase4_convergencia_energia.png")
fig.savefig(plot_path, dpi=200)
plt.close(fig)
print(f"  -> Gráfica guardada: {plot_path}")

fig, ax = plt.subplots(figsize=(8, 5.5))
if perm_data is not None and len(perm_data) > 0:
    L_vals = perm_data[:, 0].astype(int)
    fracs = perm_data[:, 2]
    bars = ax.bar(L_vals, fracs, color="teal", edgecolor="black", alpha=0.8, width=0.55)
    for bar, fraction in zip(bars, fracs):
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2.0, yval + 0.02, f"{fraction * 100:.1f}%",
                ha="center", va="bottom", fontsize=10)
    ax.set_xticks(L_vals)
    ax.set_ylim(0, max(fracs) * 1.25)
ax.set_xlabel("Longitud del Ciclo de Permutación L", fontsize=11)
ax.set_ylabel("Fracción de Partículas P(L)", fontsize=11)
ax.set_title("Distribución de Ciclos de Permutación P(L) (Sector Z)", fontsize=12, fontweight="bold")
ax.grid(True, alpha=0.3, axis="y")
fig.tight_layout()
plot_path = os.path.join(OUTPUT_DIR, "verificacion_fase4_ciclos_permutacion.png")
fig.savefig(plot_path, dpi=200)
plt.close(fig)
print(f"  -> Gráfica guardada: {plot_path}")

fig, ax = plt.subplots(figsize=(8, 5.5))
outrd_file = os.path.join(OUTPUT_DIR, "outrd.dat")
density_profile_count = 0
if os.path.exists(outrd_file):
    try:
        data_rd = np.loadtxt(outrd_file)
        if data_rd.ndim == 2 and data_rd.shape[1] >= 2:
            rd_splits = np.where(np.diff(data_rd[:, 0]) < 0)[0] + 1
            rd_blocks = np.split(data_rd, rd_splits)
            valid_blocks = [b for b in rd_blocks if len(b) > 0 and b[:, 1].sum() > 0]
            selected_rd_blocks = valid_blocks[len(valid_blocks) // 4:] if len(valid_blocks) >= 4 else valid_blocks
            r_sim = selected_rd_blocks[-1][:, 0]
            profiles = []
            for b in selected_rd_blocks:
                p = b[:, 1]
                norm = np.trapz(2.0 * np.pi * r_sim * p, r_sim)
                if norm > 0:
                    profiles.append(p / norm)
            if len(profiles) > 0:
                profiles = np.array(profiles)
                density_profile_count = len(profiles)
                dens_sim = profiles.mean(axis=0)
                dens_err = (profiles.std(axis=0, ddof=1) / np.sqrt(len(profiles))
                            if len(profiles) > 1 else np.zeros_like(dens_sim))
                mask = r_sim <= 3.5
                ax.errorbar(r_sim[mask], dens_sim[mask], yerr=dens_err[mask], fmt='o', color='purple',
                            ecolor='black', elinewidth=1.1, markersize=3.5, capsize=3,
                            label='PIMC (error estándar entre bloques)', alpha=0.85, zorder=4)

                r_grid = np.linspace(0, 3.5, 200)
                n_bose = analytical_bose_density(r_grid, N_PARTICLES, TEMP)
                n_dist = analytical_distinguishable_density(r_grid, TEMP)
                ax.plot(r_grid, n_bose, color='red', linewidth=2.0, label='Bose Canónico (exacto)', zorder=2)
                ax.plot(r_grid, n_dist, color='#2b2b2b', linestyle='--', linewidth=2.0, alpha=0.9, label='Distinguible (referencia Boltzmann)', zorder=1)
                ax.set_xlim(0, 3.5)
                print(f"  -> Perfil radial promediado sobre {density_profile_count} bloques.")
                print(f"  -> Error estándar estimado entre bloques: {np.mean(dens_err):.6e}")
    except Exception as ex:
        print(f"  [Aviso] No se pudo graficar densidad radial: {ex}")
ax.set_xlabel("Radio r", fontsize=11)
ax.set_ylabel("Densidad por unidad de área n(r)", fontsize=11)
ax.set_title(f"Perfil de Densidad Radial n(r) (N={N_PARTICLES}, T={TEMP}, {density_profile_count} bloques)",
             fontsize=12, fontweight="bold")
ax.legend(fontsize=9, loc="best")
ax.grid(True, alpha=0.3)
fig.tight_layout()
plot_path = os.path.join(OUTPUT_DIR, "verificacion_fase4_densidad_radial.png")
fig.savefig(plot_path, dpi=200)
plt.close(fig)
print(f"  -> Gráfica guardada: {plot_path}")
overleaf_img = os.path.join(os.path.dirname(SCRIPT_DIR), "overleaf", "imágenes", "verificacion_fase4_densidad_radial.png")
if os.path.exists(os.path.dirname(overleaf_img)):
    shutil.copy(plot_path, overleaf_img)

fig, ax = plt.subplots(figsize=(8, 5.5))
ax.plot(plot_blocks, plot_epot, color="forestgreen", label=f"Epot = {epot_mean:.4f}")
ax.plot(plot_blocks, plot_ekin, color="coral", linestyle="--", label=f"Ekin (Virial) = {ekin_mean:.4f}")
ax.plot(plot_blocks, plot_epot - plot_ekin, color="crimson", linewidth=1.2,
        label=f"Epot - Ekin = {epot_mean - ekin_mean:.4f}")
ax.axhline(0, color="black", linestyle=":", alpha=0.7)
ax.set_xlabel("Bloque de Medición", fontsize=11)
ax.set_ylabel("Energía por Partícula", fontsize=11)
ax.set_title("Teorema del Virial <Epot> vs <Ekin>", fontsize=12, fontweight="bold")
ax.legend(fontsize=9, loc="best")
ax.grid(True, alpha=0.3)
fig.tight_layout()
plot_path = os.path.join(OUTPUT_DIR, "verificacion_fase4_teorema_virial.png")
fig.savefig(plot_path, dpi=200)
plt.close(fig)
print(f"  -> Gráfica guardada: {plot_path}")

legacy_plot_path = os.path.join(OUTPUT_DIR, "verificacion_fase4_worm.png")
if os.path.exists(legacy_plot_path):
    os.remove(legacy_plot_path)

print(f"\n=================================================================")
print(f"  RESULTADO DE LA VERIFICACIÓN DE FASE 4:                       ")
diff_pct = abs(etot_mean - e_bose) / e_bose * 100.0
virial_diff = abs(epot_mean - ekin_mean)
if diff_pct < 2.0 and virial_diff < 0.05 and swap_rate > 0.0 and max_cycle_z <= 3.0:
    print(f"  >>> ¡ÉXITO COMPLETO! <<<")
    print(f"  La simulación con Movimiento Swap reproduce la energía bosónica")
    print(f"  analítica dentro del {diff_pct:.2f}% de error.")
else:
    print(f"  Discrepancia vs Bose: {diff_pct:.2f}%, Virial diff: {virial_diff:.4f}")
print(f"=================================================================\n")
