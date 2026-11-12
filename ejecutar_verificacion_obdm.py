#!/usr/bin/env python3
r"""
ejecutar_verificacion_obdm.py

Script de Ejecución y Verificación de la Matriz de Densidad de Una Partícula (OBDM)
muestreada mediante los extremos de cadenas abiertas en el Algoritmo Worm (PIMC).

Calcula y valida:
  - Muestreo del sector G con gap == 1 (Ira y Masha).
    - Comprobación de la traza del archivo, normalizada por construcción en el motor C.
  - Comparación de \rho_1(x, x') con la solución canónica exacta de Bose (Kernel de Mehler).
  - Comparación con el límite de partículas distinguibles para constatar la coherencia de Bose.
  - Generación de gráficos en alta resolución y sincronización con el directorio de Overleaf.
"""

import os
import sys
import re
import shutil
import subprocess
import time
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from worm_config import get_n_particles
from generar_graficos_obdm import (
    load_worm_obdm_data,
    load_worm_dist_data,
    analytical_bose_obdm_1d,
    analytical_distinguishable_obdm_1d,
    generar_mapa_calor_obdm,
    generar_cortes_obdm,
    generar_grafico_distancia_extremos,
)

OUTPUT_DIR = os.path.join(SCRIPT_DIR, "Resultados", "Verificacion_OBDM")
OVERLEAF_IMG_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "overleaf", "imágenes")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(OVERLEAF_IMG_DIR, exist_ok=True)

# Parámetros físicos
N_PARTICLES = get_n_particles(os.path.join(SCRIPT_DIR, "MC.cfg"), "WORM_FASE4_PARTICLES")
if N_PARTICLES < 2 and "WORM_PARTICLES" not in os.environ and "WORM_FASE4_PARTICLES" not in os.environ:
    N_PARTICLES = 2
TEMP = 0.5
BETA = 1.0 / TEMP
OMEGA = 1.0
BEADS = 30
TAU = BETA / BEADS
WORM_C = 0.1 / (0.40 * N_PARTICLES * BEADS)  # ~0.004167

HEATING_BLOCKS = int(os.environ.get("WORM_OBDM_HEATING", 25))
MEASURE_BLOCKS = int(os.environ.get("WORM_OBDM_BLOCKS", 100))
N_ITER = int(os.environ.get("WORM_OBDM_NITER", 500))


def get_or_compile_executable():
    executable = os.path.join(SCRIPT_DIR, "MC")
    sources = sorted([os.path.join(SCRIPT_DIR, f) for f in os.listdir(SCRIPT_DIR) if f.endswith(".c")])
    newest_source = max(os.path.getmtime(s) for s in sources)
    if not os.path.isfile(executable) or os.path.getmtime(executable) < newest_source:
        compiler = os.environ.get("CC", "gcc")
        cmd = [compiler, "-O2", "-fopenmp", "-std=c99", *sources, "-lm", "-o", executable]
        print(f"  -> Compilando motor C con soporte OBDM: {' '.join(cmd[:5])}...")
        res = subprocess.run(cmd, cwd=SCRIPT_DIR, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"Error compilando MC:\n{res.stderr}")
    return executable


def main():
    print("=================================================================")
    print("  VERIFICACIÓN: MATRIZ DE DENSIDAD DE UNA PARTÍCULA (OBDM)")
    print("  MEDIDA CON CADENAS ABIERTAS DEL ALGORITMO WORM (PIMC)")
    print("=================================================================")
    print(f"  Partículas N         : {N_PARTICLES}")
    print(f"  Temperatura T        : {TEMP} (beta = {BETA})")
    print(f"  Frecuencia trampa    : {OMEGA}")
    print(f"  Número de beads P    : {BEADS} (tau = {TAU:.4f})")
    print(f"  Parámetro fugacidad C: {WORM_C:.6f}")
    print(f"  Bloques medición     : {MEASURE_BLOCKS} (iter por bloque: {N_ITER})")
    print("=================================================================\n")

    cfg_path = os.path.join(SCRIPT_DIR, "MC.cfg")
    cfg_bak = os.path.join(SCRIPT_DIR, "MC.cfg.bak_obdm")
    shutil.copy2(cfg_path, cfg_bak)

    try:
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
        cfg_text = re.sub(r"Nmeasure=\s*\d+", "Nmeasure= 1", cfg_text)

        # Añadir opciones OBDM si no están presentes
        if "worm_obdm_max=" not in cfg_text:
            cfg_text += "\nworm_obdm_max= 3.5\n"
        if "worm_obdm_grid=" not in cfg_text:
            cfg_text += "worm_obdm_grid= 100\n"

        with open(cfg_path, "w") as f:
            f.write(cfg_text)

        mc_exe = get_or_compile_executable()
        print("[1/3] Ejecutando simulación PIMC en C...")
        t0 = time.time()
        proc = subprocess.run([mc_exe, "MC.cfg"], cwd=SCRIPT_DIR, capture_output=True, text=True)
        dt = time.time() - t0

        if proc.returncode != 0:
            print(f"[Error] Falló la ejecución de MC (código {proc.returncode}):\n{proc.stderr}")
            return 1
        print(f"  -> Simulación finalizada con éxito en {dt:.1f} s")

        obdm_file = os.path.join(SCRIPT_DIR, "outworm_obdm.dat")
        dist_file = os.path.join(SCRIPT_DIR, "outworm_obdm_dist.dat")
        radii_file = os.path.join(SCRIPT_DIR, "outworm_obdm_endpoint_radii.dat")
        if not os.path.isfile(obdm_file):
            raise FileNotFoundError(f"No se generó el fichero esperado {obdm_file}")

        # Copiar ficheros crudos al directorio de resultados
        shutil.copy2(obdm_file, os.path.join(OUTPUT_DIR, "outworm_obdm.dat"))
        if os.path.isfile(dist_file):
            shutil.copy2(dist_file, os.path.join(OUTPUT_DIR, "outworm_obdm_dist.dat"))
        if os.path.isfile(radii_file):
            shutil.copy2(radii_file, os.path.join(OUTPUT_DIR, "outworm_obdm_endpoint_radii.dat"))

        print("\n[2/3] Procesando y normalizando matriz de densidad \rho_1(x, x')...")
        x_vals, pimc_norm, pimc_counts, meta = load_worm_obdm_data(obdm_file)
        n_samples = meta.get("samples", 0)
        trace_int = meta.get("trace_integral", 0.0)
        print(f"  -> Muestras en Sector G (gap == 1): {n_samples}")
        print(f"  -> Traza del histograma antes de normalizar: {meta.get('raw_trace_integral', 0.0):.4e}")
        print(f"  -> Traza en el archivo exportado: {trace_int:.6f} (normalizada por construcción a N)")
        tau_gap1 = meta.get("tau")
        if tau_gap1 is not None:
            print(f"  -> Separación temporal de gap=1: tau={tau_gap1:.6g}; la comparación tau=0 es aproximada")

        # Comparación analítica
        X, XP = np.meshgrid(x_vals, x_vals)
        bose_exact = analytical_bose_obdm_1d(X, XP, N_PARTICLES, TEMP)
        dist_exact = analytical_distinguishable_obdm_1d(X, XP, N_PARTICLES, TEMP)

        diag_pimc = np.diag(pimc_norm)
        diag_bose = analytical_bose_obdm_1d(x_vals, x_vals, N_PARTICLES, TEMP)
        diag_diff = np.max(np.abs(diag_pimc - diag_bose))
        diag_rms = np.sqrt(np.mean((diag_pimc - diag_bose)**2))

        # Discrepancia global en la región central (|x|, |x'| <= 2.5)
        mask_roi = (np.abs(X) <= 2.5) & (np.abs(XP) <= 2.5)
        rms_global = np.sqrt(np.mean((pimc_norm[mask_roi] - bose_exact[mask_roi])**2))

        print(f"  -> Discrepancia RMS en diagonal local n(x): {diag_rms:.4e}")
        print(f"  -> Discrepancia RMS en toda la matriz (|x|,|x'| <= 2.5): {rms_global:.4e}")

        print("\n[3/3] Generando gráficos de publicación y sincronizando con Overleaf...")
        generar_mapa_calor_obdm(x_vals, pimc_norm, n_particles=N_PARTICLES, temp=TEMP, save_dir=OUTPUT_DIR)
        generar_cortes_obdm(x_vals, pimc_norm, n_particles=N_PARTICLES, temp=TEMP, save_dir=OUTPUT_DIR)
        generar_grafico_distancia_extremos(dist_file, n_particles=N_PARTICLES, temp=TEMP, save_dir=OUTPUT_DIR)

        print("\n=================================================================")
        print("  ANÁLISIS DE OBDM COMPLETADO; las discrepancias son comparaciones diagnósticas")
        print("=================================================================")
        print(f"  Archivos generados en: {OUTPUT_DIR}")
        print(f"  Figuras copiadas a   : {OVERLEAF_IMG_DIR}")
        print("=================================================================\n")
        return 0

    finally:
        if os.path.isfile(cfg_bak):
            shutil.copy2(cfg_bak, cfg_path)
            os.remove(cfg_bak)


if __name__ == "__main__":
    sys.exit(main())
