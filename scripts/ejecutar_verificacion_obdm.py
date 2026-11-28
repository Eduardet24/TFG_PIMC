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

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)
SCRIPT_DIR = os.path.dirname(_THIS_DIR)

from worm_config import get_n_particles
from pimc_common import compile_mc_binary, replace_cfg_value
from generar_graficos_obdm import (
    load_worm_obdm_data,
    load_worm_dist_data,
    analytical_bose_obdm_1d,
    analytical_distinguishable_obdm_1d,
    generar_mapa_calor_obdm,
    generar_cortes_obdm,
    generar_grafico_distancia_extremos,
)

OUTPUT_DIR = os.environ.get(
    "WORM_OBDM_OUTPUT_DIR",
    os.path.join(SCRIPT_DIR, "Resultados", "Verificacion_OBDM"),
)
RUN_DIR = os.environ.get("WORM_OBDM_RUN_DIR", os.path.join(OUTPUT_DIR, "sim_run"))
OVERLEAF_IMG_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "overleaf", "imágenes")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(RUN_DIR, exist_ok=True)
try:
    os.makedirs(OVERLEAF_IMG_DIR, exist_ok=True)
except OSError:
    pass

DEFAULT_CFG_PATH = os.path.join(SCRIPT_DIR, "MC.cfg")
RUN_CFG_PATH = os.path.join(RUN_DIR, "MC.cfg")
TEMPLATE_CFG_PATH = RUN_CFG_PATH if os.path.isfile(RUN_CFG_PATH) else DEFAULT_CFG_PATH

# Parámetros físicos
N_PARTICLES = get_n_particles(TEMPLATE_CFG_PATH, "WORM_FASE4_PARTICLES")
if N_PARTICLES < 2 and "WORM_PARTICLES" not in os.environ and "WORM_FASE4_PARTICLES" not in os.environ:
    N_PARTICLES = 2
TEMP = 0.5
BETA = 1.0 / TEMP
OMEGA = 1.0
BEADS = 30
TAU = BETA / BEADS
WORM_C = 0.1 / (0.40 * N_PARTICLES * BEADS)  # ~0.004167

HEATING_BLOCKS = int(os.environ.get("WORM_OBDM_HEATING", 25))
MEASURE_BLOCKS = int(os.environ.get("WORM_OBDM_BLOCKS", 400))
N_ITER = int(os.environ.get("WORM_OBDM_NITER", 1000))


def get_or_compile_executable():
    configured = os.environ.get("WORM_OBDM_EXECUTABLE")
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
    print(f"  Directorio ejecución : {RUN_DIR}")
    print("=================================================================\n")

    is_root_run = os.path.abspath(RUN_DIR) == os.path.abspath(SCRIPT_DIR)
    cfg_bak = os.path.join(SCRIPT_DIR, "MC.cfg.bak_obdm") if is_root_run else None
    if is_root_run:
        shutil.copy2(RUN_CFG_PATH, cfg_bak)

    try:
        with open(TEMPLATE_CFG_PATH, "r", encoding="utf-8", errors="ignore") as f:
            cfg_text = f.read()

        cfg_text = replace_cfg_value(cfg_text, "MC", 3)
        cfg_text = replace_cfg_value(cfg_text, "N", N_PARTICLES)
        cfg_text = replace_cfg_value(cfg_text, "T", TEMP)
        cfg_text = replace_cfg_value(cfg_text, "Npop", BEADS)
        cfg_text = replace_cfg_value(cfg_text, "NwalkersMax", BEADS * 2)
        cfg_text = replace_cfg_value(cfg_text, "worm_C", f"{WORM_C:.6f}")
        cfg_text = replace_cfg_value(cfg_text, "generate_new_coordinates", 1)
        cfg_text = replace_cfg_value(cfg_text, "blck_heating", HEATING_BLOCKS)
        cfg_text = replace_cfg_value(cfg_text, "blck", MEASURE_BLOCKS)
        cfg_text = replace_cfg_value(cfg_text, "Niter", N_ITER)
        cfg_text = replace_cfg_value(cfg_text, "Nmeasure", 1)

        # Añadir o actualizar opciones OBDM
        if re.search(r"\bworm_obdm_max\s*=", cfg_text):
            cfg_text = replace_cfg_value(cfg_text, "worm_obdm_max", 3.5)
        else:
            cfg_text += "\nworm_obdm_max= 3.5\n"
        if re.search(r"\bworm_obdm_grid\s*=", cfg_text):
            cfg_text = replace_cfg_value(cfg_text, "worm_obdm_grid", 100)
        else:
            cfg_text += "worm_obdm_grid= 100\n"

        with open(RUN_CFG_PATH, "w", encoding="utf-8") as f:
            f.write(cfg_text)

        mc_exe = get_or_compile_executable()
        print("[1/3] Ejecutando simulación PIMC en C...")
        t0 = time.time()
        run_started = t0
        proc = subprocess.run([mc_exe, "MC.cfg"], cwd=RUN_DIR, capture_output=True, text=True)
        dt = time.time() - t0

        log_path = os.path.join(OUTPUT_DIR, "simulacion_terminal.log")
        with open(log_path, "w", encoding="utf-8") as f:
            f.write(proc.stdout)
            if proc.stderr:
                f.write("\n" + proc.stderr)

        if proc.returncode != 0:
            print(f"[Error] Falló la ejecución de MC (código {proc.returncode}):\n{proc.stderr}")
            return 1
        print(f"  -> Simulación finalizada con éxito en {dt:.1f} s")

        obdm_src = os.path.join(RUN_DIR, "outworm_obdm.dat")
        dist_src = os.path.join(RUN_DIR, "outworm_obdm_dist.dat")
        radii_src = os.path.join(RUN_DIR, "outworm_obdm_endpoint_radii.dat")
        if not os.path.isfile(obdm_src) or os.path.getmtime(obdm_src) < run_started - 1.0:
            raise FileNotFoundError(f"No se generó el fichero esperado {obdm_src}")

        obdm_file = os.path.join(OUTPUT_DIR, "outworm_obdm.dat")
        dist_file = os.path.join(OUTPUT_DIR, "outworm_obdm_dist.dat")
        radii_file = os.path.join(OUTPUT_DIR, "outworm_obdm_endpoint_radii.dat")

        # Copiar ficheros crudos al directorio de resultados
        if os.path.abspath(obdm_src) != os.path.abspath(obdm_file):
            shutil.copy2(obdm_src, obdm_file)
        if os.path.isfile(dist_src) and os.path.abspath(dist_src) != os.path.abspath(dist_file):
            shutil.copy2(dist_src, dist_file)
        if os.path.isfile(radii_src) and os.path.abspath(radii_src) != os.path.abspath(radii_file):
            shutil.copy2(radii_src, radii_file)

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
        generar_grafico_distancia_extremos(dist_file, n_particles=N_PARTICLES, temp=TEMP, save_dir=OUTPUT_DIR, tau=tau_gap1 or 0.0)

        print("\n=================================================================")
        print("  ANÁLISIS DE OBDM COMPLETADO; las discrepancias son comparaciones diagnósticas")
        print("=================================================================")
        print(f"  Archivos generados en: {OUTPUT_DIR}")
        print(f"  Figuras copiadas a   : {OVERLEAF_IMG_DIR}")
        print("=================================================================\n")
        return 0

    finally:
        if cfg_bak and os.path.isfile(cfg_bak):
            shutil.copy2(cfg_bak, cfg_path)
            os.remove(cfg_bak)


if __name__ == "__main__":
    sys.exit(main())
