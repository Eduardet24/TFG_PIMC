"""Utilidades compartidas para modelos analíticos 2D, estadística de bloques y compilación C."""

import os
import re
import subprocess
from pathlib import Path
import numpy as np

DEFAULT_OMEGA = 1.0


def canonical_bose_partition(n_particles, temp, omega=DEFAULT_OMEGA):
    """Calcula las sumas de partición canónicas exactas para bosones en oscilador armónico 2D."""
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
            one_cycle_partition = 1.0 / (1.0 - q) ** 2
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


def analytical_bose_energy(n_particles, temp, omega=DEFAULT_OMEGA):
    """Energía canónica exacta por partícula y fracción de condensado N0/N en 2D."""
    partition, excitation_energy_sum, excited_partition = canonical_bose_partition(
        n_particles, temp, omega
    )
    e_per_particle = omega + excitation_energy_sum[n_particles] / (
        n_particles * partition[n_particles]
    )
    ground_state_occupancy = sum(
        (n_particles - count) * excited_partition[count]
        for count in range(n_particles)
    )
    ground_state_occupancy /= partition[n_particles]
    return e_per_particle, ground_state_occupancy / n_particles


def analytical_distinguishable_energy(temp, omega=DEFAULT_OMEGA):
    """Energía por partícula para partículas distinguibles en oscilador armónico 2D."""
    beta = 1.0 / temp
    return omega / np.tanh(0.5 * beta * omega)


def analytical_bose_cycle_fractions(n_particles, temp, omega=DEFAULT_OMEGA):
    """Fracción canónica de partículas en ciclos de permutación de longitud L."""
    partition, _, _ = canonical_bose_partition(n_particles, temp, omega)
    beta = 1.0 / temp
    return np.array([
        (1.0 / (1.0 - np.exp(-cycle * beta * omega)) ** 2)
        * partition[n_particles - cycle]
        / (n_particles * partition[n_particles])
        for cycle in range(1, n_particles + 1)
    ])


def analytical_bose_density(r_grid, n_particles, temp, omega=DEFAULT_OMEGA):
    """Densidad radial canónica por partícula para N bosones ideales en 2D."""
    partition, _, _ = canonical_bose_partition(n_particles, temp, omega)
    beta = 1.0 / temp
    normalization = n_particles * partition[n_particles]

    density = np.zeros_like(r_grid)
    for cycle in range(1, n_particles + 1):
        q = np.exp(-cycle * beta * omega)
        cycle_partition = 1.0 / (1.0 - q) ** 2
        weight = cycle_partition * partition[n_particles - cycle] / normalization
        sigma_k2 = (0.5 / omega) / np.tanh(0.5 * cycle * beta * omega)
        gauss = (1.0 / (2.0 * np.pi * sigma_k2)) * np.exp(-r_grid ** 2 / (2.0 * sigma_k2))
        density += weight * gauss
    return density


def analytical_distinguishable_density(r_grid, temp, omega=DEFAULT_OMEGA):
    """Densidad radial n(r) para partículas distinguibles (gaussiana térmica 2D)."""
    beta = 1.0 / temp
    sigma2 = (0.5 / omega) / np.tanh(0.5 * beta * omega)
    return (1.0 / (2.0 * np.pi * sigma2)) * np.exp(-r_grid ** 2 / (2.0 * sigma2))


def blocked_ratio_stats(values, counts):
    """Media ponderada y error estándar asintótico para estimador de razón por bloques."""
    values = np.asarray(values, dtype=float)
    counts = np.asarray(counts, dtype=float)
    if len(values) < 2 or np.any(counts < 0) or counts.sum() <= 0:
        raise ValueError("Se necesitan al menos dos bloques válidos para la estadística")
    mean = np.average(values, weights=counts)
    residuals = counts * (values - mean)
    error = np.std(residuals, ddof=1) / (np.sqrt(len(residuals)) * counts.mean())
    return mean, error


def compile_mc_binary(project_dir, output_path, only_if_stale=False, compiler=None):
    """Compila los fuentes C del repositorio en output_path usando gcc/CC."""
    project_path = Path(project_dir)
    output = Path(output_path)
    sources = sorted(project_path.glob("*.c"))
    if not sources:
        raise RuntimeError(f"No se encontraron fuentes C en {project_path}")

    if only_if_stale and output.is_file():
        headers = list(project_path.glob("*.h"))
        newest_src = max(p.stat().st_mtime for p in [*sources, *headers])
        if output.stat().st_mtime >= newest_src:
            return output

    cc = compiler or os.environ.get("CC", "gcc")
    cmd = [cc, "-O2", "-fopenmp", "-std=c99", *(str(s) for s in sources), "-lm", "-o", str(output)]
    result = subprocess.run(cmd, cwd=project_path, capture_output=True, text=True)
    if result.returncode != 0:
        detail = result.stderr or result.stdout
        raise RuntimeError(f"Error compilando {output.name}:\n{detail[-6000:]}")
    return output


_CFG_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


def read_cfg_value(config_text, key):
    """Lee un parámetro numérico `key= valor` de un texto MC.cfg."""
    match = re.search(rf"\b{re.escape(key)}\s*=\s*({_CFG_NUMBER})", config_text)
    if match is None:
        raise ValueError(f"No se pudo leer {key} de MC.cfg")
    return float(match.group(1))


def replace_cfg_value(config_text, key, value):
    """Sustituye un parámetro numérico `key= valor` exactamente una vez en MC.cfg."""
    pattern = rf"(\b{re.escape(key)}\s*=\s*){_CFG_NUMBER}"
    updated, replacements = re.subn(
        pattern, lambda match: f"{match.group(1)}{value}", config_text, count=1
    )
    if replacements != 1:
        raise ValueError(f"No se pudo sustituir {key} en MC.cfg")
    return updated


def replace_cfg_line(config_text, key, new_val):
    """Sustituye el valor (numérico o texto) de `key=` en una línea de MC.cfg."""
    pattern = rf"(\b{re.escape(key)}\s*=\s*)[^\s³]+"
    return re.sub(pattern, rf"\g<1>{new_val}", config_text)
