"""Benchmark PIMC de N partículas ideales en un oscilador armónico 2D.

Compara la simulación PIMC en C (cadenas cerradas, partículas distinguibles / boltzmannons)
con las soluciones analíticas exactas para:
  1. Partículas distinguibles en 2D: E/N = coth(beta/2), n(r) gaussiana.
    2. Bosones idénticos en 2D: solución canónica exacta a N fijo.

Genera:
  - Grafica 1: E/N vs Temperatura (PIMC vs Distinguibles vs Bosones + Fraccion condensada).
  - Grafica 2: Perfil de densidad radial n(r) (PIMC vs Distinguibles vs Bosones).
    - Tabla de discrepancias: el PIMC actual no incluye permutaciones bosónicas.
"""

import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import matplotlib.pyplot as plt
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR
RESULTS_BASE = PROJECT_DIR / 'Resultados' / 'NParticulas2D'
TEMPLATE_CFG = PROJECT_DIR / 'MC.cfg'

DEFAULT_N = 5
DEFAULT_BEADS = 30
DEFAULT_TEMPS = (0.2, 0.5, 1.0, 1.5, 2.0, 3.0)
OMEGA = 1.0


# ==============================================================================
# Modelos Analíticos Exactos
# ==============================================================================

def analytical_distinguishable_energy(temp, omega=OMEGA):
    """Energía por partícula para oscilador armónico 2D (partículas distinguibles)."""
    beta = 1.0 / temp
    return omega / np.tanh(0.5 * beta * omega)


def analytical_distinguishable_density(r_grid, temp, omega=OMEGA):
    """Densidad radial n(r) para partículas distinguibles (gaussiana térmica)."""
    beta = 1.0 / temp
    sigma2 = (0.5 / omega) / np.tanh(0.5 * beta * omega)
    return (1.0 / (2.0 * np.pi * sigma2)) * np.exp(-r_grid**2 / (2.0 * sigma2))


def canonical_bose_partition(n_particles, temp, omega=OMEGA):
    """Return canonical partition sums and excitation-energy-weighted sums."""
    if n_particles < 1 or temp <= 0 or omega <= 0:
        raise ValueError('N, T y omega deben ser positivos')

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
    """Canonical energy per particle for a fixed number of 2D ideal bosons."""
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


def analytical_bose_density(r_grid, n_particles, temp, omega=OMEGA):
    """Canonical per-particle 2D density for a fixed number of ideal bosons."""
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


# ==============================================================================
# Ejecución y Lectura de PIMC (C Engine)
# ==============================================================================

def build_mc_binary(output_path):
    """Compile the current C sources into an isolated executable."""
    sources = sorted(SCRIPT_DIR.glob('*.c'))
    command = [os.environ.get('CC', 'gcc'), '-O2', '-fopenmp', '-std=c99']
    command.extend(str(source) for source in sources)
    command.extend(['-lm', '-o', str(output_path)])
    result = subprocess.run(command, cwd=SCRIPT_DIR, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f'No se pudo compilar MC:\n{result.stdout}\n{result.stderr}')
    return output_path


def run_simulation(binary_path, n_particles, temp, beads, run_dir,
                   blck_heating=20, blck=60, niter=200):
    """Ejecuta una corrida del ejecutable MC en C en un directorio aislado."""
    run_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(binary_path, run_dir / 'MC')
    cfg_text = TEMPLATE_CFG.read_text(encoding='utf-8', errors='ignore')

    def replace_line(text, key, new_val):
        pattern = rf'(\b{key}\s*=\s*)[^\s³]+'
        return re.sub(pattern, rf'\g<1>{new_val}', text)

    cfg_text = replace_line(cfg_text, 'MC', '3')
    cfg_text = replace_line(cfg_text, 'N', str(n_particles))
    cfg_text = replace_line(cfg_text, 'T', f'{temp:g}')
    cfg_text = replace_line(cfg_text, 'Npop', str(beads))
    cfg_text = replace_line(cfg_text, 'NwalkersMax', str(beads + 20))
    cfg_text = replace_line(cfg_text, 'omega_x', f'{OMEGA:g}')
    cfg_text = replace_line(cfg_text, 'omega_y', f'{OMEGA:g}')
    cfg_text = replace_line(cfg_text, 'alpha_x', '1')
    cfg_text = replace_line(cfg_text, 'alpha_y', '1')
    cfg_text = replace_line(cfg_text, 'blck_heating', str(blck_heating))
    cfg_text = replace_line(cfg_text, 'blck', str(blck))
    cfg_text = replace_line(cfg_text, 'Niter', str(niter))
    cfg_text = replace_line(cfg_text, 'gridRD', '60')
    cfg_text = replace_line(cfg_text, 'RDmax', '3.5')
    cfg_text = replace_line(cfg_text, 'RDwidth', '3.5')

    (run_dir / 'MC.cfg').write_text(cfg_text, encoding='utf-8')

    log_file = run_dir / 'run.log'
    with log_file.open('w', encoding='utf-8') as f_log:
        res = subprocess.run(['./MC'], cwd=run_dir, stdout=f_log, stderr=subprocess.STDOUT, text=True)

    if res.returncode != 0:
        raise RuntimeError(f'Error en simulación MC en {run_dir}. Revisa {log_file}')


def validate_run_output(run_dir):
    """Reject archived runs produced by the pre-benchmark executable."""
    log_file = run_dir / 'run.log'
    if not log_file.is_file():
        raise FileNotFoundError(f'No se encontró {log_file}; ejecuta con --run')
    log = log_file.read_text(encoding='utf-8', errors='ignore')
    expected_mode = 'exact harmonic normal modes (all noninteracting particles)'
    if expected_mode not in log or 'Interaction potential (not included in Eloc part): ABSENT' not in log:
        raise RuntimeError(
            f'{run_dir} no contiene una corrida compatible con este benchmark. '
            'Vuelve a ejecutarlo con --run para compilar y simular desde los fuentes actuales.'
        )


def read_energy_results(run_dir):
    """Lee y promedia la energía de oute.dat (# Etot Epot Ekin Ekin_prim)."""
    energy_file = run_dir / 'oute.dat'
    if not energy_file.is_file():
        raise FileNotFoundError(f'No se encontró {energy_file}')

    data = np.loadtxt(energy_file, skiprows=1)
    if data.ndim == 1:
        data = data[np.newaxis, :]

    useful_start = len(data) // 2
    useful = data[useful_start:]
    config = (run_dir / 'MC.cfg').read_text(encoding='utf-8', errors='ignore')
    niter_match = re.search(r'\bNiter=\s*(\d+)', config)
    nmeasure_match = re.search(r'\bNmeasure=\s*(\d+)', config)
    niter = int(niter_match.group(1)) if niter_match else len(data)
    nmeasure = int(nmeasure_match.group(1)) if nmeasure_match else 1
    block_size = max(niter // nmeasure, 1)
    first_full_block = ((useful_start + block_size - 1) // block_size) * block_size
    block_means = np.array([
        data[start:start + block_size].mean(axis=0)
        for start in range(first_full_block, len(data) - block_size + 1, block_size)
    ])

    e_tot_mean = np.mean(useful[:, 0])
    block_errors = (np.std(block_means, axis=0, ddof=1) / np.sqrt(len(block_means))
                    if len(block_means) > 1 else np.full(3, np.nan))
    e_tot_err, e_pot_err, e_kin_vir_err = block_errors[:3]
    e_pot_mean = np.mean(useful[:, 1])
    e_kin_vir_mean = np.mean(useful[:, 2])

    return {
        'etot': e_tot_mean,
        'etot_err': e_tot_err,
        'epot': e_pot_mean,
        'epot_err': e_pot_err,
        'ekin_vir': e_kin_vir_mean,
        'ekin_vir_err': e_kin_vir_err,
    }


def read_density_results(run_dir):
    """Average the latter half of saved density blocks."""
    rd_file = run_dir / 'outrd.dat'
    if not rd_file.is_file():
        raise FileNotFoundError(f'No se encontró {rd_file}')

    data = np.loadtxt(rd_file)
    saltos = np.where(np.diff(data[:, 0]) < 0)[0] + 1
    bloques = np.split(data, saltos)
    selected_blocks = bloques[len(bloques) // 2:]
    r = selected_blocks[-1][:, 0]
    profiles = np.stack([block[:, 1] for block in selected_blocks])
    errors = (profiles.std(axis=0, ddof=1) / np.sqrt(len(profiles))
              if len(profiles) > 1 else np.full(profiles.shape[1], np.nan))
    return r, profiles.mean(axis=0), errors


# ==============================================================================
# Gráficas y Comparación
# ==============================================================================

def plot_energy_comparison(n_particles, temps, sim_energies, out_path):
    """Genera la figura de Energía vs Temperatura comparando PIMC, Distinguibles y Bosones."""
    t_fine = np.linspace(min(temps) * 0.8, max(temps) * 1.05, 120)
    e_dist_fine = [analytical_distinguishable_energy(t) for t in t_fine]
    bose_values = [analytical_bose_energy(n_particles, t) for t in t_fine]
    e_bose_fine = [value[0] for value in bose_values]
    n0_fine = [value[1] for value in bose_values]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7.5), sharex=True, gridspec_kw={'height_ratios': [2, 1]})

    # Panel 1: Energía por partícula
    ax1.plot(t_fine, e_dist_fine, 'k-', linewidth=2.0, label='Analítico: Distinguibles (Boltzmannons)')
    ax1.plot(t_fine, e_bose_fine, 'tab:red', linestyle='--', linewidth=2.0,
             label=f'Analítico: Bosones Idénticos ($N={n_particles}$)')

    t_sim = [res['temp'] for res in sim_energies]
    e_sim = [res['etot'] for res in sim_energies]
    e_sim_err = [res['etot_err'] for res in sim_energies]
    ax1.errorbar(t_sim, e_sim, yerr=e_sim_err, fmt='o', color='tab:blue', ecolor='tab:blue',
                 elinewidth=1.5, capsize=4, markersize=6.5,
                 label='PIMC sin permutaciones ($N=%d$)' % n_particles)

    ax1.set_ylabel(r'Energía por partícula $E / N$ [$\hbar\omega$]', fontsize=11)
    ax1.set_title(f'Energía vs Temperatura en Oscilador Armónico 2D ($N = {n_particles}$)',
                  fontsize=12, fontweight='bold', pad=10)
    ax1.grid(True, linestyle=':', alpha=0.5)
    ax1.legend(loc='upper left', fontsize=10, framealpha=0.9)

    # Panel 2: Fracción de Condensado N0/N y Discrepancia
    ax2.plot(t_fine, n0_fine, color='tab:red', linewidth=1.8, label=r'Fracción Condensada $N_0 / N$ (Bose)')
    diff_fine = (np.array(e_dist_fine) - np.array(e_bose_fine)) / np.array(e_dist_fine) * 100.0
    ax2.plot(t_fine, diff_fine / 100.0, color='tab:purple', linestyle=':', linewidth=1.8,
             label=r'Diferencia Relativa $\Delta E / E_{\mathrm{dist}}$')

    ax2.set_xlabel(r'Temperatura $T$ [$k_B / \hbar\omega$]', fontsize=11)
    ax2.set_ylabel(r'Fracción / Discrepancia', fontsize=10.5)
    ax2.set_ylim(-0.05, 1.05)
    ax2.grid(True, linestyle=':', alpha=0.5)
    ax2.legend(loc='upper right', fontsize=9.5, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f'-> Figura de energía guardada en:\n   {out_path}')


def plot_density_comparison(n_particles, selected_temp, r_sim, n_sim, n_sim_err, out_path):
    """Genera la figura de densidad radial n(r) para una temperatura dada."""
    r_fine = np.linspace(0.001, max(r_sim), 200)
    n_dist_fine = analytical_distinguishable_density(r_fine, selected_temp)
    n_bose_fine = analytical_bose_density(r_fine, n_particles, selected_temp)

    fig, ax = plt.subplots(figsize=(7.5, 5))

    ax.plot(r_fine, n_dist_fine, 'k-', linewidth=2.0, label='Analítico: Distinguibles (Gaussiana)')
    ax.plot(r_fine, n_bose_fine, 'tab:red', linestyle='--', linewidth=2.0,
            label=f'Analítico: Bosones (canónico, $N={n_particles}$)')
    ax.errorbar(r_sim, n_sim, yerr=n_sim_err, fmt='o', color='tab:blue', ecolor='tab:blue',
                elinewidth=1.0, capsize=2, markersize=4.0, alpha=0.85,
                label='PIMC sin permutaciones')

    ax.set_xlabel(r'Distancia radial $r = \sqrt{x^2 + y^2}$ [$a_{\mathrm{ho}}$]', fontsize=11)
    ax.set_ylabel(r'Densidad radial por partícula $n(r)$', fontsize=11)
    ax.set_title(f'Perfil de Densidad $n(r)$ en 2D ($N = {n_particles}$, $T = {selected_temp:g}$)',
                 fontsize=12, fontweight='bold', pad=10)
    ax.grid(True, linestyle=':', alpha=0.5)
    ax.legend(loc='upper right', fontsize=10, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f'-> Figura de densidad radial guardada en:\n   {out_path}')


def print_comparison_table(n_particles, sim_energies):
    """Imprime una tabla detallada comparando los resultados numéricos y analíticos."""
    print('\n' + '=' * 82)
    print(f'TABLA COMPARATIVA: $N = {n_particles}$ PARTÍCULAS EN OSCILADOR ARMÓNICO 2D')
    print('=' * 82)
    print(f"{'T':>6} | {'E/N PIMC':>14} | {'E/N Disting':>12} | {'E/N Bose':>10} | {'(PIMC-Bose)/Edist':>18} | {'N0/N (Cond)':>10}")
    print('-' * 82)

    for item in sim_energies:
        temp = item['temp']
        e_pimc = item['etot']
        e_dist = analytical_distinguishable_energy(temp)
        e_bose, n0 = analytical_bose_energy(n_particles, temp)
        diff_pct = (e_pimc - e_bose) / e_dist * 100.0
        print(f"{temp:6.2f} | {e_pimc:8.4f} +/- {item['etot_err']:.4f} | {e_dist:12.4f} | {e_bose:10.4f} | {diff_pct:16.2f}% | {n0:10.4f}")

    print('=' * 82)


# ==============================================================================
# Función Principal
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--particles', type=int, default=DEFAULT_N,
                        help=f'Número de partículas N (por defecto: {DEFAULT_N})')
    parser.add_argument('--beads', type=int, default=DEFAULT_BEADS,
                        help=f'Número de cuentas/beads (por defecto: {DEFAULT_BEADS})')
    parser.add_argument('--temperatures', default='0.2,0.5,1.0,1.5,2.0,3.0',
                        help='Temperaturas separadas por comas')
    parser.add_argument('--temp-density', type=float, default=0.5,
                        help='Temperatura seleccionada para la gráfica de densidad n(r) (por defecto: 0.5)')
    parser.add_argument('--blocks', type=int, default=60, help='Bloques de medida en MC')
    parser.add_argument('--iterations', type=int, default=200, help='Iteraciones por bloque en MC')
    parser.add_argument('--run', action='store_true', help='Compilar el código C actual y ejecutar las simulaciones')
    args = parser.parse_args()

    temps = [float(t.strip()) for t in args.temperatures.split(',')]
    if args.particles < 1 or args.beads < 2 or args.blocks < 1 or args.iterations < 1:
        parser.error('N>=1, beads>=2, blocks>=1 e iterations>=1 son obligatorios')
    if not temps or any(temp <= 0 for temp in temps):
        parser.error('Todas las temperaturas deben ser positivas')

    n_dir = RESULTS_BASE / f'N{args.particles}'
    n_dir.mkdir(parents=True, exist_ok=True)

    if args.run:
        with tempfile.TemporaryDirectory(prefix='pimc-build-') as build_dir:
            binary_path = build_mc_binary(Path(build_dir) / 'MC')
            for temp in temps:
                sim_dir = n_dir / f'sim_T{temp:g}'
                print(f'-> Simulando N={args.particles}, T={temp:g}, beads={args.beads}...')
                run_simulation(
                    binary_path=binary_path,
                    n_particles=args.particles,
                    temp=temp,
                    beads=args.beads,
                    run_dir=sim_dir,
                    blck_heating=20,
                    blck=args.blocks,
                    niter=args.iterations,
                )

    sim_energies = []
    r_density_sample = None
    n_density_sample = None
    n_density_error = None

    for temp in temps:
        sim_dir = n_dir / f'sim_T{temp:g}'
        validate_run_output(sim_dir)
        res = read_energy_results(sim_dir)
        res['temp'] = temp
        sim_energies.append(res)

        if np.isclose(temp, args.temp_density):
            r_density_sample, n_density_sample, n_density_error = read_density_results(sim_dir)

    # Imprimir tabla
    print_comparison_table(args.particles, sim_energies)

    # Guardar figuras
    energy_fig = n_dir / f'energia_vs_T_N{args.particles}.png'
    plot_energy_comparison(args.particles, temps, sim_energies, energy_fig)

    if r_density_sample is not None and n_density_sample is not None:
        density_fig = n_dir / f'densidad_radial_N{args.particles}_T{args.temp_density:g}.png'
        plot_density_comparison(
            args.particles, args.temp_density, r_density_sample, n_density_sample,
            n_density_error, density_fig
        )


if __name__ == '__main__':
    main()

