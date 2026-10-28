"""Comparación de PIMC (N > 1) en Oscilador Armónico 2D con Resultados Analíticos.

Compara la simulación PIMC en C (cadenas cerradas, partículas distinguibles / boltzmannons)
con las soluciones analíticas exactas para:
  1. Partículas distinguibles en 2D: E/N = coth(beta/2), n(r) gaussiana.
  2. Bosones idénticos en 2D (Bose-Einstein): cálculo canónico/macrocanónico
     con resolución del potencial químico mu(T, N) y expansión en ciclos de permutación.

Genera:
  - Grafica 1: E/N vs Temperatura (PIMC vs Distinguibles vs Bosones + Fraccion condensada).
  - Grafica 2: Perfil de densidad radial n(r) (PIMC vs Distinguibles vs Bosones).
  - Tabla cuantitativa de discrepancias para justificar la necesidad del Worm Algorithm.
"""

import argparse
from pathlib import Path
import re
import shutil
import subprocess

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR
RESULTS_BASE = PROJECT_DIR / 'Resultados' / 'NParticulas2D'
BINARY_MC = PROJECT_DIR / 'MC'
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


def solve_bose_einstein(n_particles, temp, omega=OMEGA, max_terms=250):
    """Resuelve la estadística de Bose-Einstein para N bosones en oscilador 2D."""
    beta = 1.0 / temp

    # Resuelve el parámetro w = exp(beta*(mu - epsilon_0)) in (0, 1)
    # tal que sum_{k=1}^infty w^k / (1 - exp(-k * beta * omega))^2 = N
    def f_w(w):
        val = 0.0
        for k in range(1, max_terms):
            denom = (1.0 - np.exp(-k * beta * omega))**2
            term = (w**k) / denom
            val += term
            if term < 1e-13:
                break
        return val - n_particles

    w_min = 1e-15
    w_max = 1.0 - 1e-11
    w = brentq(f_w, w_min, w_max)

    # Cálculo de energía total por partícula
    # E = sum_{k=1}^infty w^k * [omega * (1 + exp(-k*beta*omega)) / (1 - exp(-k*beta*omega))^3]
    e_tot = 0.0
    for k in range(1, max_terms):
        q = np.exp(-k * beta * omega)
        term = (w**k) * omega * (1.0 + q) / ((1.0 - q)**3)
        e_tot += term
        if term < 1e-13:
            break
    e_per_particle = e_tot / n_particles

    # Fracción de condensado en el estado fundamental (n=0):
    # N0 = w / (1 - w)
    n0 = w / (1.0 - w)
    n0_fraction = min(max(n0 / n_particles, 0.0), 1.0)

    return w, e_per_particle, n0_fraction


def analytical_bose_density(r_grid, n_particles, temp, omega=OMEGA, max_terms=250):
    """Densidad radial exacta de Bose-Einstein mediante expansión en ciclos de permutación."""
    w, _, _ = solve_bose_einstein(n_particles, temp, omega, max_terms)
    beta = 1.0 / temp

    density = np.zeros_like(r_grid)
    for k in range(1, max_terms):
        weight = (w**k) / ((1.0 - np.exp(-k * beta * omega))**2)
        sigma_k2 = (0.5 / omega) / np.tanh(0.5 * k * beta * omega)
        gauss = (1.0 / (2.0 * np.pi * sigma_k2)) * np.exp(-r_grid**2 / (2.0 * sigma_k2))
        density += weight * gauss
        if weight < 1e-11:
            break

    return density / n_particles


# ==============================================================================
# Ejecución y Lectura de PIMC (C Engine)
# ==============================================================================

def run_simulation(n_particles, temp, beads, run_dir, blck_heating=20, blck=60, niter=200):
    """Ejecuta una corrida del ejecutable MC en C en un directorio aislado."""
    run_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(BINARY_MC, run_dir / 'MC')
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


def read_energy_results(run_dir):
    """Lee y promedia la energía de oute.dat (# Etot Epot Ekin Ekin_prim)."""
    energy_file = run_dir / 'oute.dat'
    if not energy_file.is_file():
        raise FileNotFoundError(f'No se encontró {energy_file}')

    data = np.loadtxt(energy_file, skiprows=1)
    if data.ndim == 1:
        data = data[np.newaxis, :]

    # Descarta la primera mitad de mediciones para asegurar termalización
    n_used = max(len(data) // 2, 1)
    useful = data[-n_used:]

    e_tot_mean = np.mean(useful[:, 0])
    e_tot_err = np.std(useful[:, 0]) / np.sqrt(n_used)
    e_pot_mean = np.mean(useful[:, 1])
    e_pot_err = np.std(useful[:, 1]) / np.sqrt(n_used)
    e_kin_vir_mean = np.mean(useful[:, 2])
    e_kin_vir_err = np.std(useful[:, 2]) / np.sqrt(n_used)

    return {
        'etot': e_tot_mean,
        'etot_err': e_tot_err,
        'epot': e_pot_mean,
        'epot_err': e_pot_err,
        'ekin_vir': e_kin_vir_mean,
        'ekin_vir_err': e_kin_vir_err,
    }


def read_density_results(run_dir):
    """Lee el perfil acumulativo final de outrd.dat."""
    rd_file = run_dir / 'outrd.dat'
    if not rd_file.is_file():
        raise FileNotFoundError(f'No se encontró {rd_file}')

    data = np.loadtxt(rd_file)
    saltos = np.where(np.diff(data[:, 0]) < 0)[0] + 1
    bloques = np.split(data, saltos)
    last_block = bloques[-1]
    r = last_block[:, 0]
    n_sim = last_block[:, 1]
    return r, n_sim


# ==============================================================================
# Gráficas y Comparación
# ==============================================================================

def plot_energy_comparison(n_particles, temps, sim_energies, out_path):
    """Genera la figura de Energía vs Temperatura comparando PIMC, Distinguibles y Bosones."""
    t_fine = np.linspace(min(temps) * 0.8, max(temps) * 1.05, 120)
    e_dist_fine = [analytical_distinguishable_energy(t) for t in t_fine]
    e_bose_fine = [solve_bose_einstein(n_particles, t)[1] for t in t_fine]
    n0_fine = [solve_bose_einstein(n_particles, t)[2] for t in t_fine]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7.5), sharex=True, gridspec_kw={'height_ratios': [2, 1]})

    # Panel 1: Energía por partícula
    ax1.plot(t_fine, e_dist_fine, 'k-', linewidth=2.0, label='Analítico: Distinguibles (Boltzmannons)')
    ax1.plot(t_fine, e_bose_fine, 'tab:red', linestyle='--', linewidth=2.0,
             label=f'Analítico: Bosones Idénticos ($N={n_particles}$)')

    t_sim = [res['temp'] for res in sim_energies]
    e_sim = [res['etot'] for res in sim_energies]
    e_sim_err = [res['etot_err'] for res in sim_energies]
    ax1.errorbar(t_sim, e_sim, yerr=e_sim_err, fmt='o', color='tab:blue', ecolor='tab:blue',
                 elinewidth=1.5, capsize=4, markersize=6.5, label='PIMC (Simulación C actual, $N=%d$)' % n_particles)

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


def plot_density_comparison(n_particles, selected_temp, r_sim, n_sim, out_path):
    """Genera la figura de densidad radial n(r) para una temperatura dada."""
    r_fine = np.linspace(0.001, max(r_sim), 200)
    n_dist_fine = analytical_distinguishable_density(r_fine, selected_temp)
    n_bose_fine = analytical_bose_density(r_fine, n_particles, selected_temp)

    fig, ax = plt.subplots(figsize=(7.5, 5))

    ax.plot(r_fine, n_dist_fine, 'k-', linewidth=2.0, label='Analítico: Distinguibles (Gaussiana)')
    ax.plot(r_fine, n_bose_fine, 'tab:red', linestyle='--', linewidth=2.0,
            label=f'Analítico: Bosones de Bose ($N={n_particles}$)')
    ax.plot(r_sim, n_sim, 'o', color='tab:blue', markersize=4.5, alpha=0.85,
            label='PIMC (Simulación C actual)')

    ax.set_xlabel(r'Distancia radial $r = \sqrt{x^2 + y^2}$ [$a_{\mathrm{ho}}$]', fontsize=11)
    ax.set_ylabel(r'Densidad radial por partícula $n(r)$', fontsize=11)
    ax.set_title(f'Perfil de Densidad $n(r)$ en 2D ($N = {n_particles}$, $T = {selected_temp:g}$)',
                 fontsize=12, fontweight='bold', pad=10)
    ax.grid(True, linestyle=':', alpha=0.5)
    ax.legend(loc='upper right', fontsize=10, framealpha=0.9)

    # Nota explicativa en recuadro
    text_info = (
        'Efecto del Worm Algorithm:\n'
        'El código PIMC actual muestrea boltzmannons\n'
        '(coincide exactamente con la curva negra).\n'
        'Al implementar el Worm (permutaciones),\n'
        'la simulación transicionará hacia la curva roja.'
    )
    ax.text(0.04, 0.06, text_info, transform=ax.transAxes, fontsize=9.0,
            verticalalignment='bottom', bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.9, edgecolor='gray'))

    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f'-> Figura de densidad radial guardada en:\n   {out_path}')


def print_comparison_table(n_particles, sim_energies):
    """Imprime una tabla detallada comparando los resultados numéricos y analíticos."""
    print('\n' + '=' * 82)
    print(f'TABLA COMPARATIVA: $N = {n_particles}$ PARTÍCULAS EN OSCILADOR ARMÓNICO 2D')
    print('=' * 82)
    print(f"{'T':>6} | {'E/N PIMC':>14} | {'E/N Disting':>12} | {'E/N Bose':>10} | {'Discrep. PIMC-Bose':>18} | {'N0/N (Cond)':>10}")
    print('-' * 82)

    for item in sim_energies:
        temp = item['temp']
        e_pimc = item['etot']
        e_dist = analytical_distinguishable_energy(temp)
        _, e_bose, n0 = solve_bose_einstein(n_particles, temp)
        diff_pct = (e_pimc - e_bose) / e_pimc * 100.0
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
    parser.add_argument('--run', action='store_true', help='Ejecutar las simulaciones en C')
    args = parser.parse_args()

    temps = [float(t.strip()) for t in args.temperatures.split(',')]
    n_dir = RESULTS_BASE / f'N{args.particles}'
    n_dir.mkdir(parents=True, exist_ok=True)

    sim_energies = []
    r_density_sample = None
    n_density_sample = None

    for temp in temps:
        sim_dir = n_dir / f'sim_T{temp:g}'
        if args.run:
            print(f'-> Simulando N={args.particles}, T={temp:g}, beads={args.beads}...')
            run_simulation(
                n_particles=args.particles,
                temp=temp,
                beads=args.beads,
                run_dir=sim_dir,
                blck_heating=20,
                blck=args.blocks,
                niter=args.iterations,
            )

        # Leer resultados
        res = read_energy_results(sim_dir)
        res['temp'] = temp
        sim_energies.append(res)

        if np.isclose(temp, args.temp_density):
            r_density_sample, n_density_sample = read_density_results(sim_dir)

    # Imprimir tabla
    print_comparison_table(args.particles, sim_energies)

    # Guardar figuras
    energy_fig = n_dir / f'energia_vs_T_N{args.particles}.png'
    plot_energy_comparison(args.particles, temps, sim_energies, energy_fig)

    if r_density_sample is not None and n_density_sample is not None:
        density_fig = n_dir / f'densidad_radial_N{args.particles}_T{args.temp_density:g}.png'
        plot_density_comparison(args.particles, args.temp_density, r_density_sample, n_density_sample, density_fig)


if __name__ == '__main__':
    main()

