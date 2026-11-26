"""Visualización de caminos en tiempo imaginario para el oscilador armónico 2D en PIMC.

Ejecuta la simulación de Monte Carlo en C (MC) para el oscilador armónico 2D
a temperaturas altas y bajas, extrae las coordenadas de las cuentas (beads)
del polímero cerrado desde in3Dprev.in y genera la figura comparativa.

Muestra:
  - Curvas de nivel del pozo de potencial 2D: V(x, y) = 1/2 (x^2 + y^2).
  - El camino cerrado (polímero cuántico) con gradiente continuo de color en tau/beta.
  - El centroide r_c y el punto de inicio/cierre tau = 0.
  - Parámetros físicos: beta, lambda_th, radio de giro R_gyr y dispersion media.
"""

import argparse
from pathlib import Path
import re
import shutil
import subprocess

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np


OMEGA = 1.0
DEFAULT_BEADS = 64
DEFAULT_T_HIGH = 3.0
DEFAULT_T_LOW = 0.2
DEFAULT_INTERMEDIATE_TEMPERATURES = (1.0, 0.5)
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR
OUTPUT_DIR = PROJECT_DIR / 'Resultados' / '2DArmonico'
BINARY_MC = PROJECT_DIR / 'MC'
TEMPLATE_CFG = PROJECT_DIR / 'MC.cfg'


def load_coords_from_in3dprev(filepath):
    """Carga las coordenadas (x, y) de las cuentas desde el archivo in3Dprev.in de MC."""
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = [line.strip() for line in f if line.strip()]

    n_beads = int(lines[0])
    coords = []
    idx = 1
    for _ in range(n_beads):
        # Cada walker tiene: línea de energía y línea de coordenadas x y z
        idx += 1  # saltar línea de energía
        parts = [float(val) for val in lines[idx].split()]
        coords.append([parts[0], parts[1]])
        idx += 1

    return np.array(coords)


def run_c_simulation(temperature, n_beads, run_dir, binary_path=BINARY_MC, template_cfg_path=TEMPLATE_CFG,
                       blck_heating=25, blck=60, niter=100, omega=OMEGA):
    """Ejecuta una corrida aislada del ejecutable en C (MC) para obtener el polímero térmico."""
    run_dir.mkdir(parents=True, exist_ok=True)
    if not binary_path.is_file():
        raise FileNotFoundError(f'No se encontró el binario MC compilado en: {binary_path}')

    shutil.copy2(binary_path, run_dir / 'MC')
    cfg_text = template_cfg_path.read_text(encoding='utf-8', errors='ignore')

    def replace_line(text, key, new_val):
        pattern = rf'(\b{key}\s*=\s*)[^\s³]+'
        return re.sub(pattern, rf'\g<1>{new_val}', text)

    cfg_text = replace_line(cfg_text, 'MC', '3')
    cfg_text = replace_line(cfg_text, 'T', f'{temperature:g}')
    cfg_text = replace_line(cfg_text, 'Npop', str(n_beads))
    cfg_text = replace_line(cfg_text, 'NwalkersMax', str(n_beads + 20))
    cfg_text = replace_line(cfg_text, 'N', '1')
    cfg_text = replace_line(cfg_text, 'omega_x', f'{omega:g}')
    cfg_text = replace_line(cfg_text, 'omega_y', f'{omega:g}')
    cfg_text = replace_line(cfg_text, 'alpha_x', '1')
    cfg_text = replace_line(cfg_text, 'alpha_y', '1')
    cfg_text = replace_line(cfg_text, 'blck_heating', str(blck_heating))
    cfg_text = replace_line(cfg_text, 'blck', str(blck))
    cfg_text = replace_line(cfg_text, 'Niter', str(niter))
    cfg_text = replace_line(cfg_text, 'file_particles', 'in3Dprev.in')

    (run_dir / 'MC.cfg').write_text(cfg_text, encoding='utf-8')

    log_path = run_dir / 'run.log'
    print(f'-> Ejecutando simulación en C (MC) para T={temperature:g}, beads={n_beads}...')
    with log_path.open('w', encoding='utf-8') as log_file:
        res = subprocess.run(['./MC'], cwd=run_dir, stdout=log_file, stderr=subprocess.STDOUT, text=True)

    if res.returncode != 0:
        raise RuntimeError(f'Error al ejecutar ./MC en {run_dir}. Revisa el log: {log_path}')

    coords_file = run_dir / 'in3Dprev.in'
    if not coords_file.is_file():
        raise FileNotFoundError(f'No se generó el archivo de coordenadas {coords_file}')

    coords = load_coords_from_in3dprev(coords_file)
    if coords.shape != (n_beads, 2):
        raise RuntimeError(
            f'El camino guardado tiene forma {coords.shape}; se esperaban '
            f'{n_beads} beads con coordenadas x,y.'
        )
    if not np.isfinite(coords).all() or np.any(np.std(coords, axis=0) < 1e-10):
        raise RuntimeError(
            'El camino no presenta fluctuaciones válidas en x e y. '
            'Recompila MC con TRIAL_2D y el muestreador armónico 2D de pimc.c.'
        )
    beta = 1.0 / temperature
    return beta, coords


def compute_path_observables(path, beta):
    """Calcula centroide, radio de giro cuántico, distancia media y lambda térmica."""
    centroid = np.mean(path, axis=0)
    fluctuations = path - centroid
    r_gyr = np.sqrt(np.mean(np.sum(fluctuations**2, axis=1)))
    r_mean = np.sqrt(np.mean(np.sum(path**2, axis=1)))
    lambda_th = np.sqrt(2.0 * np.pi * beta)
    return centroid, r_gyr, r_mean, lambda_th


def plot_paths(samples, output_file, n_beads, omega=OMEGA):
    """Genera caminos representativos en 2D sobre el potencial armónico."""
    extent = max(np.abs(path).max() for _, _, path in samples) * 1.30
    extent = max(extent, 2.8)  # margen mínimo para apreciar el pozo

    n_columns = 2
    n_rows = (len(samples) + n_columns - 1) // n_columns
    figure, axes = plt.subplots(n_rows, n_columns, figsize=(12, 5.2 * n_rows),
                                sharex=True, sharey=True, squeeze=False)
    axes = axes.ravel()

    cmap = plt.get_cmap('plasma')

    grid_coords = np.linspace(-extent, extent, 250)
    xx, yy = np.meshgrid(grid_coords, grid_coords)
    vv = 0.5 * omega**2 * (xx**2 + yy**2)

    for axis, (temp, beta, path) in zip(axes, samples):
        # 1. Contornos del potencial V(x, y) = 1/2 omega^2 (x^2 + y^2)
        contours = axis.contour(xx, yy, vv, levels=7, colors='gray', alpha=0.35, linewidths=0.9)
        axis.clabel(contours, inline=True, fontsize=7, fmt='V=%.1f')

        # 2. Cerrar el camino cíclico
        closed_path = np.vstack((path, path[0]))
        segments = np.stack((closed_path[:-1], closed_path[1:]), axis=1)
        colors = np.linspace(0.0, 1.0, len(path))

        # 3. Dibujar el polímero continuo con gradiente en tiempo imaginario
        collection = LineCollection(segments, cmap=cmap, linewidth=2.0, alpha=0.9, zorder=3)
        collection.set_array(colors)
        axis.add_collection(collection)

        # 4. Beads individuales
        axis.scatter(path[:, 0], path[:, 1], c=colors, cmap=cmap, s=18, alpha=0.85, zorder=4)

        # 5. Punto inicial tau = 0
        axis.scatter(path[0, 0], path[0, 1], marker='o', s=60, color='tab:red',
                     edgecolor='white', linewidth=1.0, zorder=6, label=r'Inicio/Fin ($\tau=0, \beta$)')

        # 6. Centroide
        centroid, r_gyr, r_mean, lambda_th = compute_path_observables(path, beta)
        axis.scatter(centroid[0], centroid[1], marker='X', s=90, color='black',
                     edgecolor='yellow', linewidth=0.8, zorder=7, label=r'Centroide $\mathbf{r}_c$')

        # 7. Formato y leyendas
        if temp / omega >= 2.0:
            regimen = f'Alta T (T/ω = {temp / omega:g})'
            desc = 'Camino compacto alrededor del centroide'
        elif temp / omega <= 0.5:
            regimen = f'Baja T (T/ω = {temp / omega:g})'
            desc = 'Fluctuaciones cuánticas internas más amplias'
        else:
            regimen = f'Zona de cruce (T/ω = {temp / omega:g})'
            desc = 'Contribuyen efectos térmicos y cuánticos'

        axis.set_title(f'{regimen}\n{desc}', fontsize=10.5, pad=8)
        axis.set_xlabel('x', fontsize=11)
        axis.set_xlim(-extent, extent)
        axis.set_ylim(-extent, extent)
        axis.set_aspect('equal', adjustable='box')
        axis.grid(True, linestyle=':', alpha=0.3)
        axis.legend(loc='upper right', fontsize=8.5, framealpha=0.9)

        info_text = (
            f'$\\beta = 1/T = {beta:.2f}$\n'
            f'$\\lambda_{{\\mathrm{{th}}}} = {lambda_th:.2f}$\n'
            f'$R_{{\\mathrm{{gyr}}}}$ (camino) $= {r_gyr:.2f}$\n'
            f'$|\\mathbf{{r}}_c| = {np.linalg.norm(centroid):.2f}$\n'
            f'$\\sqrt{{\\langle r^2 \\rangle}}$ (camino) $= {r_mean:.2f}$'
        )
        axis.text(0.04, 0.05, info_text, transform=axis.transAxes, fontsize=8.5,
                  verticalalignment='bottom', bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.85, edgecolor='gray'))

    for axis in axes[len(samples):]:
        axis.set_visible(False)
    for axis in axes[:len(samples)]:
        axis.set_ylabel('y', fontsize=11)
    figure.suptitle('Caminos representativos en tiempo imaginario: oscilador armónico 2D',
                    fontsize=13, fontweight='bold')

    figure.text(0.5, 0.015,
                'Una trayectoria por temperatura. El centroide también fluctúa; '
                'las magnitudes indicadas no son promedios de conjunto.',
                ha='center', fontsize=9)
    figure.tight_layout(rect=(0.02, 0.05, 0.88, 0.95))
    colorbar = figure.colorbar(collection, ax=axes[:len(samples)], fraction=0.035, pad=0.03)
    colorbar.set_label(r'Tiempo imaginario relativo $\tau / \beta$', fontsize=10)

    figure.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(figure)
    print(f'Figura guardada exitosamente en:\n  {output_file}')


def plot_centroid_fluctuations(samples, output_file):
    """Compara la distancia interna al centroide a lo largo del camino."""
    figure, axis = plt.subplots(figsize=(8, 4.8))
    colors = plt.get_cmap('viridis')(np.linspace(0.08, 0.92, len(samples)))

    for color, (temperature, beta, path) in zip(colors, samples):
        centroid = np.mean(path, axis=0)
        radius = np.linalg.norm(path - centroid, axis=1)
        tau_fraction = np.arange(len(path)) / len(path)
        r_gyr = np.sqrt(np.mean(radius**2))
        axis.plot(tau_fraction, radius, color=color, linewidth=1.2,
                  label=f'T={temperature:g}; $R_{{gyr}}$={r_gyr:.2f}')
        axis.axhline(r_gyr, color=color, linestyle=':', linewidth=1.0, alpha=0.8)

    axis.set_title('Fluctuación interna a lo largo del tiempo imaginario')
    axis.set_xlabel(r'Tiempo imaginario relativo $\tau / \beta$')
    axis.set_ylabel(r'$|\mathbf{r}(\tau)-\mathbf{r}_c|$')
    axis.set_xlim(0.0, 1.0)
    axis.grid(True, alpha=0.3)
    axis.legend(title='Una trayectoria por T', ncol=2)
    figure.tight_layout()
    figure.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(figure)
    print(f'Figura de fluctuaciones guardada en:\n  {output_file}')


def main():
    parser = argparse.ArgumentParser(
        description='Generar y graficar caminos en tiempo imaginario en 2D a partir de la simulación PIMC en C.'
    )
    parser.add_argument('--t-high', type=float, default=DEFAULT_T_HIGH,
                        help=f'Temperatura alta (por defecto: {DEFAULT_T_HIGH})')
    parser.add_argument('--t-low', type=float, default=DEFAULT_T_LOW,
                        help=f'Temperatura baja (por defecto: {DEFAULT_T_LOW})')
    parser.add_argument('--intermediate-temperatures', type=float, nargs='*',
                        default=DEFAULT_INTERMEDIATE_TEMPERATURES,
                        help='Temperaturas intermedias (por defecto: 1.0 0.5)')
    parser.add_argument('--beads', type=int, default=DEFAULT_BEADS,
                        help=f'Número de cuentas/beads (por defecto: {DEFAULT_BEADS})')
    parser.add_argument('--heating', type=int, default=25,
                        help='Bloques de calentamiento en MC (por defecto: 25)')
    parser.add_argument('--blocks', type=int, default=60,
                        help='Bloques de medida en MC (por defecto: 60)')
    parser.add_argument('--iterations', type=int, default=100,
                        help='Iteraciones por bloque en MC (por defecto: 100)')
    parser.add_argument('--use-existing', action='store_true',
                        help='Usa coordenadas existentes en Resultados/2DArmonico/ si ya están calculadas')
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    temperatures = sorted(
        set((args.t_high, *args.intermediate_temperatures, args.t_low)), reverse=True
    )
    if any(temperature <= 0.0 for temperature in temperatures):
        parser.error('Todas las temperaturas deben ser positivas.')
    samples = []

    for temp in temperatures:
        sim_dir = OUTPUT_DIR / f'sim_T{temp:g}'
        csv_path = OUTPUT_DIR / f'camino_T{temp:g}.csv'
        beta = 1.0 / temp

        if args.use_existing and csv_path.is_file():
            print(f'-> Cargando camino existente para T={temp:g} desde {csv_path}')
            data = np.loadtxt(csv_path, delimiter=',', skiprows=1)
            path = data[:, 1:3]
        else:
            beta, path = run_c_simulation(
                temperature=temp,
                n_beads=args.beads,
                run_dir=sim_dir,
                blck_heating=args.heating,
                blck=args.blocks,
                niter=args.iterations,
            )
            # Guardar coordenadas limpias en CSV
            tau_grid = np.linspace(0.0, beta, len(path), endpoint=False)
            np.savetxt(
                csv_path,
                np.column_stack((tau_grid, path[:, 0], path[:, 1])),
                delimiter=',',
                header='tau,x,y',
                comments='',
            )
            print(f'   Coordenadas (T={temp:g}) guardadas en: {csv_path}')

        samples.append((temp, beta, path))

    output_file = OUTPUT_DIR / 'caminos_armonico_2D.png'
    plot_paths(samples, output_file, args.beads)
    plot_centroid_fluctuations(
        samples, OUTPUT_DIR / 'fluctuaciones_centroide_2D.png'
    )


if __name__ == '__main__':
    main()
