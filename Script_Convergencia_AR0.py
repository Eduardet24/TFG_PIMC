"""Plot AR0 and AR50 observables with quadratic bead-count extrapolations."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


RESULTS_DIR = Path(__file__).resolve().parent / 'Resultados' / '1DArmonico'


def load_points(series):
    series_dir = RESULTS_DIR / series
    points = []
    for run_dir in series_dir.glob('N*'):
        csv_file = run_dir / 'convergencia_beads_1D.csv'
        if not csv_file.is_file():
            continue

        rows = np.atleast_2d(np.loadtxt(csv_file, delimiter=',', skiprows=1))
        match = run_dir.name.removeprefix('N')
        if not match.isdigit():
            continue
        n_beads = int(match)
        matching_rows = rows[rows[:, 0].astype(int) == n_beads]
        if matching_rows.shape[0] != 1:
            raise ValueError(f'{csv_file} debe contener exactamente una fila para N={n_beads}')

        row = matching_rows[0]
        points.append((n_beads, row[2], row[3], row[4], row[5]))

    points.sort(key=lambda point: point[0], reverse=True)
    if len(points) < 3:
        raise ValueError(f'Se necesitan al menos 3 valores {series} para el ajuste parabólico; hay {len(points)}')
    return np.asarray(points, dtype=float)


def plot_series(series):
    points = load_points(series)
    n_beads = points[:, 0]
    inverse_n = 1.0 / n_beads
    x2_values, x2_theory = points[:, 1], points[:, 2]
    energy_values, energy_theory = points[:, 3], points[:, 4]

    x2_fit = np.polyfit(inverse_n, x2_values, 2)
    energy_fit = np.polyfit(inverse_n, energy_values, 2)
    x_curve = np.linspace(0.0, inverse_n.max() * 1.05, 400)

    figure, axes = plt.subplots(1, 2, figsize=(12, 5))
    acceptance_label = f'AR={series.removeprefix("AR")}'
    plot_specs = (
        (axes[0], x2_values, x2_theory[0], x2_fit, r'$\langle x^2 \rangle$', '⟨x²⟩'),
        (axes[1], energy_values, energy_theory[0], energy_fit, 'Energía total', 'Energía'),
    )

    for axis, values, theory, fit, ylabel, label in plot_specs:
        extrapolated = np.polyval(fit, 0.0)
        axis.scatter(inverse_n, values, s=65, zorder=3, label=f'Simulación {acceptance_label}')
        for beads, x_value, y_value in zip(n_beads, inverse_n, values):
            axis.annotate(f'N={int(beads)}', (x_value, y_value), xytext=(5, 7),
                          textcoords='offset points', fontsize=8)
        axis.plot(x_curve, np.polyval(fit, x_curve), color='tab:orange',
                  label='Ajuste cuadrático')
        axis.axhline(theory, color='black', linestyle='--', label='Analítico')
        axis.scatter([0.0], [extrapolated], marker='X', s=90, color='tab:red', zorder=4,
                     label=f'Extrapolación N→∞: {extrapolated:.3f}')
        axis.set_xlabel(r'$1/N_{beads}$  (el límite continuo es 0)')
        axis.set_ylabel(ylabel)
        axis.set_title(f'Convergencia de {label}')
        axis.grid(True, alpha=0.3)
        axis.legend()

        print(f'{label}: extrapolación N→∞ = {extrapolated:.8f}; analítico = {theory:.8f}')

    figure.suptitle(f'Convergencia PIMC con acceptance_rate = {series.removeprefix("AR")}')
    figure.tight_layout()
    output_file = RESULTS_DIR / series / f'convergencia_{series}.png'
    figure.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(figure)
    print(f'Gráfica guardada en {output_file}')


def main():
    for series in ('AR0', 'AR50'):
        plot_series(series)


if __name__ == '__main__':
    main()