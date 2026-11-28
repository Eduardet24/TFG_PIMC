## @file Script_Armonico1D_basico.py
#  @brief Analiza la densidad de PIMC del oscilador armónico 1D.
#
#  Este script compara la distribución de posiciones obtenida mediante PIMC
#  con la solución analítica del oscilador armónico cuántico a temperatura
#  finita. También calcula el radio cuadrático medio y estima las energías
#  cinética, potencial y total usando el teorema del virial.
#
#  @details
#  El programa MC guarda la distribución radial en `outrd.dat`. En 1D, esta
#  salida representa la rama x >= 0 de una distribución simétrica y contiene
#  varios bloques de Monte Carlo consecutivos. El script separa los bloques,
#  descarta el periodo de termalización, comprueba que tienen la misma malla
#  y calcula su promedio para reducir el ruido estadístico.
#
#  Para una partícula en un oscilador armónico con hbar = m = omega = 1, la
#  predicción analítica usada es:
#  @f[
#  n(x) = \left(\frac{\tanh(\beta/2)}{\pi}\right)^{D/2}
#         \exp\left[-x^2\tanh(\beta/2)\right].
#  @f]
#
#  El valor exacto del radio cuadrático medio y de la energía total es:
#  @f[
#  \langle x^2\rangle = E = \frac{D}{2\tanh(\beta/2)}.
#  @f]
#
#  En este problema se cumple Ekin = Epot = <x^2>/2. La salida `oute.dat`
#  contiene ahora la energía total estimada directamente por PIMC; sus
#  contribuciones cinética y potencial se siguen comparando por separado
#  mediante la densidad y el teorema del virial.
#
#  @pre `MC.cfg`, `outrd.dat` y `oute.dat` deben estar en el directorio de ejecución.
#  @post Se guarda el resumen numérico y las figuras de densidad y residuos.

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from pimc_common import read_cfg_value

PROJECT_DIR = _SCRIPTS_DIR.parent
DEFAULT_RESULTS_DIR = PROJECT_DIR / 'Resultados' / '1DArmonico' / 'ejecucion_actual'
D = 1.0


def load_density_blocks(rd_path, burn_in=0.5):
    """Carga los bloques de `outrd.dat`, valida la malla y aplica el descarte de termalización."""
    data = np.loadtxt(rd_path, comments='#', ndmin=2)
    if data.size == 0 or data.shape[1] < 2:
        raise ValueError(f'{rd_path} no contiene columnas válidas de densidad')

    saltos = np.where(np.diff(data[:, 0]) < 0)[0] + 1
    bloques = np.split(data, saltos)

    if any(len(bloque) != len(bloques[0]) for bloque in bloques):
        raise ValueError('Los bloques de outrd.dat no tienen el mismo tamaño')

    x_sim = bloques[0][:, 0]
    if x_sim.size < 2:
        raise ValueError('Se necesitan al menos 2 bins en outrd.dat')
    if not all(np.allclose(bloque[:, 0], x_sim) for bloque in bloques[1:]):
        raise ValueError('Los centros de bin cambian entre bloques')

    n_bloques = np.array([bloque[:, 1] for bloque in bloques], dtype=float)
    inicio_valido = min(int(len(n_bloques) * burn_in), len(n_bloques) - 1)
    n_bloques_utiles = n_bloques[inicio_valido:]

    n_sim = np.mean(n_bloques_utiles, axis=0)
    if len(n_bloques_utiles) > 1:
        n_error = np.std(n_bloques_utiles, axis=0, ddof=1) / np.sqrt(len(n_bloques_utiles))
    else:
        n_error = np.zeros_like(n_sim)

    return x_sim, n_sim, n_error, len(n_bloques_utiles), len(bloques)


def load_energy_stats(energy_path, burn_in=0.5):
    """Carga `oute.dat`, descarta la termalización y calcula medias y desviaciones."""
    energy_data = np.loadtxt(energy_path, comments='#', ndmin=2)
    if energy_data.size == 0:
        raise ValueError(f'{energy_path} no contiene mediciones de energía')

    equilibrio_inicio = min(int(energy_data.shape[0] * burn_in), energy_data.shape[0] - 1)
    measurements = energy_data[equilibrio_inicio:]
    ddof = 1 if measurements.shape[0] > 1 else 0

    stats = {
        'etot': (float(np.mean(measurements[:, 0])), float(np.std(measurements[:, 0], ddof=ddof))),
    }
    if measurements.shape[1] >= 3:
        stats['epot'] = (float(np.mean(measurements[:, 1])), float(np.std(measurements[:, 1], ddof=ddof)))
        stats['ekin_vir'] = (float(np.mean(measurements[:, 2])), float(np.std(measurements[:, 2], ddof=ddof)))
    if measurements.shape[1] >= 4:
        stats['ekin_prim'] = (float(np.mean(measurements[:, 3])), float(np.std(measurements[:, 3], ddof=ddof)))
    return stats


def save_single_point_convergence(output_dir, n_beads, r2_sim, r2_teo, e_pimc, e_teo, save_plot=False):
    """Exporta `convergencia_beads_1D.csv` y, si se solicita, la gráfica de un único punto."""
    inverso_n = 1.0 / n_beads
    tabla_estudio = np.array([[n_beads, inverso_n, r2_sim, r2_teo, e_pimc, e_teo]], dtype=float)
    np.savetxt(
        output_dir / 'convergencia_beads_1D.csv',
        tabla_estudio,
        delimiter=',',
        header='N_beads,1/N,x2_PIMC,x2_analitico,E_PIMC,E_analitica',
        comments='',
        fmt=['%d', '%.8g', '%.10g', '%.10g', '%.10g', '%.10g'],
    )
    if not save_plot:
        return

    figura, ejes = plt.subplots(1, 2, figsize=(11, 4.5))
    ejes[0].plot([inverso_n], [r2_sim], 'o', label=rf'PIMC $\langle x^2 \rangle$ (N={n_beads})')
    ejes[0].axhline(r2_teo, color='black', linestyle='--', label='Analítico')
    ejes[0].set_xlabel(r'$1/N_{beads}$')
    ejes[0].set_ylabel(r'$\langle x^2 \rangle$')
    ejes[0].set_title(r'Convergencia de $\langle x^2 \rangle$')
    ejes[0].legend()
    ejes[0].grid(True)

    ejes[1].plot([inverso_n], [e_pimc], 's', color='tab:orange', label=f'PIMC (N={n_beads})')
    ejes[1].axhline(e_teo, color='black', linestyle='--', label='Analítica')
    ejes[1].set_xlabel(r'$1/N_{beads}$')
    ejes[1].set_ylabel('Energía total')
    ejes[1].set_title('Convergencia de la energía')
    ejes[1].legend()
    ejes[1].grid(True)

    figura.tight_layout()
    figura.savefig(output_dir / 'convergencia_beads_1D.png', dpi=300, bbox_inches='tight')
    plt.close(figura)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description='Analiza la corrida actual de PIMC en oscilador armónico 1D.')
    parser.add_argument('-T', '--temperature', type=float, default=None,
                        help='Temperatura T (por defecto se lee de MC.cfg)')
    parser.add_argument('--burn-in', type=float, default=0.5,
                        help='Fracción inicial de bloques/muestras descartadas por termalización (por defecto: 0.5)')
    parser.add_argument('--input-dir', type=Path, default=PROJECT_DIR,
                        help='Directorio con MC.cfg, outrd.dat y oute.dat')
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_RESULTS_DIR,
                        help='Directorio donde guardar resultados y figuras')
    parser.add_argument('--plot-convergencia', action='store_true',
                        help='Genera también la gráfica de un único punto convergencia_beads_1D.png')
    parser.add_argument('--no-show', action='store_true',
                        help='No abre ventanas interactivas de matplotlib al finalizar')
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if not (0.0 <= args.burn_in < 1.0):
        raise ValueError('--burn-in debe estar en el intervalo [0.0, 1.0)')

    input_dir = args.input_dir
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    config_text = (input_dir / 'MC.cfg').read_text(encoding='utf-8', errors='ignore')
    temperature = args.temperature if args.temperature is not None else read_cfg_value(config_text, 'T')
    if temperature <= 0.0:
        raise ValueError(f'La temperatura T debe ser positiva; recibido {temperature}')
    n_beads_actual = int(read_cfg_value(config_text, 'Npop'))
    beta = 1.0 / temperature

    x_sim, n_sim, n_error, used_blocks, total_blocks = load_density_blocks(
        input_dir / 'outrd.dat', burn_in=args.burn_in
    )

    x_teo = np.linspace(x_sim.min(), x_sim.max(), 500)
    alpha = np.tanh(beta / 2.0)
    prefactor = (alpha / np.pi) ** (D / 2.0)
    n_teo = prefactor * np.exp(-(x_teo ** 2) * alpha)

    r2_teo = D / (2.0 * alpha)
    e_kin_teo = r2_teo / 2.0
    e_pot_teo = r2_teo / 2.0
    e_teo = e_kin_teo + e_pot_teo

    dx = x_sim[1] - x_sim[0]
    r2_sim = 2.0 * float(np.sum(x_sim ** 2 * n_sim)) * dx
    e_pot_sim = r2_sim / 2.0
    e_kin_sim = r2_sim / 2.0
    e_sim = e_kin_sim + e_pot_sim

    energy_stats = load_energy_stats(input_dir / 'oute.dat', burn_in=args.burn_in)
    e_pimc, e_pimc_std = energy_stats['etot']

    n_teo_sim = prefactor * np.exp(-(x_sim ** 2) * alpha)
    residuo = n_sim - n_teo_sim
    residuo_rms = float(np.sqrt(np.mean(residuo ** 2)))

    resultados_txt = output_dir / 'resultados_graficas_armonico_1D.txt'
    with resultados_txt.open('w', encoding='utf-8') as resultados:
        print(f'T = {temperature:g}, N_beads = {n_beads_actual}, bloques utiles = {used_blocks}/{total_blocks}', file=resultados)
        print(f'Ekin simulacion = {e_kin_sim:.6f}', file=resultados)
        print(f'Epot simulacion = {e_pot_sim:.6f}', file=resultados)
        print(f'E total simulacion = {e_sim:.6f}', file=resultados)
        print(f'E total PIMC directo = {e_pimc:.6f} +/- {e_pimc_std:.6f}', file=resultados)
        if 'epot' in energy_stats:
            e_pot_pimc, e_pot_pimc_std = energy_stats['epot']
            e_kin_pimc, e_kin_pimc_std = energy_stats['ekin_vir']
            print(f'Epot PIMC directo = {e_pot_pimc:.6f} +/- {e_pot_pimc_std:.6f}', file=resultados)
            print(f'Ekin (virial) PIMC = {e_kin_pimc:.6f} +/- {e_kin_pimc_std:.6f}', file=resultados)
        if 'ekin_prim' in energy_stats:
            e_kin_prim_pimc, e_kin_prim_pimc_std = energy_stats['ekin_prim']
            print(f'Ekin (primitivo) PIMC = {e_kin_prim_pimc:.6f} +/- {e_kin_prim_pimc_std:.6f}', file=resultados)
        print(f'Ekin teoria = {e_kin_teo:.6f}', file=resultados)
        print(f'Epot teoria = {e_pot_teo:.6f}', file=resultados)
        print(f'E total teoria = {e_teo:.6f}', file=resultados)
        print(f'<x^2> simulacion = {r2_sim:.6f}', file=resultados)
        print(f'<x^2> teoria = {r2_teo:.6f}', file=resultados)
        print(f'Residuo RMS = {residuo_rms:.6f}', file=resultados)

    fig_dens, ax_dens = plt.subplots()
    ax_dens.plot(x_teo, n_teo, 'r-', linewidth=2, label=f'Teoría (T={temperature:g})')
    ax_dens.plot(x_sim, n_sim, 'b.', markersize=5,
                 label=f'Simulación PIMC ({used_blocks}/{total_blocks} bloques)')
    ax_dens.set_title('Densidad Cartesiana del Oscilador Armónico 1D')
    ax_dens.set_xlabel('Posición (z)')
    ax_dens.set_ylabel('Densidad n(z)')
    ax_dens.legend()
    ax_dens.grid(True)
    fig_dens.savefig(output_dir / 'densidad_armonico_1D.png', dpi=300, bbox_inches='tight')

    fig_res, ax_res = plt.subplots()
    ax_res.axhline(0.0, color='black', linewidth=1)
    ax_res.errorbar(x_sim, residuo, yerr=n_error, fmt='g.', markersize=5,
                    capsize=2, label='Residuo +/- error estándar')
    ax_res.set_title('Residuo de la densidad PIMC respecto a la teoría')
    ax_res.set_xlabel('Posición (z)')
    ax_res.set_ylabel(r'$n_{PIMC}(z) - n_{teoria}(z)$')
    ax_res.legend()
    ax_res.grid(True)
    fig_res.savefig(output_dir / 'residuo_armonico_1D.png', dpi=300, bbox_inches='tight')

    save_single_point_convergence(
        output_dir,
        n_beads=n_beads_actual,
        r2_sim=r2_sim,
        r2_teo=r2_teo,
        e_pimc=e_pimc,
        e_teo=e_teo,
        save_plot=args.plot_convergencia,
    )

    if args.no_show:
        plt.close('all')
    else:
        plt.show()


if __name__ == '__main__':
    main()