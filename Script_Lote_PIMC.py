"""Run a serial matrix of isolated PIMC simulations."""

import argparse
import re
import shutil
import subprocess
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = PROJECT_DIR / 'Resultados' / '1DArmonico'
DEFAULT_BEAD_COUNTS = (5, 10, 15, 20, 25, 30)
DEFAULT_RUN_MATRIX = (
    (0.5, DEFAULT_BEAD_COUNTS),
    (1.0, DEFAULT_BEAD_COUNTS),
    (2.0, (5, 10, 15, 20)),
    (3.0, (5,)),
    (4.0, (5,)),
)
INPUT_FILE_KEYS = ('file_particles', 'file_wf', 'file_wf1b')


def parse_values(text, cast):
    try:
        values = tuple(cast(value.strip()) for value in text.split(','))
    except ValueError as error:
        raise argparse.ArgumentTypeError(f'Lista de valores no válida: {text}') from error
    if not values:
        raise argparse.ArgumentTypeError('La lista de valores no puede estar vacía')
    return values


def replace_setting(config_text, key, value):
    number = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
    pattern = rf'(\b{re.escape(key)}\s*=\s*){number}'
    updated, replacements = re.subn(
        pattern, lambda match: f'{match.group(1)}{value}', config_text, count=1
    )
    if replacements != 1:
        raise ValueError(f'No se pudo sustituir {key} en MC.cfg')
    return updated


def read_config_value(config_text, key):
    number = r'([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)'
    match = re.search(rf'\b{re.escape(key)}\s*=\s*{number}', config_text)
    if match is None:
        raise ValueError(f'No se pudo leer {key} de MC.cfg')
    return float(match.group(1))


def configure_run(template, temperature, beads):
    config = replace_setting(template, 'T', f'{temperature:g}')
    return replace_setting(config, 'Npop', str(beads))


def expected_measurements(config_text):
    blocks = int(read_config_value(config_text, 'blck'))
    iterations = int(read_config_value(config_text, 'Niter'))
    interval = int(read_config_value(config_text, 'Nmeasure'))
    if interval <= 0 or iterations % interval:
        raise ValueError('Niter debe ser divisible por Nmeasure')
    return blocks * (iterations // interval)


def valid_output(run_dir, required_rows):
    energy_file = run_dir / 'oute.dat'
    if not energy_file.is_file():
        return False
    with energy_file.open('r', encoding='utf-8') as energy_data:
        header = energy_data.readline().strip()
        rows = sum(1 for line in energy_data if line.strip())
    return (header in ('# Etot Epot Ekin', '# Etot Epot Ekin Ekin_prim')) and rows >= required_rows


def completed_run(run_dir, required_rows):
    return (run_dir / '.complete').is_file() and valid_output(run_dir, required_rows)


def copy_run_inputs(run_dir, config_text):
    (run_dir / 'MC.cfg').write_text(config_text, encoding='utf-8')
    shutil.copy2(PROJECT_DIR / 'MC', run_dir / 'MC')

    for key in INPUT_FILE_KEYS:
        match = re.search(rf'\b{key}\s*=\s*(\S+)', config_text)
        if match is None:
            continue
        source = PROJECT_DIR / match.group(1)
        if source.is_file():
            shutil.copy2(source, run_dir / source.name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', action='store_true',
                        help='ejecuta las simulaciones; sin esta opción solo muestra el plan')
    parser.add_argument('--temperatures',
                        help='temperaturas separadas por comas; requiere --beads')
    parser.add_argument('--beads',
                        help='números de beads separados por comas; requiere --temperatures')
    parser.add_argument('--campaign', default='lote_modos_normales_blck1000',
                        help='nombre del subdirectorio para guardar esta campaña')
    args = parser.parse_args()

    if (args.temperatures is None) != (args.beads is None):
        parser.error('--temperatures y --beads deben usarse juntos')
    if args.temperatures is None:
        run_matrix = DEFAULT_RUN_MATRIX
    else:
        temperatures = parse_values(args.temperatures, float)
        bead_counts = parse_values(args.beads, int)
        run_matrix = tuple((temperature, bead_counts) for temperature in temperatures)

    if any(temperature <= 0 for temperature, _ in run_matrix):
        parser.error('todas las temperaturas deben ser positivas')
    if any(beads <= 0 for _, counts in run_matrix for beads in counts):
        parser.error('todos los números de beads deben ser positivos')
    if not re.fullmatch(r'[A-Za-z0-9_-]+', args.campaign):
        parser.error('--campaign solo puede contener letras, números, guion y guion bajo')

    binary = PROJECT_DIR / 'MC'
    if not binary.is_file():
        raise FileNotFoundError(f'No se encontró el ejecutable {binary}; compila MC primero')
    template = (PROJECT_DIR / 'MC.cfg').read_text(encoding='utf-8')
    acceptance_rate = read_config_value(template, 'acceptance_rate')
    root = RESULTS_DIR / f'AR{acceptance_rate:g}' / args.campaign
    planned_runs = []

    for temperature, bead_counts in run_matrix:
        for beads in bead_counts:
            run_dir = root / f'Temp{temperature:g}' / f'N{beads}'
            run_config = configure_run(template, temperature, beads)
            expected_rows = expected_measurements(run_config)
            planned_runs.append((temperature, beads, run_dir, run_config, expected_rows))

    print(f'Matriz: {len(run_matrix)} temperaturas, {len(planned_runs)} corridas')
    print(f'Aceptación: AR{acceptance_rate:g}')
    print(f'Calentamiento por corrida: {int(read_config_value(template, "blck_heating"))} bloques')
    print(f'Mediciones esperadas por corrida: {expected_measurements(template):,}')
    print(f'Destino: {root}')
    if not args.run:
        print('Modo de prueba: no se ejecutó MC. Usa --run para lanzar el lote.')

    for temperature, beads, run_dir, run_config, expected_rows in planned_runs:
        if completed_run(run_dir, expected_rows):
            print(f'OMITIR completo: T={temperature:g}, N={beads}')
            continue
        if run_dir.exists() and any(run_dir.iterdir()):
            print(f'BLOQUEAR carpeta no vacía (no se sobrescribe): {run_dir}')
            continue

        print(f'PREPARAR T={temperature:g}, N={beads}: {run_dir}')
        if not args.run:
            continue

        run_dir.mkdir(parents=True, exist_ok=True)
        copy_run_inputs(run_dir, run_config)
        with (run_dir / 'run.log').open('w', encoding='utf-8') as log_file:
            result = subprocess.run(
                ['./MC'], cwd=run_dir, stdout=log_file, stderr=subprocess.STDOUT,
                check=False,
            )

        if result.returncode != 0:
            print(f'FALLÓ T={temperature:g}, N={beads}; revisa {run_dir / "run.log"}')
            continue
        if not valid_output(run_dir, expected_rows):
            print(f'INCOMPLETA T={temperature:g}, N={beads}; revisa {run_dir / "run.log"}')
            continue

        (run_dir / '.complete').write_text('ok\n', encoding='ascii')
        print(f'COMPLETA T={temperature:g}, N={beads}')


if __name__ == '__main__':
    main()