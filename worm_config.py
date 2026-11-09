"""Shared Worm simulation configuration helpers."""

import os
from pathlib import Path
import re


def get_n_particles(config_path, legacy_env_var=None):
    value = os.environ.get("WORM_PARTICLES")
    if value is None and legacy_env_var is not None:
        value = os.environ.get(legacy_env_var)

    if value is not None:
        particle_count = int(value)
    else:
        config_text = Path(config_path).read_text(encoding="utf-8", errors="ignore")
        match = re.search(r"(?<![A-Za-z0-9_])N\s*=\s*(\d+)", config_text)
        if match is None:
            raise ValueError(f"No se encontró N= en {config_path}")
        particle_count = int(match.group(1))

    if particle_count < 1:
        raise ValueError("El número de partículas N debe ser >= 1")
    return particle_count