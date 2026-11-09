#!/usr/bin/env python3
"""
Generador de esquemas didácticos y gráficos de transiciones aceptadas del Worm.

El modo --from-simulation lee snapshots before/after registrados por el código C.

Cada ejecución reemplaza automáticamente todos los archivos PNG (300 DPI) y PDF en:
  - Resultados/Figuras_Worm_Dornheim/
  - ../overleaf/imágenes/
"""

import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# Rutas de salida
SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR / "Resultados" / "Figuras_Worm_Dornheim"
OVERLEAF_DIR = SCRIPT_DIR.parent / "overleaf" / "imágenes"


# ==============================================================================
# CONFIGURACIÓN Y FUNCIONES COMUNES
# ==============================================================================

TIME = np.arange(21, dtype=float)
PATH_A = "#c84c4c"
PATH_B = "#4d8b78"
IRA = "#2458c6"
MASHA = "#bd6a20"
MUTED = "#78808d"


def save_figure(fig, filename):
    """Guarda y reemplaza automáticamente la figura en formato PNG y PDF en ambas carpetas."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    OVERLEAF_DIR.mkdir(parents=True, exist_ok=True)

    png_out = OUTPUT_DIR / filename
    pdf_out = png_out.with_suffix(".pdf")
    png_overleaf = OVERLEAF_DIR / filename
    pdf_overleaf = png_overleaf.with_suffix(".pdf")

    # Guardar en carpeta local de Resultados
    fig.savefig(png_out, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_out, bbox_inches="tight")

    # Reemplazar directamente en la carpeta de imágenes de Overleaf
    fig.savefig(png_overleaf, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_overleaf, bbox_inches="tight")

    plt.close(fig)
    print(f"  [REEMPLAZADO] {filename:<40} -> Resultados/ y overleaf/imágenes/")


# ==============================================================================
# ==============================================================================
# FORMATO 2: ESTILO DIDÁCTICO SUAVE (CURVAS CONTINUAS CON IRA Y MASHA)
# ==============================================================================

def setup_axis(ax, title):
    """Configuración didáctica de ejes con fondo cuadriculado suave."""
    ax.set_xlim(-2.2, 2.2)
    ax.set_ylim(0, 20)
    ax.set_xticks([-2, -1, 0, 1, 2])
    ax.set_yticks([0, 5, 10, 15, 20])
    ax.set_xlabel("x")
    ax.set_ylabel(r"$\tau/\epsilon$")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.grid(True, color="#d9dde3", linewidth=0.7, alpha=0.7)
    ax.set_axisbelow(True)
    ax.tick_params(direction="in", top=True, right=True)


def draw_track(ax, x, start=0, stop=None, color=PATH_A, linewidth=1.5,
               marker_size=4.2, zorder=2):
    if stop is None:
        stop = len(x) - 1
    tau = TIME[start:stop + 1]
    x_values = x[start:stop + 1]
    ax.plot(x_values, tau, color=color, linewidth=linewidth, zorder=zorder)
    ax.scatter(x_values, tau, s=marker_size**2, color=color, zorder=zorder + 1)


def mark_endpoint(ax, x, tau, label, color, offset=(7, 5), hollow=False):
    face = "white" if hollow else color
    ax.scatter([x], [tau], s=68, facecolor=face, edgecolor=color, linewidth=1.7,
               zorder=5)
    ax.annotate(label, (x, tau), xytext=offset, textcoords="offset points",
                color=color, fontsize=9, fontweight="bold", zorder=6)


def draw_open_path(ax, x, head, masha):
    draw_track(ax, x, stop=head, color=PATH_A)
    draw_track(ax, x, start=masha, color=PATH_A)
    mark_endpoint(ax, x[head], head, "Ira", IRA, offset=(7, 7))
    mark_endpoint(ax, x[masha], masha, "Masha", MASHA, offset=(7, -13))


def make_open_close():
    x_a = -1.05 + 0.28 * np.sin(2 * np.pi * TIME / 20 + 0.2)
    x_b = 1.05 + 0.27 * np.sin(2 * np.pi * TIME / 20 + 1.5)
    cut = 10

    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5), sharey=True)
    setup_axis(axes[0], "Z: caminos cerrados")
    draw_track(axes[0], x_a, color=PATH_A)
    draw_track(axes[0], x_b, color=PATH_B)

    setup_axis(axes[1], "G: enlace abierto")
    draw_track(axes[1], x_a, stop=cut, color=PATH_A)
    draw_track(axes[1], x_a, start=cut + 1, color=PATH_A)
    draw_track(axes[1], x_b, color=PATH_B)
    mark_endpoint(axes[1], x_a[cut], cut, "Ira", IRA, offset=(9, 12))
    mark_endpoint(axes[1], x_a[cut + 1], cut + 1, "Masha", MASHA, offset=(-58, -14))

    fig.suptitle("Open / Close: cortar o reconectar un enlace", fontweight="bold")
    fig.text(0.5, 0.015, "Open: Z → G     |     Close: G → Z (extremos adyacentes; τ periódico)",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.06, 1, 0.88))
    save_figure(fig, "worm_open_close.png")


def make_advance_recede():
    x = 0.1 + 0.35 * np.sin(2 * np.pi * TIME / 20 + 0.5)
    old_head = 7
    masha = 15
    advanced_head = 10
    receded_head = 4

    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5), sharey=True)
    setup_axis(axes[0], "Advance: Ira avanza en τ")
    draw_open_path(axes[0], x, advanced_head, masha)
    draw_track(axes[0], x, start=old_head, stop=advanced_head, color=IRA,
               linewidth=2.2, marker_size=5.0, zorder=3)
    mark_endpoint(axes[0], x[old_head], old_head, "Ira anterior", MUTED,
                  offset=(-64, -15), hollow=True)
    axes[0].annotate("beads añadidos", xy=(x[9], 9), xytext=(1.15, 8.5),
                     arrowprops={"arrowstyle": "->", "color": IRA}, color=IRA,
                     fontsize=8, ha="center")

    setup_axis(axes[1], "Recede: Ira retrocede en τ")
    draw_open_path(axes[1], x, receded_head, masha)
    mark_endpoint(axes[1], x[old_head], old_head, "Ira anterior", MUTED,
                  offset=(7, 7), hollow=True)
    axes[1].annotate("el gap se amplía", xy=(x[6], 6), xytext=(1.15, 8.5),
                     arrowprops={"arrowstyle": "->", "color": IRA}, color=IRA,
                     fontsize=8, ha="center")

    fig.suptitle("Advance / Recede: desplazar Ira y cambiar el gap", fontweight="bold")
    fig.text(0.5, 0.015, "En esta implementación, Advance/Recede se selecciona para N = 1.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.06, 1, 0.88))
    save_figure(fig, "worm_advance_recede.png")


def make_ira_spatial_move():
    x_before = 0.15 + 0.3 * np.sin(2 * np.pi * TIME / 20 + 0.3)
    head = 10
    masha = 15
    x_after = x_before.copy()
    x_after[head] += 0.72

    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5), sharey=True)
    setup_axis(axes[0], "Antes: posición actual de Ira")
    draw_open_path(axes[0], x_before, head, masha)

    setup_axis(axes[1], "Después: propuesta aceptada")
    draw_open_path(axes[1], x_after, head, masha)
    axes[1].scatter([x_before[head]], [head], s=68, facecolor="white",
                    edgecolor=MUTED, linewidth=1.7, zorder=5)
    axes[1].annotate("desplazamiento espacial", xy=(x_after[head], head),
                     xytext=(x_before[head] - 0.9, head + 2.0),
                     arrowprops={"arrowstyle": "->", "color": IRA}, color=IRA,
                     fontsize=8, ha="center")

    fig.suptitle("IraSpatialMove: mover Ira en el mismo slice", fontweight="bold")
    fig.text(0.5, 0.015, "Se propone un desplazamiento gaussiano en x e y; aquí se muestra la proyección x–τ.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.06, 1, 0.88))
    save_figure(fig, "worm_ira_spatial_move.png")


def make_swap():
    x_a = -1.0 + 0.2 * np.sin(2 * np.pi * TIME / 20 + 0.2)
    x_b = 1.0 + 0.2 * np.sin(2 * np.pi * TIME / 20 + 1.4)
    current_slice = 10
    next_slice = current_slice + 1

    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5), sharey=True)
    setup_axis(axes[0], "Antes: cabeza Ira en la trayectoria A")
    draw_track(axes[0], x_a, stop=current_slice, color=PATH_A)
    draw_track(axes[0], x_b, color=PATH_B)
    axes[0].plot(x_b[current_slice:next_slice + 1],
                 TIME[current_slice:next_slice + 1], color=MUTED,
                 linestyle="--", linewidth=2.0, zorder=3)
    mark_endpoint(axes[0], x_a[current_slice], current_slice, "Ira", IRA)
    axes[0].scatter([x_b[next_slice]], [next_slice], s=62, facecolor="white",
                    edgecolor=IRA, linewidth=1.6, zorder=5)
    axes[0].annotate("candidato k", (x_b[next_slice], next_slice),
                     xytext=(8, 8), textcoords="offset points", color=IRA, fontsize=8)

    setup_axis(axes[1], "Después: reconexión Swap")
    draw_track(axes[1], x_b, stop=current_slice, color=PATH_B)
    draw_track(axes[1], x_a, stop=current_slice, color=PATH_A)
    draw_track(axes[1], x_b, start=next_slice, color=PATH_A)
    axes[1].plot([x_a[current_slice], x_b[next_slice]],
                 [current_slice, next_slice], color=IRA, linewidth=2.4,
                 zorder=3)
    mark_endpoint(axes[1], x_b[current_slice], current_slice, "Ira nueva", IRA,
                  offset=(7, -14))
    axes[1].scatter([x_a[current_slice]], [current_slice], s=48, facecolor="white",
                    edgecolor=MUTED, linewidth=1.3, zorder=5)
    axes[1].annotate("enlace reconectado", xy=(x_a[current_slice] * 0.35 + x_b[next_slice] * 0.65,
                                                current_slice + 0.65),
                     xytext=(0.2, current_slice + 2.7),
                     arrowprops={"arrowstyle": "->", "color": IRA}, color=IRA,
                     fontsize=8, ha="center")

    fig.suptitle("Swap: cambiar el enlace de Ira a un candidato del siguiente slice",
                 fontweight="bold")
    fig.text(0.5, 0.015, "Esquema local; Masha no cambia y las demás conexiones del camino se omiten.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.06, 1, 0.88))
    save_figure(fig, "worm_swap.png")


def read_worm_snapshots(snapshot_file):
    events = {}
    current_move = None
    current_phase = None
    current_state = None

    with snapshot_file.open(encoding="utf-8") as source:
        for line in source:
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            if fields[0] == "EVENT":
                current_move, current_phase = fields[1:3]
                n_beads, n_particles = map(int, fields[3:5])
                current_state = {
                    "n_beads": n_beads,
                    "n_particles": n_particles,
                    "worm_open": int(fields[5]),
                    "ira": tuple(map(int, fields[6:8])),
                    "masha": tuple(map(int, fields[8:10])),
                    "x": np.full((n_beads, n_particles), np.nan),
                    "y": np.full((n_beads, n_particles), np.nan),
                    "next": np.full((n_beads, n_particles), -1, dtype=int),
                    "prev": np.full((n_beads, n_particles), -1, dtype=int),
                }
            elif fields[0] == "BEAD" and current_state is not None:
                bead, particle = map(int, fields[1:3])
                current_state["x"][bead, particle] = float(fields[3])
                current_state["y"][bead, particle] = float(fields[4])
                current_state["next"][bead, particle] = int(fields[5])
                current_state["prev"][bead, particle] = int(fields[6])
            elif fields[0] == "END" and current_state is not None:
                if not np.isfinite(current_state["x"]).all():
                    raise ValueError(f"Snapshot incompleto: {current_move} {current_phase}")
                events.setdefault(current_move, {})[current_phase] = current_state
                current_state = None

    if not events:
        raise ValueError(f"No se encontraron eventos Worm en {snapshot_file}")
    return events


def open_worm_edges(state):
    if not state["worm_open"]:
        return set()

    n_beads = state["n_beads"]
    bead, particle = state["masha"]
    ira = state["ira"]
    edges = set()
    visited = set()

    for _ in range(n_beads * state["n_particles"]):
        current = (bead, particle)
        if current == ira or current in visited:
            break
        visited.add(current)
        next_particle = state["next"][bead, particle]
        if next_particle < 0:
            break
        edges.add(current)
        bead = (bead + 1) % n_beads
        particle = next_particle
        if (bead, particle) == ira:
            break
    return edges


def draw_simulation_state(ax, state, title):
  n_beads = state["n_beads"]
  tau = np.arange(n_beads + 1) / n_beads
  open_edges = open_worm_edges(state)

  for bead in range(n_beads):
    next_bead = (bead + 1) % n_beads
    for particle in range(state["n_particles"]):
      next_particle = state["next"][bead, particle]
      if next_particle < 0:
        continue
      color = IRA if (bead, particle) in open_edges else PATH_A
      ax.plot([state["x"][bead, particle], state["x"][next_bead, next_particle]],
          [tau[bead], tau[bead + 1]], color=color,
          linewidth=1.35 if color == PATH_A else 2.0, zorder=2)

  bead_times = np.arange(n_beads) / n_beads
  ax.scatter(state["x"].ravel(), np.repeat(bead_times, state["n_particles"]),
         color=PATH_A, s=12, zorder=3)

  if state["worm_open"]:
    ira_bead, ira_particle = state["ira"]
    masha_bead, masha_particle = state["masha"]
    ax.scatter([state["x"][ira_bead, ira_particle]], [bead_times[ira_bead]],
           color=IRA, edgecolor="white", linewidth=0.7, s=54, zorder=5)
    ax.annotate("Ira", (state["x"][ira_bead, ira_particle], bead_times[ira_bead]),
          xytext=(8, 9), textcoords="offset points", color=IRA,
          fontsize=8, fontweight="bold")
    ax.scatter([state["x"][masha_bead, masha_particle]], [bead_times[masha_bead]],
           color=MASHA, edgecolor="white", linewidth=0.7, s=54, zorder=5)
    ax.annotate("Masha", (state["x"][masha_bead, masha_particle], bead_times[masha_bead]),
          xytext=(-52, -13), textcoords="offset points", color=MASHA,
          fontsize=8, fontweight="bold")

  xmin = float(np.min(state["x"]))
  xmax = float(np.max(state["x"]))
  padding = max(0.35, 0.12 * (xmax - xmin))
  ax.set_xlim(xmin - padding, xmax + padding)
  ax.set_ylim(0, 1)
  ax.set_xlabel("x (coordenada de la simulación)")
  ax.set_ylabel(r"$\tau/\beta$")
  ax.set_title(title, fontsize=10, fontweight="bold")
  ax.grid(True, color="#d9dde3", linewidth=0.7, alpha=0.7)
  ax.set_axisbelow(True)
  ax.tick_params(direction="in", top=True, right=True)


def make_simulation_figures(snapshot_file):
    events = read_worm_snapshots(snapshot_file)
    generated = 0
    for move, phases in events.items():
        if "before" not in phases or "after" not in phases:
            print(f"[omitido] {move}: falta el estado before/after")
            continue

        fig, axes = plt.subplots(1, 2, figsize=(9, 4.7), sharey=True)
        draw_simulation_state(axes[0], phases["before"], "Antes de la aceptación")
        draw_simulation_state(axes[1], phases["after"], "Después de la aceptación")
        fig.suptitle(f"{move}: transición aceptada de la simulación Worm", fontweight="bold")
        fig.text(0.5, 0.015,
                 "Coordenadas x y enlaces next/prev guardados por el algoritmo; se muestra la proyección x–τ.",
                 ha="center", fontsize=8)
        fig.tight_layout(rect=(0, 0.06, 1, 0.9))
        filename = f"worm_simulacion_{move.lower()}.png"
        save_figure(fig, filename)
        generated += 1

    if generated == 0:
        raise ValueError(f"No hay pares before/after completos en {snapshot_file}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--from-simulation",
        type=Path,
        metavar="SNAPSHOTS",
        help="Genera figuras before/after a partir de outworm_snapshots.dat",
    )
    args = parser.parse_args()

    if args.from_simulation is not None:
        make_simulation_figures(args.from_simulation)
        return

    print("=" * 76)
    print("GENERADOR DE ESQUEMAS DIDÁCTICOS DEL WORM ALGORITHM")
    print("=" * 76)

    print("\nGenerando figuras didácticas suaves con anotaciones Ira/Masha...")
    make_open_close()
    make_advance_recede()
    make_ira_spatial_move()
    make_swap()

    print("\n" + "=" * 76)
    print("¡Éxito! Todas las figuras han sido generadas y reemplazadas.")
    print(f"Ubicaciones actualizadas (PNG 300 DPI y PDF):")
    print(f"  1. {OUTPUT_DIR}")
    print(f"  2. {OVERLEAF_DIR}")
    print("=" * 76)


if __name__ == "__main__":
    main()
