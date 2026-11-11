#!/usr/bin/env python3
"""
Generador de figuras y diagramas de transiciones discretas del Worm Algorithm
basados en la formulación y esquemas de Tim Dornheim (PhD 2014).

Soporta dos modalidades:
  1. Por defecto: reconstruye con alta resolución y paleta diferencial los esquemas
     discretos de Dornheim (Open/Close, Advance/Recede, Bisección, Swap y Panel Completo)
     garantizando que la cadena completa esté siempre definida y que las anotaciones
     sean nítidas y legibles.
  2. Modo --from-simulation: lee snapshots before/after registrados por el código C.

Cada ejecución actualiza automáticamente los archivos PNG (300 DPI) y PDF en:
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
# PALETA DE COLORES DIFERENCIAL Y CONFIGURACIÓN VISUAL
# ==============================================================================

# Partícula 1 (cerrada / espectadora)
COLOR_P1 = "#1b7837"           # Verde esmeralda profundo
COLOR_P1_BEAD = "#2ca02c"      # Verde vivo para las cuentas de P1

# Partícula 2 (trayectoria de la partícula donde opera el gusano)
COLOR_P2 = "#d95f02"           # Naranja terracota intenso
COLOR_P2_BEAD = "#fd8d3c"      # Naranja cálido para cuentas de P2

# Resto de la cadena definida de la Partícula 2 (tramos fuera del worm activo)
COLOR_P2_REST = "#b25d00"      # Ocre tostado distintivo
COLOR_P2_REST_BEAD = "#e69f00" # Ámbar para cuentas del resto de cadena

# Modificaciones topológicas y remuestreos (Advance, Bisección, enlace cortado/añadido)
COLOR_MODIFIED = "#0066cc"     # Azul real vibrante
COLOR_MODIFIED_BEAD = "#3182bd"

# Permutación cuántica Swap (enlace de intercambio entre partículas)
COLOR_SWAP = "#7b1fa2"         # Púrpura intenso
COLOR_SWAP_BEAD = "#9c27b0"    # Violeta para cuentas del enlace swap

# Extremos cuánticos del gusano
COLOR_IRA = "#d90429"          # Rojo carmesí intenso (cabeza Ira)
COLOR_MASHA = "#008080"        # Teal / Verde azulado profundo (cola Masha)


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
    print(f"  [REEMPLAZADO] {filename:<42} -> Resultados/ y overleaf/imágenes/")


def clean_obsolete_didactic_files():
    """Elimina los archivos del antiguo 'estilo didáctico' en ambas carpetas de salida."""
    obsolete = [
        "worm_open_close",
        "worm_advance_recede",
        "worm_ira_spatial_move",
        "worm_swap",
    ]
    removed_count = 0
    for name in obsolete:
        for folder in [OUTPUT_DIR, OVERLEAF_DIR]:
            for ext in [".png", ".pdf"]:
                fpath = folder / f"{name}{ext}"
                if fpath.exists():
                    try:
                        fpath.unlink()
                        removed_count += 1
                    except OSError:
                        pass
    if removed_count > 0:
        print(f"  [LIMPIEZA] Eliminados {removed_count} archivos obsoletos del estilo didáctico.")


# ==============================================================================
# PARSEO DE ESQUEMAS DISCRETOS DE TIM DORNHEIM (PHD 2014)
# ==============================================================================

BOXES = [
    (50.300011, 22.603624, 196.052818, 168.40049),   # Subplot Izquierdo
    (266.299671, 22.603624, 412.052478, 168.40049)   # Subplot Derecho
]


def ensure_dornheim_svg(page_num):
    svg_path = Path(f"/tmp/page_{page_num}.svg")
    if not svg_path.exists():
        pdf_path = SCRIPT_DIR.parent / "Bibliografía" / "PhD Dornheim 2014.pdf"
        if pdf_path.exists():
            import subprocess
            subprocess.run(["pdftocairo", "-svg", "-f", str(page_num), "-l", str(page_num),
                            str(pdf_path), str(svg_path)], check=True)
    return svg_path


def to_dornheim_coords(px, py, box_idx):
    x0, y0, x1, y1 = BOXES[box_idx]
    x = -3.0 + 6.0 * (px - x0) / (x1 - x0)
    tau = 20.0 * (y1 - py) / (y1 - y0)
    return x, tau


def parse_dornheim_svg(svg_file):
    import re
    with open(svg_file) as f:
        text = f.read()

    matches = re.findall(r'<path[^>]+stroke=\"(rgb\([^)]+\))\"[^>]*d=\"([^\"]+)\"', text)
    paths = []
    for stroke, d in matches:
        if stroke == "rgb(0%, 0%, 0%)":
            continue
        if d.count('L') < 1:
            continue
        tokens = d.replace('M', ' ').replace('L', ' ').split()
        coords = [float(t) for t in tokens]
        pts = list(zip(coords[0::2], coords[1::2]))
        avg_x = sum(p[0] for p in pts) / len(pts)
        box_idx = 0 if avg_x < 230 else 1
        col = 'red' if '100%, 0%, 0%' in stroke else 'blue'
        converted = [to_dornheim_coords(px, py, box_idx) for px, py in pts]
        paths.append((box_idx, col, converted))

    circle_matches = re.findall(r'<path[^>]+fill=\"(rgb\([^)]+\))\"[^>]+d=\"M\s*([\d\.]+)\s+([\d\.]+)\s+C[^\"]+\"', text)
    beads = {0: {'red': [], 'blue': []}, 1: {'red': [], 'blue': []}}
    for fill, mx, my in circle_matches:
        px = float(mx) - 3.0
        py = float(my)
        box_idx = 0 if px < 230 else 1
        x, tau = to_dornheim_coords(px, py, box_idx)
        if '100%, 0%, 0%' in fill:
            col = 'red'
        elif '0%, 0%, 100%' in fill:
            col = 'blue'
        else:
            continue
        beads[box_idx][col].append((x, tau))

    return paths, beads


def setup_dornheim_axis(ax, xlabel="x", ylabel=r"$\tau/\epsilon$"):
    """Configuración tipográfica clara con marcas internas y rejilla suave de orientación."""
    ax.set_xlim(-3.0, 3.0)
    ax.set_ylim(0.0, 20.0)
    ax.set_xticks([-3.0, -1.5, 0.0, 1.5, 3.0])
    ax.set_yticks([0, 5, 10, 15, 20])
    ax.set_xlabel(xlabel, fontsize=13, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=13, fontweight="bold")
    ax.tick_params(direction='in', length=6, width=1.1, top=True, right=True, labelsize=11)
    ax.grid(True, linestyle=":", alpha=0.45, color="#b0b0b0")


def draw_annot_box(ax, x, y, text, color, fontsize=11.5, weight='bold', ha='left', va='center'):
    """Dibuja una anotación con fondo blanco redondeado para máxima legibilidad."""
    return ax.text(x, y, text, color=color, fontsize=fontsize, weight=weight, ha=ha, va=va,
                   bbox=dict(boxstyle="round,pad=0.24", facecolor="white", edgecolor=color, alpha=0.92, lw=0.85),
                   zorder=10)


# ==============================================================================
# FIGURAS DISCRETAS INDIVIDUALES (PÁGINAS 22, 24 Y 26)
# ==============================================================================

def make_dornheim_fig7_open_close():
    svg22 = ensure_dornheim_svg(22)
    paths22, _ = parse_dornheim_svg(svg22)

    fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.4), dpi=300)

    # Subplot 0: Insert / Remove
    ax = axs[0]
    setup_dornheim_axis(ax)
    xs1, ts1 = zip(*paths22[0][2])
    ax.plot(xs1, ts1, color=COLOR_P1, lw=2.2, zorder=2)
    ax.scatter(xs1, ts1, color=COLOR_P1_BEAD, edgecolor=COLOR_P1, s=46, zorder=3)

    xs_mod0, ts_mod0 = zip(*paths22[1][2])
    ax.plot(xs_mod0, ts_mod0, color=COLOR_MODIFIED, lw=3.0, zorder=4)
    ax.scatter(xs_mod0, ts_mod0, color=COLOR_MODIFIED_BEAD, edgecolor=COLOR_MODIFIED, s=58, zorder=5)

    ax.scatter([xs_mod0[-1]], [ts_mod0[-1]], color=COLOR_IRA, edgecolor='black', s=130, zorder=6, marker='o')
    ax.scatter([xs_mod0[0]], [ts_mod0[0]], color=COLOR_MASHA, edgecolor='black', s=130, zorder=6, marker='s')

    ax.set_title(r"$\mathbf{Insert\ /\ Remove}$", fontsize=14, fontweight="bold", pad=10)
    draw_annot_box(ax, -2.4, 18.2, "Partícula 1\n(Cerrada)", COLOR_P1, fontsize=11, ha="center")
    draw_annot_box(ax, 2.2, 18.2, "Gusano insertado\n(Sector G)", COLOR_MODIFIED, fontsize=11, ha="center")
    draw_annot_box(ax, xs_mod0[-1] + 0.35, ts_mod0[-1], r"$\mathcal{I}$ (Ira)", COLOR_IRA, fontsize=11.5)
    draw_annot_box(ax, xs_mod0[0] + 0.35, ts_mod0[0], r"$\mathcal{M}$ (Masha)", COLOR_MASHA, fontsize=11.5)

    # Subplot 1: Open / Close
    ax = axs[1]
    setup_dornheim_axis(ax)
    xs1, ts1 = zip(*paths22[2][2])
    ax.plot(xs1, ts1, color=COLOR_P1, lw=2.2, zorder=2)
    ax.scatter(xs1, ts1, color=COLOR_P1_BEAD, edgecolor=COLOR_P1, s=46, zorder=3)

    xs2, ts2 = zip(*paths22[3][2])
    ax.plot(xs2, ts2, color=COLOR_P2, lw=2.2, zorder=2)
    ax.scatter(xs2, ts2, color=COLOR_P2_BEAD, edgecolor=COLOR_P2, s=46, zorder=3)

    xs_mod1, ts_mod1 = zip(*paths22[4][2])
    ax.plot(xs_mod1, ts_mod1, color=COLOR_MODIFIED, lw=3.0, linestyle='--', zorder=4)
    ax.scatter(xs_mod1, ts_mod1, color=COLOR_MODIFIED_BEAD, edgecolor=COLOR_MODIFIED, s=58, zorder=5)

    ax.scatter([xs_mod1[-1]], [ts_mod1[-1]], color=COLOR_IRA, edgecolor='black', s=130, zorder=6, marker='o')
    ax.scatter([xs_mod1[0]], [ts_mod1[0]], color=COLOR_MASHA, edgecolor='black', s=130, zorder=6, marker='s')

    ax.set_title(r"$\mathbf{Open\ /\ Close}$ ($Z \leftrightarrow G$)", fontsize=14, fontweight="bold", pad=10)
    draw_annot_box(ax, -2.4, 18.2, "Partícula 1\n(Cerrada)", COLOR_P1, fontsize=11, ha="center")
    draw_annot_box(ax, 2.2, 18.2, "Partícula 2\n(Sector Z)", COLOR_P2, fontsize=11, ha="center")
    draw_annot_box(ax, xs_mod1[-1] + 0.35, ts_mod1[-1], r"$\mathcal{I}$ (Ira)", COLOR_IRA, fontsize=11.5)
    draw_annot_box(ax, xs_mod1[0] + 0.35, ts_mod1[0], r"$\mathcal{M}$ (Masha)", COLOR_MASHA, fontsize=11.5)
    draw_annot_box(ax, 2.2, 8.5, "Enlace cortado (Open)\no cerrado (Close)", COLOR_MODIFIED, fontsize=10.5, ha="center")

    plt.tight_layout()
    save_figure(fig, "dornheim_fig7_open_close.png")


def make_dornheim_fig8_advance_recede_bisection():
    svg24 = ensure_dornheim_svg(24)
    paths24, _ = parse_dornheim_svg(svg24)
    p2_full = paths24[4][2] # Trayectoria completa de la Partícula 2 (tau 0..20)

    fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.4), dpi=300)

    # Subplot 0: Advance / Recede (con resto de cadena definida)
    ax = axs[0]
    setup_dornheim_axis(ax)
    xs1, ts1 = zip(*paths24[0][2])
    ax.plot(xs1, ts1, color=COLOR_P1, lw=2.2, zorder=2)
    ax.scatter(xs1, ts1, color=COLOR_P1_BEAD, edgecolor=COLOR_P1, s=46, zorder=3)

    # Resto de la cadena definida de Partícula 2 (tau 0..2 y 12..20)
    rest_1 = p2_full[0:3]   # tau 0..2
    rest_2 = p2_full[12:21] # tau 12..20
    xs_r1, ts_r1 = zip(*rest_1)
    xs_r2, ts_r2 = zip(*rest_2)
    ax.plot(xs_r1, ts_r1, color=COLOR_P2_REST, lw=2.0, linestyle=':', zorder=2)
    ax.scatter(xs_r1, ts_r1, color=COLOR_P2_REST_BEAD, edgecolor=COLOR_P2_REST, s=38, zorder=3)
    ax.plot(xs_r2, ts_r2, color=COLOR_P2_REST, lw=2.0, linestyle=':', zorder=2)
    ax.scatter(xs_r2, ts_r2, color=COLOR_P2_REST_BEAD, edgecolor=COLOR_P2_REST, s=38, zorder=3)

    # Gusano preexistente (tau 2..8)
    worm_base = p2_full[2:9]
    xs_wb, ts_wb = zip(*worm_base)
    ax.plot(xs_wb, ts_wb, color=COLOR_P2, lw=2.6, zorder=3)
    ax.scatter(xs_wb, ts_wb, color=COLOR_P2_BEAD, edgecolor=COLOR_P2, s=50, zorder=4)

    # Segmento de Advance nuevo (tau 8..12)
    xs_adv, ts_adv = zip(*paths24[2][2])
    ax.plot(xs_adv, ts_adv, color=COLOR_MODIFIED, lw=3.2, zorder=4)
    ax.scatter(xs_adv, ts_adv, color=COLOR_MODIFIED_BEAD, edgecolor=COLOR_MODIFIED, s=60, zorder=5)

    ira_old = (xs_adv[0], ts_adv[0])
    ira_new = (xs_adv[-1], ts_adv[-1])
    masha = (xs_wb[0], ts_wb[0])

    ax.scatter([ira_new[0]], [ira_new[1]], color=COLOR_MODIFIED, edgecolor='black', s=140, zorder=6, marker='o')
    ax.scatter([ira_old[0]], [ira_old[1]], color=COLOR_IRA, edgecolor='black', s=120, zorder=6, marker='o')
    ax.scatter([masha[0]], [masha[1]], color=COLOR_MASHA, edgecolor='black', s=130, zorder=6, marker='s')

    ax.set_title(r"$\mathbf{Advance\ /\ Recede}$ (Sector $G$)", fontsize=14, fontweight="bold", pad=10)
    draw_annot_box(ax, -2.4, 18.2, "Partícula 1\n(Cerrada)", COLOR_P1, fontsize=11, ha="center")
    draw_annot_box(ax, ira_new[0] + 0.35, ira_new[1], r"$\mathcal{I}_{\mathrm{new}}$ (Advance $\uparrow$)", COLOR_MODIFIED, fontsize=11.5)
    draw_annot_box(ax, ira_old[0] + 0.35, ira_old[1] - 0.5, r"$\mathcal{I}_{\mathrm{old}}$ (Recede $\downarrow$)", COLOR_IRA, fontsize=11.5)
    draw_annot_box(ax, masha[0] - 0.4, masha[1] - 0.2, r"$\mathcal{M}$ (Masha)", COLOR_MASHA, fontsize=11.5, ha="right")
    draw_annot_box(ax, 2.25, 16.5, "Resto de cadena\ndefinida (Partícula 2)", COLOR_P2_REST, fontsize=10.5, ha="center")
    draw_annot_box(ax, 2.35, 10.0, "Cuentas añadidas\n(Advance azul)", COLOR_MODIFIED, fontsize=10.5, ha="center")

    # Subplot 1: Bisection
    ax = axs[1]
    setup_dornheim_axis(ax)
    xs1, ts1 = zip(*paths24[3][2])
    ax.plot(xs1, ts1, color=COLOR_P1, lw=2.2, zorder=2)
    ax.scatter(xs1, ts1, color=COLOR_P1_BEAD, edgecolor=COLOR_P1, s=46, zorder=3)

    xs2, ts2 = zip(*paths24[4][2])
    ax.plot(xs2, ts2, color=COLOR_P2, lw=2.2, zorder=2)
    ax.scatter(xs2, ts2, color=COLOR_P2_BEAD, edgecolor=COLOR_P2, s=46, zorder=3)

    xs_bis, ts_bis = zip(*paths24[5][2])
    ax.plot(xs_bis, ts_bis, color=COLOR_MODIFIED, lw=3.0, linestyle='--', zorder=4)
    ax.scatter(xs_bis, ts_bis, color=COLOR_MODIFIED_BEAD, edgecolor=COLOR_MODIFIED, s=58, zorder=5)

    ax.set_title(r"$\mathbf{Bisección\ /\ IraSpatialMove}$", fontsize=14, fontweight="bold", pad=10)
    draw_annot_box(ax, -2.4, 18.2, "Partícula 1", COLOR_P1, fontsize=11, ha="center")
    draw_annot_box(ax, 2.2, 18.2, "Partícula 2", COLOR_P2, fontsize=11, ha="center")
    draw_annot_box(ax, 2.25, 9.5, "Segmento interno\nremuestreado (azul)", COLOR_MODIFIED, fontsize=10.5, ha="center")

    plt.tight_layout()
    save_figure(fig, "dornheim_fig8_advance_recede_bisection.png")


def make_dornheim_fig9_swap():
    svg24 = ensure_dornheim_svg(24)
    svg26 = ensure_dornheim_svg(26)
    paths24, _ = parse_dornheim_svg(svg24)
    paths26, _ = parse_dornheim_svg(svg26)
    p2_full = paths24[4][2] # Trayectoria completa Partícula 2

    fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.4), dpi=300)

    # Subplot 0: Antes del Swap
    ax = axs[0]
    setup_dornheim_axis(ax)
    xs1, ts1 = zip(*paths26[0][2])
    ax.plot(xs1, ts1, color=COLOR_P1, lw=2.2, zorder=2)
    ax.scatter(xs1, ts1, color=COLOR_P1_BEAD, edgecolor=COLOR_P1, s=46, zorder=3)

    rest_1 = p2_full[0:3]   # tau 0..2
    rest_2 = p2_full[12:21] # tau 12..20
    xs_r1, ts_r1 = zip(*rest_1)
    xs_r2, ts_r2 = zip(*rest_2)
    ax.plot(xs_r1, ts_r1, color=COLOR_P2_REST, lw=2.0, linestyle=':', zorder=2)
    ax.scatter(xs_r1, ts_r1, color=COLOR_P2_REST_BEAD, edgecolor=COLOR_P2_REST, s=38, zorder=3)
    ax.plot(xs_r2, ts_r2, color=COLOR_P2_REST, lw=2.0, linestyle=':', zorder=2)
    ax.scatter(xs_r2, ts_r2, color=COLOR_P2_REST_BEAD, edgecolor=COLOR_P2_REST, s=38, zorder=3)

    worm_base = p2_full[2:13] # tau 2..12
    xs_wb, ts_wb = zip(*worm_base)
    ax.plot(xs_wb, ts_wb, color=COLOR_P2, lw=2.6, zorder=3)
    ax.scatter(xs_wb, ts_wb, color=COLOR_P2_BEAD, edgecolor=COLOR_P2, s=50, zorder=4)

    ax.scatter([xs_wb[-1]], [ts_wb[-1]], color=COLOR_IRA, edgecolor='black', s=130, zorder=6, marker='o')
    ax.scatter([xs_wb[0]], [ts_wb[0]], color=COLOR_MASHA, edgecolor='black', s=130, zorder=6, marker='s')

    ax.set_title(r"$\mathbf{Antes\ del\ Swap}$ (Sector $G$)", fontsize=14, fontweight="bold", pad=10)
    draw_annot_box(ax, -2.4, 18.2, "Partícula 1\n(Cerrada)", COLOR_P1, fontsize=11, ha="center")
    draw_annot_box(ax, xs_wb[-1] + 0.35, ts_wb[-1], r"$\mathcal{I}$ (Ira en Part. 2)", COLOR_IRA, fontsize=11.5)
    draw_annot_box(ax, xs_wb[0] + 0.35, ts_wb[0], r"$\mathcal{M}$ (Masha)", COLOR_MASHA, fontsize=11.5)
    draw_annot_box(ax, 2.25, 16.5, "Resto de cadena\ndefinida (Partícula 2)", COLOR_P2_REST, fontsize=10.5, ha="center")

    # Subplot 1: Después del Swap
    ax = axs[1]
    setup_dornheim_axis(ax)
    xs_p1_a, ts_p1_a = zip(*paths26[3][2]) # tau 0..12
    ax.plot(xs_p1_a, ts_p1_a, color=COLOR_P1, lw=2.2, zorder=2)
    ax.scatter(xs_p1_a, ts_p1_a, color=COLOR_P1_BEAD, edgecolor=COLOR_P1, s=46, zorder=3)

    xs_p1_b, ts_p1_b = zip(*paths26[4][2]) # tau 16..20
    ax.plot(xs_p1_b, ts_p1_b, color=COLOR_P1, lw=2.2, zorder=2)
    ax.scatter(xs_p1_b, ts_p1_b, color=COLOR_P1_BEAD, edgecolor=COLOR_P1, s=46, zorder=3)

    p5 = paths26[5][2]
    worm_p2 = p5[0:11]    # tau 2..12
    swap_link = p5[10:15] # tau 12..16 cross
    xs_w2, ts_w2 = zip(*worm_p2)
    xs_sw, ts_sw = zip(*swap_link)

    ax.plot(xs_w2, ts_w2, color=COLOR_P2, lw=2.6, zorder=3)
    ax.scatter(xs_w2, ts_w2, color=COLOR_P2_BEAD, edgecolor=COLOR_P2, s=50, zorder=4)

    ax.plot(xs_sw, ts_sw, color=COLOR_SWAP, lw=3.4, zorder=5)
    ax.scatter(xs_sw, ts_sw, color=COLOR_SWAP_BEAD, edgecolor=COLOR_SWAP, s=62, zorder=6)

    # Resto de cadena Partícula 2
    ax.plot(xs_r1, ts_r1, color=COLOR_P2_REST, lw=2.0, linestyle=':', zorder=2)
    ax.scatter(xs_r1, ts_r1, color=COLOR_P2_REST_BEAD, edgecolor=COLOR_P2_REST, s=38, zorder=3)
    ax.plot(xs_r2, ts_r2, color=COLOR_P2_REST, lw=2.0, linestyle=':', zorder=2)
    ax.scatter(xs_r2, ts_r2, color=COLOR_P2_REST_BEAD, edgecolor=COLOR_P2_REST, s=38, zorder=3)

    ira_new_sw = (xs_p1_a[-1], ts_p1_a[-1])
    ira_old_sw = (xs_w2[-1], ts_w2[-1])
    masha_sw = (xs_w2[0], ts_w2[0])

    ax.scatter([ira_new_sw[0]], [ira_new_sw[1]], color=COLOR_IRA, edgecolor='black', s=140, zorder=7, marker='o')
    ax.scatter([ira_old_sw[0]], [ira_old_sw[1]], color=COLOR_SWAP, edgecolor='black', s=120, zorder=7, marker='o')
    ax.scatter([masha_sw[0]], [masha_sw[1]], color=COLOR_MASHA, edgecolor='black', s=130, zorder=7, marker='s')

    ax.set_title(r"$\mathbf{Después\ del\ Swap}$ (Permutación Cuántica)", fontsize=14, fontweight="bold", pad=10)
    draw_annot_box(ax, 0.3, 14.5, r"Enlace Swap $\longrightarrow$", COLOR_SWAP, fontsize=12, ha="center")
    draw_annot_box(ax, ira_new_sw[0] - 0.4, ira_new_sw[1] + 0.8, r"$\mathcal{I}_{\mathrm{new}}$ (Nueva Ira)", COLOR_IRA, fontsize=11.5, ha="right")
    draw_annot_box(ax, ira_old_sw[0] + 0.35, ira_old_sw[1] - 0.5, r"$\mathcal{I}_{\mathrm{old}}$ (Vieja Ira)", COLOR_SWAP, fontsize=11.5)
    draw_annot_box(ax, masha_sw[0] + 0.35, masha_sw[1], r"$\mathcal{M}$ (Masha)", COLOR_MASHA, fontsize=11.5)

    plt.tight_layout()
    save_figure(fig, "dornheim_fig9_swap.png")


# ==============================================================================
# PANEL COMPLETO 2x2: TODOS LOS MOVIMIENTOS WORM
# ==============================================================================

def make_dornheim_panel_completo():
    svg22 = ensure_dornheim_svg(22)
    svg24 = ensure_dornheim_svg(24)
    svg26 = ensure_dornheim_svg(26)
    paths22, _ = parse_dornheim_svg(svg22)
    paths24, _ = parse_dornheim_svg(svg24)
    paths26, _ = parse_dornheim_svg(svg26)
    p2_full = paths24[4][2] # Trayectoria completa Partícula 2 (tau 0..20)

    fig, axs = plt.subplots(2, 2, figsize=(13.5, 11.5), dpi=300)

    # --------------------------------------------------------------------------
    # (0, 0): Open / Close
    # --------------------------------------------------------------------------
    ax = axs[0, 0]
    setup_dornheim_axis(ax)
    xs1, ts1 = zip(*paths22[2][2])
    ax.plot(xs1, ts1, color=COLOR_P1, lw=2.2, zorder=2)
    ax.scatter(xs1, ts1, color=COLOR_P1_BEAD, edgecolor=COLOR_P1, s=44, zorder=3)

    xs2, ts2 = zip(*paths22[3][2])
    ax.plot(xs2, ts2, color=COLOR_P2, lw=2.2, zorder=2)
    ax.scatter(xs2, ts2, color=COLOR_P2_BEAD, edgecolor=COLOR_P2, s=44, zorder=3)

    xs_mod, ts_mod = zip(*paths22[4][2])
    ax.plot(xs_mod, ts_mod, color=COLOR_MODIFIED, lw=3.0, linestyle='--', zorder=4)
    ax.scatter(xs_mod, ts_mod, color=COLOR_MODIFIED_BEAD, edgecolor=COLOR_MODIFIED, s=56, zorder=5)

    ax.scatter([xs_mod[-1]], [ts_mod[-1]], color=COLOR_IRA, edgecolor='black', s=130, zorder=6, marker='o')
    ax.scatter([xs_mod[0]], [ts_mod[0]], color=COLOR_MASHA, edgecolor='black', s=130, zorder=6, marker='s')

    ax.set_title(r"$\mathbf{(a)\ Movimiento\ Open\ /\ Close}$ ($Z \leftrightarrow G$)", fontsize=14, fontweight="bold", pad=10)
    draw_annot_box(ax, xs_mod[-1] + 0.35, ts_mod[-1], r"$\mathcal{I}$ (Ira)", COLOR_IRA, fontsize=11.5)
    draw_annot_box(ax, xs_mod[0] + 0.35, ts_mod[0], r"$\mathcal{M}$ (Masha)", COLOR_MASHA, fontsize=11.5)
    draw_annot_box(ax, -2.4, 18.2, "Partícula 1\n(Cerrada)", COLOR_P1, fontsize=11, ha="center")
    draw_annot_box(ax, 2.2, 18.2, "Partícula 2\n(Sector Z/G)", COLOR_P2, fontsize=11, ha="center")
    draw_annot_box(ax, 2.25, 8.5, "Enlace cortado (Open)\no cerrado (Close)", COLOR_MODIFIED, fontsize=10.5, ha="center")

    # --------------------------------------------------------------------------
    # (0, 1): Advance / Recede (con cadena completa definida)
    # --------------------------------------------------------------------------
    ax = axs[0, 1]
    setup_dornheim_axis(ax)
    xs1, ts1 = zip(*paths24[0][2])
    ax.plot(xs1, ts1, color=COLOR_P1, lw=2.2, zorder=2)
    ax.scatter(xs1, ts1, color=COLOR_P1_BEAD, edgecolor=COLOR_P1, s=44, zorder=3)

    rest_1 = p2_full[0:3]
    rest_2 = p2_full[12:21]
    xs_r1, ts_r1 = zip(*rest_1)
    xs_r2, ts_r2 = zip(*rest_2)
    ax.plot(xs_r1, ts_r1, color=COLOR_P2_REST, lw=2.0, linestyle=':', zorder=2)
    ax.scatter(xs_r1, ts_r1, color=COLOR_P2_REST_BEAD, edgecolor=COLOR_P2_REST, s=38, zorder=3)
    ax.plot(xs_r2, ts_r2, color=COLOR_P2_REST, lw=2.0, linestyle=':', zorder=2)
    ax.scatter(xs_r2, ts_r2, color=COLOR_P2_REST_BEAD, edgecolor=COLOR_P2_REST, s=38, zorder=3)

    worm_base = p2_full[2:9]
    xs_wb, ts_wb = zip(*worm_base)
    ax.plot(xs_wb, ts_wb, color=COLOR_P2, lw=2.6, zorder=3)
    ax.scatter(xs_wb, ts_wb, color=COLOR_P2_BEAD, edgecolor=COLOR_P2, s=50, zorder=4)

    xs_adv, ts_adv = zip(*paths24[2][2])
    ax.plot(xs_adv, ts_adv, color=COLOR_MODIFIED, lw=3.2, zorder=4)
    ax.scatter(xs_adv, ts_adv, color=COLOR_MODIFIED_BEAD, edgecolor=COLOR_MODIFIED, s=60, zorder=5)

    ira_old = (xs_adv[0], ts_adv[0])
    ira_new = (xs_adv[-1], ts_adv[-1])
    masha = (xs_wb[0], ts_wb[0])

    ax.scatter([ira_new[0]], [ira_new[1]], color=COLOR_MODIFIED, edgecolor='black', s=140, zorder=6, marker='o')
    ax.scatter([ira_old[0]], [ira_old[1]], color=COLOR_IRA, edgecolor='black', s=120, zorder=6, marker='o')
    ax.scatter([masha[0]], [masha[1]], color=COLOR_MASHA, edgecolor='black', s=130, zorder=6, marker='s')

    ax.set_title(r"$\mathbf{(b)\ Movimiento\ Advance\ /\ Recede}$ (Sector $G$)", fontsize=14, fontweight="bold", pad=10)
    draw_annot_box(ax, ira_new[0] + 0.35, ira_new[1], r"$\mathcal{I}_{\mathrm{new}}$ (Advance $\uparrow$)", COLOR_MODIFIED, fontsize=11.5)
    draw_annot_box(ax, ira_old[0] + 0.35, ira_old[1] - 0.5, r"$\mathcal{I}_{\mathrm{old}}$ (Recede $\downarrow$)", COLOR_IRA, fontsize=11.5)
    draw_annot_box(ax, masha[0] - 0.4, masha[1] - 0.2, r"$\mathcal{M}$ (Masha)", COLOR_MASHA, fontsize=11.5, ha="right")
    draw_annot_box(ax, 2.25, 16.5, "Resto de cadena\ndefinida (Partícula 2)", COLOR_P2_REST, fontsize=10.5, ha="center")
    draw_annot_box(ax, 2.35, 10.0, "Cuentas añadidas\n(Advance azul)", COLOR_MODIFIED, fontsize=10.5, ha="center")

    # --------------------------------------------------------------------------
    # (1, 0): Bisección / IraSpatialMove
    # --------------------------------------------------------------------------
    ax = axs[1, 0]
    setup_dornheim_axis(ax)
    xs1, ts1 = zip(*paths24[3][2])
    ax.plot(xs1, ts1, color=COLOR_P1, lw=2.2, zorder=2)
    ax.scatter(xs1, ts1, color=COLOR_P1_BEAD, edgecolor=COLOR_P1, s=44, zorder=3)

    xs2, ts2 = zip(*paths24[4][2])
    ax.plot(xs2, ts2, color=COLOR_P2, lw=2.2, zorder=2)
    ax.scatter(xs2, ts2, color=COLOR_P2_BEAD, edgecolor=COLOR_P2, s=44, zorder=3)

    xs_bis, ts_bis = zip(*paths24[5][2])
    ax.plot(xs_bis, ts_bis, color=COLOR_MODIFIED, lw=3.0, linestyle='--', zorder=4)
    ax.scatter(xs_bis, ts_bis, color=COLOR_MODIFIED_BEAD, edgecolor=COLOR_MODIFIED, s=58, zorder=5)

    ax.set_title(r"$\mathbf{(c)\ Bisección\ /\ IraSpatialMove}$", fontsize=14, fontweight="bold", pad=10)
    draw_annot_box(ax, -2.4, 18.2, "Partícula 1", COLOR_P1, fontsize=11, ha="center")
    draw_annot_box(ax, 2.2, 18.2, "Partícula 2", COLOR_P2, fontsize=11, ha="center")
    draw_annot_box(ax, 2.25, 9.5, "Segmento interno\nremuestreado (azul)", COLOR_MODIFIED, fontsize=10.5, ha="center")

    # --------------------------------------------------------------------------
    # (1, 1): Movimiento Swap (Permutación Cuántica)
    # --------------------------------------------------------------------------
    ax = axs[1, 1]
    setup_dornheim_axis(ax)
    xs_p1_a, ts_p1_a = zip(*paths26[3][2])
    ax.plot(xs_p1_a, ts_p1_a, color=COLOR_P1, lw=2.2, zorder=2)
    ax.scatter(xs_p1_a, ts_p1_a, color=COLOR_P1_BEAD, edgecolor=COLOR_P1, s=44, zorder=3)

    xs_p1_b, ts_p1_b = zip(*paths26[4][2])
    ax.plot(xs_p1_b, ts_p1_b, color=COLOR_P1, lw=2.2, zorder=2)
    ax.scatter(xs_p1_b, ts_p1_b, color=COLOR_P1_BEAD, edgecolor=COLOR_P1, s=44, zorder=3)

    p5 = paths26[5][2]
    worm_p2 = p5[0:11]
    swap_link = p5[10:15]
    xs_w2, ts_w2 = zip(*worm_p2)
    xs_sw, ts_sw = zip(*swap_link)

    ax.plot(xs_w2, ts_w2, color=COLOR_P2, lw=2.6, zorder=3)
    ax.scatter(xs_w2, ts_w2, color=COLOR_P2_BEAD, edgecolor=COLOR_P2, s=50, zorder=4)

    ax.plot(xs_sw, ts_sw, color=COLOR_SWAP, lw=3.4, zorder=5)
    ax.scatter(xs_sw, ts_sw, color=COLOR_SWAP_BEAD, edgecolor=COLOR_SWAP, s=62, zorder=6)

    # Resto de cadena Partícula 2
    ax.plot(xs_r1, ts_r1, color=COLOR_P2_REST, lw=2.0, linestyle=':', zorder=2)
    ax.scatter(xs_r1, ts_r1, color=COLOR_P2_REST_BEAD, edgecolor=COLOR_P2_REST, s=38, zorder=3)
    ax.plot(xs_r2, ts_r2, color=COLOR_P2_REST, lw=2.0, linestyle=':', zorder=2)
    ax.scatter(xs_r2, ts_r2, color=COLOR_P2_REST_BEAD, edgecolor=COLOR_P2_REST, s=38, zorder=3)

    ira_new_sw = (xs_p1_a[-1], ts_p1_a[-1])
    ira_old_sw = (xs_w2[-1], ts_w2[-1])
    masha_sw = (xs_w2[0], ts_w2[0])

    ax.scatter([ira_new_sw[0]], [ira_new_sw[1]], color=COLOR_IRA, edgecolor='black', s=140, zorder=7, marker='o')
    ax.scatter([ira_old_sw[0]], [ira_old_sw[1]], color=COLOR_SWAP, edgecolor='black', s=120, zorder=7, marker='o')
    ax.scatter([masha_sw[0]], [masha_sw[1]], color=COLOR_MASHA, edgecolor='black', s=130, zorder=7, marker='s')

    ax.set_title(r"$\mathbf{(d)\ Movimiento\ Swap}$ (Permutación Cuántica)", fontsize=14, fontweight="bold", pad=10)
    draw_annot_box(ax, 0.3, 14.5, r"Enlace Swap $\longrightarrow$", COLOR_SWAP, fontsize=12, ha="center")
    draw_annot_box(ax, ira_new_sw[0] - 0.4, ira_new_sw[1] + 0.8, r"$\mathcal{I}_{\mathrm{new}}$ (Nueva Ira)", COLOR_IRA, fontsize=11.5, ha="right")
    draw_annot_box(ax, ira_old_sw[0] + 0.35, ira_old_sw[1] - 0.5, r"$\mathcal{I}_{\mathrm{old}}$ (Vieja Ira)", COLOR_SWAP, fontsize=11.5)
    draw_annot_box(ax, masha_sw[0] + 0.35, masha_sw[1], r"$\mathcal{M}$ (Masha)", COLOR_MASHA, fontsize=11.5)

    plt.tight_layout()
    save_figure(fig, "figura_panel_completo_worm.png")


# ==============================================================================
# MODALIDAD SIMULACIÓN C (PROYECCIONES REALES BEFORE/AFTER)
# ==============================================================================

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
            is_open = (bead, particle) in open_edges
            color = COLOR_MODIFIED if is_open else COLOR_P1
            lw = 2.4 if is_open else 1.5
            ax.plot([state["x"][bead, particle], state["x"][next_bead, next_particle]],
                    [tau[bead], tau[bead + 1]], color=color, linewidth=lw, zorder=2)

    bead_times = np.arange(n_beads) / n_beads
    ax.scatter(state["x"].ravel(), np.repeat(bead_times, state["n_particles"]),
               color=COLOR_P1_BEAD, edgecolor=COLOR_P1, s=18, zorder=3)

    if state["worm_open"]:
        ira_bead, ira_particle = state["ira"]
        masha_bead, masha_particle = state["masha"]
        ax.scatter([state["x"][ira_bead, ira_particle]], [bead_times[ira_bead]],
                   color=COLOR_IRA, edgecolor="black", linewidth=1.0, s=70, zorder=5)
        draw_annot_box(ax, state["x"][ira_bead, ira_particle] + 0.1, bead_times[ira_bead],
                       r"$\mathcal{I}$ (Ira)", COLOR_IRA, fontsize=10.5)

        ax.scatter([state["x"][masha_bead, masha_particle]], [bead_times[masha_bead]],
                   color=COLOR_MASHA, edgecolor="black", linewidth=1.0, s=70, zorder=5, marker="s")
        draw_annot_box(ax, state["x"][masha_bead, masha_particle] - 0.1, bead_times[masha_bead],
                       r"$\mathcal{M}$ (Masha)", COLOR_MASHA, fontsize=10.5, ha="right")

    xmin = float(np.min(state["x"]))
    xmax = float(np.max(state["x"]))
    padding = max(0.35, 0.12 * (xmax - xmin))
    ax.set_xlim(xmin - padding, xmax + padding)
    ax.set_ylim(0, 1)
    ax.set_xlabel("x (coordenada de la simulación)", fontsize=12, fontweight="bold")
    ax.set_ylabel(r"$\tau/\beta$", fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.45, color="#b0b0b0")
    ax.set_axisbelow(True)
    ax.tick_params(direction="in", top=True, right=True, labelsize=10.5)


def make_simulation_figures(snapshot_file):
    events = read_worm_snapshots(snapshot_file)
    generated = 0
    for move, phases in events.items():
        if "before" not in phases or "after" not in phases:
            print(f"[omitido] {move}: falta el estado before/after")
            continue

        fig, axes = plt.subplots(1, 2, figsize=(10, 5.0), sharey=True)
        draw_simulation_state(axes[0], phases["before"], "Antes de la aceptación")
        draw_simulation_state(axes[1], phases["after"], "Después de la aceptación")
        fig.suptitle(f"{move}: transición aceptada de la simulación Worm", fontweight="bold", fontsize=13)
        fig.text(0.5, 0.015,
                 "Coordenadas x y enlaces next/prev guardados por el algoritmo; se muestra la proyección x–τ.",
                 ha="center", fontsize=9)
        fig.tight_layout(rect=(0, 0.06, 1, 0.92))
        filename = f"worm_simulacion_{move.lower()}.png"
        save_figure(fig, filename)
        generated += 1

    if generated == 0:
        raise ValueError(f"No hay pares before/after completos en {snapshot_file}")


# ==============================================================================
# PUNTO DE ENTRADA PRINCIPAL
# ==============================================================================

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
    print("GENERADOR DE ESQUEMAS WORM ALGORITHM: FORMATO TIM DORNHEIM")
    print("=" * 76)

    # Limpiar figuras obsoletas del antiguo estilo didáctico
    clean_obsolete_didactic_files()

    print("\nGenerando figuras discretas de Tim Dornheim con paleta diferencial y alta resolución...")
    make_dornheim_fig7_open_close()
    make_dornheim_fig8_advance_recede_bisection()
    make_dornheim_fig9_swap()
    make_dornheim_panel_completo()

    print("\n" + "=" * 76)
    print("¡Éxito! Todas las figuras han sido generadas y reemplazadas.")
    print("Ubicaciones actualizadas (PNG 300 DPI y PDF):")
    print(f"  1. {OUTPUT_DIR}")
    print(f"  2. {OVERLEAF_DIR}")
    print("=" * 76)


if __name__ == "__main__":
    main()
