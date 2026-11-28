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
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt
import numpy as np

# Rutas de salida
SCRIPT_DIR = Path(__file__).resolve().parent.parent
OUTPUT_DIR = SCRIPT_DIR / "Resultados" / "Figuras_Worm_Dornheim"
OVERLEAF_DIR = SCRIPT_DIR.parent / "overleaf" / "imágenes"
DEFAULT_DORNHEIM_PDF = SCRIPT_DIR.parent / "Bibliografía" / "PhD Dornheim 2014.pdf"

_DORNHEIM_CACHE = {}
_DORNHEIM_FALLBACK_WARNED = False


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
    """Guarda en Resultados (PNG 300 DPI y PDF) y copia a Overleaf si es accesible."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    png_out = OUTPUT_DIR / filename
    pdf_out = png_out.with_suffix(".pdf")

    # Renderizar una sola vez y copiar a Overleaf si es posible
    fig.savefig(png_out, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_out, bbox_inches="tight")

    copied_overleaf = False
    try:
        OVERLEAF_DIR.mkdir(parents=True, exist_ok=True)
        shutil.copy2(png_out, OVERLEAF_DIR / png_out.name)
        shutil.copy2(pdf_out, OVERLEAF_DIR / pdf_out.name)
        copied_overleaf = True
    except OSError as exc:
        print(f"  [AVISO] No se pudo exportar a {OVERLEAF_DIR}: {exc}")

    plt.close(fig)
    dest = "Resultados/ y overleaf/imágenes/" if copied_overleaf else "Resultados/"
    print(f"  [REEMPLAZADO] {filename:<42} -> {dest}")
    return copied_overleaf


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


def _build_fallback_dornheim_paths(page_num):
    """Genera trayectorias esquemáticas autocontenidas si falta el PDF externo."""
    tau = np.linspace(0.0, 20.0, 21)
    p1 = list(zip(-1.25 + 0.30 * np.sin(0.25 * np.pi * tau / 2.5 + 0.3), tau))
    p2_full = list(zip(1.10 + 0.32 * np.cos(0.22 * np.pi * tau / 2.2), tau))
    beads = {0: {"red": [], "blue": []}, 1: {"red": [], "blue": []}}

    if page_num == 22:
        p_ins = list(zip(1.05 + 0.22 * np.sin(0.5 * (tau[4:13] - 4.0)), tau[4:13]))
        cut_link = p2_full[9:12]
        return [(0, "blue", p1), (0, "red", p_ins), (1, "blue", p1), (1, "blue", p2_full), (1, "red", cut_link)], beads

    if page_num == 24:
        adv_seg = p2_full[8:13]
        bis_x = np.array([p2_full[k][0] for k in range(6, 14)], dtype=float)
        bis_x[1:-1] += np.array([-0.24, -0.35, -0.28, -0.10, 0.14, 0.22])
        bis_seg = list(zip(bis_x, tau[6:14]))
        return [
            (0, "blue", p1),
            (0, "blue", p2_full[2:9]),
            (0, "red", adv_seg),
            (1, "blue", p1),
            (1, "blue", p2_full),
            (1, "red", bis_seg),
        ], beads

    if page_num == 26:
        p1_low = p1[0:13]
        p1_high = p1[16:21]
        bridge_tau = tau[13:17]
        bridge_x = np.linspace(p2_full[12][0], p1_high[0][0], 5)[1:]
        p5 = list(p2_full[2:13]) + list(zip(bridge_x, bridge_tau))
        return [
            (0, "blue", p1),
            (0, "blue", p2_full),
            (0, "red", p2_full[2:13]),
            (1, "blue", p1_low),
            (1, "blue", p1_high),
            (1, "red", p5),
        ], beads

    raise ValueError(f"Página Dornheim no soportada: {page_num}")


def ensure_dornheim_svg(page_num, pdf_path=None):
    """Localiza o genera el SVG de la página de Dornheim si existen PDF y pdftocairo."""
    legacy_svg = Path(tempfile.gettempdir()) / f"page_{page_num}.svg"
    if legacy_svg.is_file():
        return legacy_svg

    cached_svg = Path(tempfile.gettempdir()) / f"dornheim_page_{page_num}.svg"
    if cached_svg.is_file():
        return cached_svg

    pdf = Path(pdf_path) if pdf_path is not None else Path(
        os.environ.get("DORNHEIM_PDF", str(DEFAULT_DORNHEIM_PDF))
    )
    if pdf.is_file() and shutil.which("pdftocairo"):
        try:
            subprocess.run(
                ["pdftocairo", "-svg", "-f", str(page_num), "-l", str(page_num),
                 str(pdf), str(cached_svg)],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            if cached_svg.is_file():
                return cached_svg
        except (OSError, subprocess.SubprocessError):
            pass
    return None


def to_dornheim_coords(px, py, box_idx):
    x0, y0, x1, y1 = BOXES[box_idx]
    x = -3.0 + 6.0 * (px - x0) / (x1 - x0)
    tau = 20.0 * (y1 - py) / (y1 - y0)
    return x, tau


def parse_dornheim_svg(svg_file):
    text = Path(svg_file).read_text(encoding="utf-8", errors="ignore")

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


def load_dornheim_paths(page_num, pdf_path=None):
    """Carga y cachea las trayectorias SVG de Dornheim o usa fallback autocontenido."""
    global _DORNHEIM_FALLBACK_WARNED
    cache_key = (page_num, str(pdf_path) if pdf_path else "")
    if cache_key in _DORNHEIM_CACHE:
        return _DORNHEIM_CACHE[cache_key]

    svg_path = ensure_dornheim_svg(page_num, pdf_path=pdf_path)
    if svg_path is not None and svg_path.is_file():
        parsed = parse_dornheim_svg(svg_path)
    else:
        if not _DORNHEIM_FALLBACK_WARNED:
            print("  [AVISO] No se encontró el PDF de Dornheim o pdftocairo; usando trayectorias esquemáticas autocontenidas.")
            _DORNHEIM_FALLBACK_WARNED = True
        parsed = _build_fallback_dornheim_paths(page_num)

    _DORNHEIM_CACHE[cache_key] = parsed
    return parsed


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


def draw_annot_box(ax, x, y, text, color, fontsize=11.5, weight='bold', ha='left', va='center', xycoords='data'):
    """Dibuja una anotación con fondo blanco redondeado para máxima legibilidad."""
    transform = ax.transAxes if xycoords == 'axes fraction' else ax.transData
    return ax.text(x, y, text, color=color, fontsize=fontsize, weight=weight, ha=ha, va=va,
                   transform=transform,
                   bbox=dict(boxstyle="round,pad=0.24", facecolor="white", edgecolor=color, alpha=0.92, lw=0.85),
                   zorder=10)


# ==============================================================================
# FIGURAS DISCRETAS INDIVIDUALES (PÁGINAS 22, 24 Y 26)
# ==============================================================================

def make_dornheim_fig7_open_close(pdf_path=None):
    paths22, _ = load_dornheim_paths(22, pdf_path=pdf_path)

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
    return save_figure(fig, "dornheim_fig7_open_close.png")


def make_dornheim_fig8_advance_recede_bisection(pdf_path=None):
    paths24, _ = load_dornheim_paths(24, pdf_path=pdf_path)
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
    return save_figure(fig, "dornheim_fig8_advance_recede_bisection.png")


def make_dornheim_fig9_swap(pdf_path=None):
    paths24, _ = load_dornheim_paths(24, pdf_path=pdf_path)
    paths26, _ = load_dornheim_paths(26, pdf_path=pdf_path)
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
    return save_figure(fig, "dornheim_fig9_swap.png")


# ==============================================================================
# PANEL COMPLETO 2x2: TODOS LOS MOVIMIENTOS WORM
# ==============================================================================

def make_dornheim_panel_completo(pdf_path=None):
    paths22, _ = load_dornheim_paths(22, pdf_path=pdf_path)
    paths24, _ = load_dornheim_paths(24, pdf_path=pdf_path)
    paths26, _ = load_dornheim_paths(26, pdf_path=pdf_path)
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

    # Resto de la cadena definida de Partícula 2 (tau 0..2 y 12..20)
    rest_1 = p2_full[0:3]   # tau 0..2
    rest_2 = p2_full[12:21] # tau 12..20
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
    draw_annot_box(ax, -2.4, 18.2, "Partícula 1\n(Cerrada)", COLOR_P1, fontsize=11, ha="center")
    draw_annot_box(ax, 2.2, 18.2, "Resto de cadena\n(Partícula 2)", COLOR_P2_REST, fontsize=10.5, ha="center")
    draw_annot_box(ax, ira_new[0] + 0.35, ira_new[1], r"$\mathcal{I}_{\mathrm{new}}$ (Advance $\uparrow$)", COLOR_MODIFIED, fontsize=10.5)
    draw_annot_box(ax, ira_old[0] + 0.35, ira_old[1] - 0.35, r"$\mathcal{I}_{\mathrm{old}}$ (Recede $\downarrow$)", COLOR_IRA, fontsize=10.5)
    draw_annot_box(ax, masha[0] + 0.35, masha[1], r"$\mathcal{M}$ (Masha)", COLOR_MASHA, fontsize=10.5)
    draw_annot_box(ax, 0.2, 10.5, "Cuentas añadidas\n(Advance azul)", COLOR_MODIFIED, fontsize=10.0, ha="center")
    draw_annot_box(ax, 0.25, 4.5, "Partícula 2\n(Gusano activo)", COLOR_P2, fontsize=10.0, ha="center")

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

    ira_bis = (xs_bis[-1], ts_bis[-1])
    masha_bis = (xs_bis[0], ts_bis[0])
    ax.scatter([ira_bis[0]], [ira_bis[1]], color=COLOR_IRA, edgecolor='black', s=125, zorder=6, marker='o')
    ax.scatter([masha_bis[0]], [masha_bis[1]], color=COLOR_MASHA, edgecolor='black', s=125, zorder=6, marker='s')

    ax.set_title(r"$\mathbf{(c)\ Bisección\ /\ IraSpatialMove}$", fontsize=14, fontweight="bold", pad=10)
    draw_annot_box(ax, -2.4, 18.2, "Partícula 1", COLOR_P1, fontsize=11, ha="center")
    draw_annot_box(ax, 2.2, 18.2, "Partícula 2", COLOR_P2, fontsize=11, ha="center")
    draw_annot_box(ax, 2.25, 9.5, "Segmento interno\nremuestreado (azul)", COLOR_MODIFIED, fontsize=10.5, ha="center")
    draw_annot_box(ax, ira_bis[0] + 0.35, ira_bis[1] + 0.35, r"$\mathcal{I}$ (Ira)", COLOR_IRA, fontsize=10.5)
    draw_annot_box(ax, masha_bis[0] - 0.35, masha_bis[1] - 0.35, r"$\mathcal{M}$ (Masha)", COLOR_MASHA, fontsize=10.5, ha="right")

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

    # Conserva visible el tramo restante de P2 abierta con un color bien distinguible.
    color_p2_open_rest = "#4b5563"
    color_p2_open_rest_bead = "#9ca3af"
    p2_legacy = p2_full[12:21]
    xs_legacy, ts_legacy = zip(*p2_legacy)
    ax.plot(xs_legacy, ts_legacy, color=color_p2_open_rest, lw=2.2, zorder=2)
    ax.scatter(xs_legacy, ts_legacy, color=color_p2_open_rest_bead, edgecolor=color_p2_open_rest, s=42, zorder=3)

    ax.plot(xs_sw, ts_sw, color=COLOR_SWAP, lw=3.4, zorder=5)
    ax.scatter(xs_sw, ts_sw, color=COLOR_SWAP_BEAD, edgecolor=COLOR_SWAP, s=62, zorder=6)

    ira_new_sw = (xs_p1_a[-1], ts_p1_a[-1])
    ira_old_sw = (xs_w2[-1], ts_w2[-1])
    masha_sw = (xs_w2[0], ts_w2[0])
    p1_open_end = (xs_p1_b[0], ts_p1_b[0])

    ax.scatter([ira_new_sw[0]], [ira_new_sw[1]], color=COLOR_IRA, edgecolor='black', s=140, zorder=7, marker='o')
    ax.scatter([ira_old_sw[0]], [ira_old_sw[1]], color=COLOR_SWAP, edgecolor='black', s=120, zorder=7, marker='o')
    ax.scatter([masha_sw[0]], [masha_sw[1]], color=COLOR_MASHA, edgecolor='black', s=130, zorder=7, marker='s')
    ax.scatter([p1_open_end[0]], [p1_open_end[1]], facecolor='white', edgecolor=COLOR_P1, linewidth=2, s=100, zorder=8, marker='o')

    ax.set_title(r"$\mathbf{(d)\ Movimiento\ Swap}$ (Permutación Cuántica)", fontsize=14, fontweight="bold", pad=10)
    draw_annot_box(ax, 0.3, 14.5, r"Enlace Swap $\longrightarrow$", COLOR_SWAP, fontsize=12, ha="center")
    draw_annot_box(ax, -2.35, 14.1, "Cadena P1 abierta\n(enlace cortado)", COLOR_IRA, fontsize=9.5, ha="center")
    draw_annot_box(ax, 2.45, 18.1, "Cadena P2 abierta\n(tramo restante)", color_p2_open_rest, fontsize=9.5, ha="center")
    draw_annot_box(ax, ira_new_sw[0] - 0.35, ira_new_sw[1] + 0.85, r"$\mathcal{I}_{\mathrm{new}}$ (Nueva Ira)", COLOR_IRA, fontsize=10.5, ha="right")
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


PARTICLE_PALETTE = [
    (COLOR_P1, COLOR_P1_BEAD),
    (COLOR_P2, COLOR_P2_BEAD),
    ("#8c564b", "#c49c94"),
    ("#e377c2", "#f7b6d2"),
    ("#7f7f7f", "#c7c7c7"),
    ("#bcbd22", "#dbdb8d"),
    ("#17becf", "#9edae5"),
]


def _particle_colors(particle, p0_color=COLOR_P1, p0_bead=COLOR_P1_BEAD,
                     p1_color=COLOR_P2, p1_bead=COLOR_P2_BEAD):
    if particle == 0:
        return p0_color, p0_bead
    if particle == 1:
        return p1_color, p1_bead
    return PARTICLE_PALETTE[particle % len(PARTICLE_PALETTE)]


def _add_segments(ax, segments, **kwargs):
    if segments:
        ax.add_collection(LineCollection(segments, **kwargs))


def draw_simulation_state(ax, state, title):
    n_beads = state["n_beads"]
    n_particles = state["n_particles"]
    tau = np.arange(n_beads + 1) / n_beads
    open_edges = open_worm_edges(state)

    open_segments = []
    closed_by_particle = [[] for _ in range(n_particles)]
    for bead in range(n_beads):
        next_bead = (bead + 1) % n_beads
        t0, t1 = tau[bead], tau[bead + 1]
        for particle in range(n_particles):
            next_particle = int(state["next"][bead, particle])
            if next_particle < 0:
                continue
            seg = [
                (state["x"][bead, particle], t0),
                (state["x"][next_bead, next_particle], t1),
            ]
            if (bead, particle) in open_edges:
                open_segments.append(seg)
            else:
                closed_by_particle[particle].append(seg)

    for particle, segs in enumerate(closed_by_particle):
        c_line, _ = _particle_colors(particle)
        _add_segments(ax, segs, colors=c_line, linewidths=1.5, zorder=2)
    _add_segments(ax, open_segments, colors=COLOR_MODIFIED, linewidths=2.4, zorder=3)

    bead_times = np.arange(n_beads) / n_beads
    for particle in range(n_particles):
        c_line, c_bead = _particle_colors(particle)
        ax.scatter(state["x"][:, particle], bead_times,
                   color=c_bead, edgecolor=c_line, s=18, zorder=4)

    xmin = float(np.min(state["x"]))
    xmax = float(np.max(state["x"]))
    span = max(xmax - xmin, 0.5)
    offset = max(0.06, 0.04 * span)

    if state["worm_open"]:
        ira_bead, ira_particle = state["ira"]
        masha_bead, masha_particle = state["masha"]
        ax.scatter([state["x"][ira_bead, ira_particle]], [bead_times[ira_bead]],
                   color=COLOR_IRA, edgecolor="black", linewidth=1.0, s=70, zorder=5)
        draw_annot_box(ax, state["x"][ira_bead, ira_particle] + offset, bead_times[ira_bead],
                       r"$\mathcal{I}$ (Ira)", COLOR_IRA, fontsize=10.5)

        ax.scatter([state["x"][masha_bead, masha_particle]], [bead_times[masha_bead]],
                   color=COLOR_MASHA, edgecolor="black", linewidth=1.0, s=70, zorder=5, marker="s")
        draw_annot_box(ax, state["x"][masha_bead, masha_particle] - offset, bead_times[masha_bead],
                       r"$\mathcal{M}$ (Masha)", COLOR_MASHA, fontsize=10.5, ha="right")

    padding = max(0.35, 0.12 * span)
    ax.set_xlim(xmin - padding, xmax + padding)
    ax.set_ylim(0, 1)
    ax.set_xlabel("x (coordenada de la simulación)", fontsize=12, fontweight="bold")
    ax.set_ylabel(r"$\tau/\beta$", fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.45, color="#b0b0b0")
    ax.set_axisbelow(True)
    ax.tick_params(direction="in", top=True, right=True, labelsize=10.5)


def _draw_sim_panel_worldlines(ax, state, p0_color=COLOR_P1, p0_bead=COLOR_P1_BEAD,
                               p1_color=COLOR_P2, p1_bead=COLOR_P2_BEAD,
                               rest_edges=None, skip_edges=None):
    """Dibuja las trayectorias reales de la simulación diferenciando partículas y tramos."""
    n_beads = state["n_beads"]
    n_particles = state["n_particles"]
    tau = np.arange(n_beads + 1) / n_beads
    rest_edges = rest_edges or set()
    skip_edges = skip_edges or set()

    rest_segments = []
    particle_segments = [[] for _ in range(n_particles)]
    for bead in range(n_beads):
        next_bead = (bead + 1) % n_beads
        t0, t1 = tau[bead], tau[bead + 1]
        for particle in range(n_particles):
            next_particle = int(state["next"][bead, particle])
            if next_particle < 0 or (bead, particle) in skip_edges:
                continue
            seg = [
                (state["x"][bead, particle], t0),
                (state["x"][next_bead, next_particle], t1),
            ]
            if (bead, particle) in rest_edges:
                rest_segments.append(seg)
            else:
                particle_segments[particle].append(seg)

    _add_segments(ax, rest_segments, colors="#4b5563", linestyles=":", linewidths=2.0, zorder=2)
    for particle, segs in enumerate(particle_segments):
        c_edge, _ = _particle_colors(particle, p0_color=p0_color, p0_bead=p0_bead,
                                     p1_color=p1_color, p1_bead=p1_bead)
        _add_segments(ax, segs, colors=c_edge, linestyles="-", linewidths=2.0, zorder=2)

    bead_times = np.arange(n_beads) / n_beads
    for particle in range(n_particles):
        c_edge, c_bead = _particle_colors(particle, p0_color=p0_color, p0_bead=p0_bead,
                                          p1_color=p1_color, p1_bead=p1_bead)
        ax.scatter(state["x"][:, particle], bead_times, color=c_bead, edgecolor=c_edge, s=26, zorder=3)


def _setup_sim_panel_axis(ax, state, title):
    xmin = float(np.min(state["x"]))
    xmax = float(np.max(state["x"]))
    span = max(xmax - xmin, 0.5)
    ax.set_xlim(xmin - 0.38 * span, xmax + 0.42 * span)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("x (simulación)", fontsize=12, fontweight="bold")
    ax.set_ylabel(r"$\tau/\beta$", fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=13.5, fontweight="bold", pad=10)
    ax.grid(True, linestyle=":", alpha=0.45, color="#b0b0b0")
    ax.set_axisbelow(True)
    ax.tick_params(direction="in", length=5, width=1.0, top=True, right=True, labelsize=10.5)
    return xmin, xmax, span


def make_simulation_panel_completo(snapshot_file):
    """Genera figura_panel_completo_worm (2x2) a partir de los snapshots reales de simulación."""
    events = read_worm_snapshots(snapshot_file)
    required = ["Open", "Advance", "IraSpatialMove", "Swap"]
    missing = [m for m in required if m not in events or "before" not in events[m] or "after" not in events[m]]
    if missing:
        raise ValueError(f"Faltan eventos completos en {snapshot_file}: {missing}")

    fig, axs = plt.subplots(2, 2, figsize=(13.5, 11.5), dpi=300)

    # --------------------------------------------------------------------------
    # (0, 0): Open / Close
    # --------------------------------------------------------------------------
    ax = axs[0, 0]
    b_open = events["Open"]["before"]
    a_open = events["Open"]["after"]
    n_beads = a_open["n_beads"]
    tau = np.arange(n_beads + 1) / n_beads
    _, _, span = _setup_sim_panel_axis(ax, a_open, r"$\mathbf{(a)\ Movimiento\ Open\ /\ Close}$ ($Z \leftrightarrow G$)")
    _draw_sim_panel_worldlines(ax, a_open)

    # Enlace cortado en Open (presente en before, -1 en after)
    cut_edges = [
        (m, p, int(b_open["next"][m, p]))
        for m in range(n_beads)
        for p in range(a_open["n_particles"])
        if b_open["next"][m, p] >= 0 and a_open["next"][m, p] < 0
    ]
    for m, p, np_b in cut_edges:
        nb = (m + 1) % n_beads
        ax.plot([b_open["x"][m, p], b_open["x"][nb, np_b]], [tau[m], tau[m + 1]],
                color=COLOR_MODIFIED, lw=3.0, linestyle="--", zorder=4)

    ira_b, ira_p = a_open["ira"]
    masha_b, masha_p = a_open["masha"]
    ix, itau = a_open["x"][ira_b, ira_p], tau[ira_b]
    mx, mtau = a_open["x"][masha_b, masha_p], tau[masha_b]
    ax.scatter([ix], [itau], color=COLOR_IRA, edgecolor="black", s=120, zorder=6, marker="o")
    ax.scatter([mx], [mtau], color=COLOR_MASHA, edgecolor="black", s=120, zorder=6, marker="s")

    dx = 0.05 * span
    draw_annot_box(ax, 0.03, 0.94, "Partícula 1\n(Sector G)", COLOR_P1, fontsize=10,
                   ha="left", va="top", xycoords="axes fraction")
    draw_annot_box(ax, 0.97, 0.94, "Partícula 2\n(Cerrada)", COLOR_P2, fontsize=10,
                   ha="right", va="top", xycoords="axes fraction")
    draw_annot_box(ax, ix - dx, float(np.clip(itau - 0.07, 0.06, 0.94)),
                   r"$\mathcal{I}$ (Ira)", COLOR_IRA, fontsize=10.5, ha="right")
    draw_annot_box(ax, mx + dx, float(np.clip(mtau + 0.05, 0.06, 0.94)),
                   r"$\mathcal{M}$ (Masha)", COLOR_MASHA, fontsize=10.5, ha="left")
    draw_annot_box(ax, 0.97, 0.07, "Enlace cortado (Open)\no cerrado (Close)", COLOR_MODIFIED,
                   fontsize=9.5, ha="right", va="bottom", xycoords="axes fraction")

    # --------------------------------------------------------------------------
    # (0, 1): Advance / Recede
    # --------------------------------------------------------------------------
    ax = axs[0, 1]
    b_adv = events["Advance"]["before"]
    a_adv = events["Advance"]["after"]
    n_beads = a_adv["n_beads"]
    tau = np.arange(n_beads + 1) / n_beads
    _, _, span = _setup_sim_panel_axis(ax, a_adv, r"$\mathbf{(b)\ Movimiento\ Advance\ /\ Recede}$ (Sector $G$)")

    open_after = open_worm_edges(a_adv)
    all_edges_adv = {
        (m, p) for m in range(n_beads) for p in range(a_adv["n_particles"]) if a_adv["next"][m, p] >= 0
    }
    rest_adv = all_edges_adv - open_after
    added_edges = [
        (m, p, int(a_adv["next"][m, p]))
        for m in range(n_beads)
        for p in range(a_adv["n_particles"])
        if b_adv["next"][m, p] < 0 and a_adv["next"][m, p] >= 0
    ]
    added_set = {(m, p) for m, p, _ in added_edges}

    _draw_sim_panel_worldlines(ax, a_adv, rest_edges=rest_adv, skip_edges=added_set)
    for m, p, np_a in added_edges:
        nb = (m + 1) % n_beads
        ax.plot([a_adv["x"][m, p], a_adv["x"][nb, np_a]], [tau[m], tau[m + 1]],
                color=COLOR_MODIFIED, lw=3.4, zorder=5)
        ax.scatter([a_adv["x"][nb, np_a]], [tau[nb]], color=COLOR_MODIFIED_BEAD, edgecolor=COLOR_MODIFIED, s=55, zorder=5)

    ib_b, ib_p = b_adv["ira"]
    ia_b, ia_p = a_adv["ira"]
    ma_b, ma_p = a_adv["masha"]
    ax.scatter([a_adv["x"][ia_b, ia_p]], [tau[ia_b]], color=COLOR_MODIFIED, edgecolor="black", s=130, zorder=6, marker="o")
    ax.scatter([b_adv["x"][ib_b, ib_p]], [tau[ib_b]], color=COLOR_IRA, edgecolor="black", s=115, zorder=6, marker="o")
    ax.scatter([a_adv["x"][ma_b, ma_p]], [tau[ma_b]], color=COLOR_MASHA, edgecolor="black", s=120, zorder=6, marker="s")

    dx = 0.05 * span
    draw_annot_box(ax, 0.03, 0.94, "Partícula 1", COLOR_P1, fontsize=10,
                   ha="left", va="top", xycoords="axes fraction")
    draw_annot_box(ax, 0.97, 0.94, "Partícula 2\n(Gusano activo)", COLOR_P2, fontsize=10,
                   ha="right", va="top", xycoords="axes fraction")
    draw_annot_box(ax, a_adv["x"][ia_b, ia_p] + dx, float(np.clip(tau[ia_b] + 0.04, 0.06, 0.94)),
                   r"$\mathcal{I}_{\mathrm{new}}$ (Advance $\uparrow$)", COLOR_MODIFIED, fontsize=10)
    draw_annot_box(ax, b_adv["x"][ib_b, ib_p] + dx, float(np.clip(tau[ib_b] - 0.04, 0.06, 0.94)),
                   r"$\mathcal{I}_{\mathrm{old}}$ (Recede $\downarrow$)", COLOR_IRA, fontsize=10)
    draw_annot_box(ax, a_adv["x"][ma_b, ma_p] + dx, float(np.clip(tau[ma_b] + 0.05, 0.06, 0.94)),
                   r"$\mathcal{M}$ (Masha)", COLOR_MASHA, fontsize=10)
    draw_annot_box(ax, 0.03, 0.07, "Enlace añadido\n(Advance azul)", COLOR_MODIFIED, fontsize=9.5,
                   ha="left", va="bottom", xycoords="axes fraction")

    # --------------------------------------------------------------------------
    # (1, 0): Bisección / IraSpatialMove
    # --------------------------------------------------------------------------
    ax = axs[1, 0]
    b_ira = events["IraSpatialMove"]["before"]
    a_ira = events["IraSpatialMove"]["after"]
    n_beads = a_ira["n_beads"]
    tau = np.arange(n_beads + 1) / n_beads
    _, _, span = _setup_sim_panel_axis(ax, a_ira, r"$\mathbf{(c)\ Bisección\ /\ IraSpatialMove}$")
    _draw_sim_panel_worldlines(ax, a_ira)

    # Resaltar enlace incidente a la cuenta desplazada (en x o y)
    moved = [
        (m, p) for m in range(n_beads) for p in range(a_ira["n_particles"])
        if np.hypot(b_ira["x"][m, p] - a_ira["x"][m, p], b_ira["y"][m, p] - a_ira["y"][m, p]) > 1e-9
    ]
    for m, p in moved:
        prev_b = (m - 1) % n_beads
        for src_p in range(a_ira["n_particles"]):
            if a_ira["next"][prev_b, src_p] == p:
                ax.plot([b_ira["x"][prev_b, src_p], b_ira["x"][m, p]], [tau[prev_b], tau[m]],
                        color="#9ca3af", lw=2.2, linestyle=":", zorder=4)
                ax.plot([a_ira["x"][prev_b, src_p], a_ira["x"][m, p]], [tau[prev_b], tau[m]],
                        color=COLOR_MODIFIED, lw=3.0, linestyle="--", zorder=5)

    ix, itau = a_ira["x"][a_ira["ira"]], tau[a_ira["ira"][0]]
    mx, mtau = a_ira["x"][a_ira["masha"]], tau[a_ira["masha"][0]]
    ax.scatter([ix], [itau], color=COLOR_IRA, edgecolor="black", s=120, zorder=6, marker="o")
    ax.scatter([mx], [mtau], color=COLOR_MASHA, edgecolor="black", s=120, zorder=6, marker="s")

    dx = 0.05 * span
    draw_annot_box(ax, 0.03, 0.94, "Partícula 1", COLOR_P1, fontsize=10,
                   ha="left", va="top", xycoords="axes fraction")
    draw_annot_box(ax, 0.97, 0.94, "Partícula 2", COLOR_P2, fontsize=10,
                   ha="right", va="top", xycoords="axes fraction")
    draw_annot_box(ax, ix - dx, float(np.clip(itau - 0.07, 0.06, 0.94)),
                   r"$\mathcal{I}$ (Ira movida)", COLOR_IRA, fontsize=10, ha="right")
    draw_annot_box(ax, mx + dx, float(np.clip(mtau + 0.05, 0.06, 0.94)),
                   r"$\mathcal{M}$ (Masha)", COLOR_MASHA, fontsize=10)
    draw_annot_box(ax, 0.97, 0.07, "Enlace remuestreado\n(IraSpatialMove azul)", COLOR_MODIFIED,
                   fontsize=9.5, ha="right", va="bottom", xycoords="axes fraction")

    # --------------------------------------------------------------------------
    # (1, 1): Swap
    # --------------------------------------------------------------------------
    ax = axs[1, 1]
    b_sw = events["Swap"]["before"]
    a_sw = events["Swap"]["after"]
    n_beads = a_sw["n_beads"]
    tau = np.arange(n_beads + 1) / n_beads
    _, _, span = _setup_sim_panel_axis(ax, a_sw, r"$\mathbf{(d)\ Movimiento\ Swap}$ (Permutación Cuántica)")

    swap_created = [
        (m, p, int(a_sw["next"][m, p]))
        for m in range(n_beads)
        for p in range(a_sw["n_particles"])
        if b_sw["next"][m, p] < 0 and a_sw["next"][m, p] >= 0
    ]
    swap_cut = [
        (m, p, int(b_sw["next"][m, p]))
        for m in range(n_beads)
        for p in range(a_sw["n_particles"])
        if b_sw["next"][m, p] >= 0 and a_sw["next"][m, p] < 0
    ]
    swap_set = {(m, p) for m, p, _ in swap_created}

    # Tramo restante de cadena abierta siguiendo next periódicamente tras el corte
    p2_rest_edges = set()
    if swap_cut:
        curr_b = (swap_cut[0][0] + 1) % n_beads
        curr_p = swap_cut[0][2]
        visited_rest = set()
        for _ in range(n_beads * a_sw["n_particles"]):
            state_key = (curr_b, curr_p)
            if state_key in visited_rest or state_key == a_sw["masha"]:
                break
            visited_rest.add(state_key)
            next_p = int(a_sw["next"][curr_b, curr_p])
            if next_p < 0:
                break
            p2_rest_edges.add(state_key)
            curr_b = (curr_b + 1) % n_beads
            curr_p = next_p

    _draw_sim_panel_worldlines(ax, a_sw, rest_edges=p2_rest_edges, skip_edges=swap_set)

    # Dibujar cuentas del tramo restante con color gris pizarra para máxima distinción
    color_p2_open_rest = "#4b5563"
    color_p2_open_rest_bead = "#9ca3af"
    if p2_rest_edges:
        rest_nodes = sorted(
            set(p2_rest_edges)
            | {((m + 1) % n_beads, int(a_sw["next"][m, p])) for m, p in p2_rest_edges}
        )
        rest_x = [a_sw["x"][m, p] for m, p in rest_nodes]
        rest_t = [tau[m] for m, _ in rest_nodes]
        ax.scatter(rest_x, rest_t, color=color_p2_open_rest_bead, edgecolor=color_p2_open_rest, s=28, zorder=4)

    for m, p, np_a in swap_created:
        nb = (m + 1) % n_beads
        ax.plot([a_sw["x"][m, p], a_sw["x"][nb, np_a]], [tau[m], tau[m + 1]],
                color=COLOR_SWAP, lw=3.4, zorder=6)
    for m, p, np_b in swap_cut:
        nb = (m + 1) % n_beads
        ax.plot([b_sw["x"][m, p], b_sw["x"][nb, np_b]], [tau[m], tau[m + 1]],
                color=COLOR_IRA, lw=2.2, linestyle="--", zorder=4)
        ax.scatter([a_sw["x"][nb, np_b]], [tau[nb]], facecolor="white", edgecolor=COLOR_P1, linewidth=2, s=90, zorder=7)

    ib_sw = b_sw["ira"]
    ia_sw = a_sw["ira"]
    ma_sw = a_sw["masha"]
    ax.scatter([a_sw["x"][ia_sw]], [tau[ia_sw[0]]], color=COLOR_IRA, edgecolor="black", s=130, zorder=7, marker="o")
    ax.scatter([b_sw["x"][ib_sw]], [tau[ib_sw[0]]], color=COLOR_SWAP, edgecolor="black", s=115, zorder=7, marker="o")
    ax.scatter([a_sw["x"][ma_sw]], [tau[ma_sw[0]]], color=COLOR_MASHA, edgecolor="black", s=120, zorder=7, marker="s")

    dx = 0.05 * span
    draw_annot_box(ax, 0.03, 0.94, "Partícula 1", COLOR_P1, fontsize=10,
                   ha="left", va="top", xycoords="axes fraction")
    draw_annot_box(ax, 0.97, 0.94, "Cadena P2 abierta\n(tramo restante)", color_p2_open_rest,
                   fontsize=9.5, ha="right", va="top", xycoords="axes fraction")
    draw_annot_box(ax, 0.03, 0.36, "Enlace Swap\n(P1 → P2)", COLOR_SWAP, fontsize=9.5,
                   ha="left", va="center", xycoords="axes fraction")
    draw_annot_box(ax, a_sw["x"][ia_sw] - dx, float(np.clip(tau[ia_sw[0]] - 0.07, 0.06, 0.94)),
                   r"$\mathcal{I}_{\mathrm{new}}$ (Nueva Ira)", COLOR_IRA, fontsize=9.5, ha="right")
    draw_annot_box(ax, b_sw["x"][ib_sw] + dx, float(np.clip(tau[ib_sw[0]] - 0.06, 0.06, 0.94)),
                   r"$\mathcal{I}_{\mathrm{old}}$ (Vieja Ira)", COLOR_SWAP, fontsize=9.5, ha="left")
    draw_annot_box(ax, a_sw["x"][ma_sw] + dx, float(np.clip(tau[ma_sw[0]] + 0.05, 0.06, 0.94)),
                   r"$\mathcal{M}$ (Masha)", COLOR_MASHA, fontsize=10)

    fig.suptitle("Transiciones reales del algoritmo Worm (outworm_snapshots.dat)", fontweight="bold", fontsize=15, y=0.99)
    plt.tight_layout(rect=(0, 0, 1, 0.97))
    save_figure(fig, "figura_panel_completo_worm.png")


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
    parser.add_argument(
        "--dornheim-pdf",
        type=Path,
        default=None,
        metavar="PDF",
        help="Ruta opcional a PhD Dornheim 2014.pdf (o usar variable DORNHEIM_PDF)",
    )
    args = parser.parse_args()

    if args.from_simulation is not None:
        make_simulation_figures(args.from_simulation)
        make_simulation_panel_completo(args.from_simulation)
        return

    print("=" * 50)
    print("Generando gráficos movimiento partículas Worm")
    print("=" * 50)

    # Limpiar figuras obsoletas del antiguo estilo didáctico
    clean_obsolete_didactic_files()

    print("\nGenerando gráficos.")
    make_dornheim_fig7_open_close(pdf_path=args.dornheim_pdf)
    make_dornheim_fig8_advance_recede_bisection(pdf_path=args.dornheim_pdf)
    make_dornheim_fig9_swap(pdf_path=args.dornheim_pdf)

    default_snapshots = SCRIPT_DIR / "outworm_snapshots.dat"
    copied_panel = False
    if default_snapshots.exists():
        try:
            make_simulation_panel_completo(default_snapshots)
            copied_panel = True
        except ValueError as exc:
            print(f"  [AVISO] {exc}. Generando panel completo desde esquemas de Dornheim...")

    if not copied_panel:
        make_dornheim_panel_completo(pdf_path=args.dornheim_pdf)

    print("\n" + "=" * 76)
    print("¡Éxito! Todas las figuras han sido generadas y reemplazadas.")
    print("Ubicaciones configuradas (PNG 300 DPI y PDF):")
    print(f"  1. {OUTPUT_DIR}")
    print(f"  2. {OVERLEAF_DIR} (si es accesible)")
    print("=" * 76)


if __name__ == "__main__":
    main()
