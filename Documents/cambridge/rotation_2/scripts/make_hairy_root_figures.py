#!/usr/bin/env python3
"""
make_hairy_root_figures.py

Generates BF/FL composite figures for hairy root transformation image sets
(lines 43224 LSH2oe and 43225 LSH1oe).

Layout: BF top row / FL bottom row, all plants per line,
        FL images contrast-enhanced to reveal DsRed signal.

Usage:
    python3 scripts/make_hairy_root_figures.py

Outputs:
    figures/hairy_root_lsh2oe_v2.png
    figures/hairy_root_lsh1oe_v2.png
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageOps

# 300 px = 6 mm  →  50 px/mm
PX_PER_MM = 50

REPO = Path(__file__).resolve().parent.parent
IMG_DIR = REPO / "source_images"
FIG_DIR = REPO / "figures"
FIG_DIR.mkdir(exist_ok=True)

LINES = {
    "LSH2oe": {
        "label": "pLjUBI::DsRed-LSH2",
        "bf": [IMG_DIR / f"43224_LSH2oe_bf_{i}.jpeg" for i in range(1, 4)],
        "fl": [IMG_DIR / f"43224_LSH2oe_fl_{i}.jpeg" for i in range(1, 4)],
    },
    "LSH1oe": {
        "label": "pLjUBI::DsRed-LSH1",
        "bf": [IMG_DIR / f"43225_LSH1oe_bf_{i}.jpeg" for i in (2, 3, 4)],
        "fl": [IMG_DIR / f"43225_LSH1oe_fl_{i}.jpeg" for i in (2, 3, 4)],
    },
}


def load_img(path, enhance=False):
    img = Image.open(path).convert("RGB")
    if enhance:
        img = ImageOps.autocontrast(img, cutoff=0.5)
    return np.array(img)


def add_scale_bar(ax, img_shape, bar_mm=5, label="5 mm"):
    """Draw a white scale bar in the bottom-right corner of ax (image coords)."""
    h, w = img_shape[:2]
    bar_px  = bar_mm * PX_PER_MM
    margin_x = max(15, int(w * 0.03))
    margin_y = max(15, int(h * 0.03))
    bar_h    = max(6, h // 70)

    x1 = w - margin_x - bar_px
    y2 = h - margin_y
    y1 = y2 - bar_h

    ax.add_patch(mpatches.Rectangle(
        (x1, y1), bar_px, bar_h,
        color="white", zorder=5, clip_on=True,
    ))
    ax.text(
        x1 + bar_px / 2, y1 - 4, label,
        color="white", ha="center", va="bottom",
        fontsize=8, fontweight="bold", zorder=6,
    )


def panel_letter(ax, letter, fontsize=13):
    ax.text(
        0.025, 0.97, letter,
        transform=ax.transAxes,
        fontsize=fontsize, fontweight="bold", color="white",
        va="top", ha="left",
        bbox=dict(facecolor="black", alpha=0.55, pad=1.5, edgecolor="none"),
    )


def clean_ax(ax):
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])


def row_label(ax, text, fontsize=11):
    ax.set_ylabel(text, fontsize=fontsize, labelpad=5, rotation=90,
                  va="center", ha="center")
    ax.yaxis.set_visible(True)
    ax.set_yticks([])


def make_figure(line_key):
    info = LINES[line_key]
    bf, fl = info["bf"], info["fl"]
    n = len(bf)

    fig, axes = plt.subplots(2, n, figsize=(n * 4.8, 2 * 3.8), facecolor="white")
    if n == 1:
        axes = axes.reshape(2, 1)
    fig.subplots_adjust(wspace=0.025, hspace=0.05, left=0.07,
                        right=0.995, top=0.93, bottom=0.015)

    alpha_seq = "ABCDEFGH"
    for col in range(n):
        bf_arr = load_img(bf[col], enhance=False)
        fl_arr = load_img(fl[col], enhance=True)
        axes[0, col].imshow(bf_arr)
        axes[1, col].imshow(fl_arr)
        panel_letter(axes[0, col], alpha_seq[col])
        panel_letter(axes[1, col], alpha_seq[n + col])
        add_scale_bar(axes[0, col], bf_arr.shape)
        add_scale_bar(axes[1, col], fl_arr.shape)
        clean_ax(axes[0, col])
        clean_ax(axes[1, col])

    row_label(axes[0, 0], "Brightfield")
    row_label(axes[1, 0], "Fluorescence")
    fig.suptitle(info["label"], fontsize=13, fontweight="bold", y=0.995)
    return fig


if __name__ == "__main__":
    for line_key in ("LSH2oe", "LSH1oe"):
        fig = make_figure(line_key)
        out = FIG_DIR / f"hairy_root_{line_key.lower()}_v2.png"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved  {out.relative_to(REPO)}")
