#!/usr/bin/env python3
"""
make_gus_figure.py

Generates figures/gus_comparison.png — 3 × 2 composite comparing GUS staining
protocols in P. sativum nodules.

Layout:
    Row 1 (8 h / 9 min):   mature_nodules_batch1      early_nodules_batch1
    Row 2 (18 h / 15 min): infection_thread_detail_batch2  mature_nodules_batch2
    Row 3 (12 h / 15 min): infection_thread_detail_batch3  mature_nodules_batch3

Usage:
    python3 scripts/make_gus_figure.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

REPO = Path(__file__).resolve().parent.parent
FIG_DIR = REPO / "figures"
FIG_DIR.mkdir(exist_ok=True)

ROWS = [
    {
        "label": "8 h / 9 min",
        "images": [
            FIG_DIR / "mature_nodules_batch1.jpg",
            FIG_DIR / "early_nodules_batch1.jpg",
        ],
    },
    {
        "label": "18 h / 15 min",
        "images": [
            FIG_DIR / "infection_thread_detail_batch2.jpg",
            FIG_DIR / "mature_nodules_batch2.jpg",
        ],
    },
    {
        "label": "12 h / 15 min",
        "images": [
            FIG_DIR / "infection_thread_detail_batch3.jpg",
            FIG_DIR / "mature_nodules_batch3.jpg",
        ],
    },
]


def load_img(path):
    return np.array(Image.open(path).convert("RGB"))


def panel_letter(ax, letter, fontsize=18):
    ax.text(
        0.025, 0.97, letter,
        transform=ax.transAxes,
        fontsize=fontsize, fontweight="bold", color="black",
        va="top", ha="left",
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


if __name__ == "__main__":
    n_rows = len(ROWS)
    n_cols = 2

    aspect = 4 / 3  # all GUS images are 1600x1200
    cell_w = 4.8
    cell_h = cell_w / aspect
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(n_cols * cell_w, n_rows * cell_h),
                             facecolor="white")
    fig.subplots_adjust(wspace=0.025, hspace=0.0,
                        left=0.1, right=0.995, top=0.995, bottom=0.015)

    alpha_seq = "ABCDEF"
    for r, row in enumerate(ROWS):
        for c, path in enumerate(row["images"]):
            ax = axes[r, c]
            ax.imshow(load_img(path))
            panel_letter(ax, alpha_seq[r * n_cols + c])
            clean_ax(ax)
        row_label(axes[r, 0], row["label"])

    out = FIG_DIR / "gus_comparison.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved  {out.relative_to(REPO)}")
