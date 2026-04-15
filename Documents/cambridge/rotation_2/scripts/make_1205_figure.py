#!/usr/bin/env python3
"""
make_1205_figure.py

Composite figure for EMS accession 1205-2, showing hypernodulation phenotype.

Layout:
  A (left, full height)  — overview of 1205-2(1) root system (JIC M2 phenotyping)
  B (top right)          — Keyence close-up, sibling 1: dense clustered nodules
  C (middle right)       — Keyence close-up, sibling 2: scattered small nodules
  D (bottom right)       — Wild-type JI2822 nodule morphology for comparison

Output: figures/ems_1205.png
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from PIL import Image

REPO     = Path(__file__).resolve().parent.parent
IMG_DIR  = REPO / "source_images"
OUT      = REPO / "figures" / "ems_1205.png"

OVERVIEW = IMG_DIR / "1205-2_1_close.jpg"
SIB1     = IMG_DIR / "1205-2-1_keyence.jpg"
SIB2     = IMG_DIR / "1205-2-3_keyence.jpg"
WT       = IMG_DIR / "cool random keyence images" / "2026_03_10_JI-5-2000005.jpg"


def load(path, rotate=0):
    img = Image.open(path).convert("RGB")
    if rotate:
        img = img.rotate(rotate, expand=True)
    return np.array(img)


def panel_letter(ax, letter, fontsize=16, color="black"):
    ax.text(
        0.025, 0.975, letter,
        transform=ax.transAxes,
        fontsize=fontsize, fontweight="bold", color=color,
        va="top", ha="left",
    )


def clean_ax(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


overview = load(OVERVIEW, rotate=90)
# Crop: remove left 30% and bottom 30%; right and top borders unchanged
_h, _w = overview.shape[:2]
overview = overview[0:int(0.7 * _h), int(0.3 * _w):]
sib1     = load(SIB1)
sib2     = load(SIB2)
wt       = load(WT)

fig = plt.figure(figsize=(12, 10), facecolor="white")
gs  = gridspec.GridSpec(3, 2, figure=fig,
                        width_ratios=[1.35, 1],
                        hspace=0.04, wspace=0.04,
                        left=0.01, right=0.99,
                        top=0.97, bottom=0.02)

ax_a = fig.add_subplot(gs[:, 0])   # full-height left
ax_b = fig.add_subplot(gs[0, 1])   # top right
ax_c = fig.add_subplot(gs[1, 1])   # middle right
ax_d = fig.add_subplot(gs[2, 1])   # bottom right (WT)

ax_a.imshow(overview)
ax_b.imshow(sib1)
ax_c.imshow(sib2)
ax_d.imshow(wt)

# Scale bar for panel A: 5 cm = 1058 px (10 cm = 2116 px in original)
h, w = overview.shape[:2]
bar_len = 1058
bar_x0  = w * 0.95 - bar_len
bar_y   = h * 0.94
bar_h   = h * 0.007
ax_a.add_patch(plt.Rectangle((bar_x0, bar_y), bar_len, bar_h,
                              color="white", zorder=5))
ax_a.text(bar_x0 + bar_len / 2, bar_y - h * 0.012, "5 cm",
          color="white", ha="center", va="bottom",
          fontsize=11, fontweight="bold", zorder=5)

def panel_label(ax, label, color="black"):
    ax.text(
        0.975, 0.975, label,
        transform=ax.transAxes,
        fontsize=13, fontweight="bold", color=color,
        va="top", ha="right",
    )

panel_letter(ax_a, "A", color="white")
panel_letter(ax_b, "B")
panel_letter(ax_c, "C")
panel_letter(ax_d, "D")
panel_label(ax_a, "1205-2 (1)", color="white")
panel_label(ax_b, "1205-2 (1)")
panel_label(ax_c, "1205-2 (2)")
panel_label(ax_d, "JI2822")

for ax in (ax_a, ax_b, ax_c, ax_d):
    clean_ax(ax)

fig.savefig(OUT, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Saved  {OUT.relative_to(REPO)}")
