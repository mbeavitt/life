#!/usr/bin/env python3
"""
make_fig1_fig3_variants.py

Generates:
  - figures/setup_composed.png  — final faceted Figure 1 (box left, windmills stacked right)
  - figures/setup_v{1,2,3}.png  — draft layout variants (kept for reference)
  - figures/curly_v{1,2,3}.png  — draft curly root variants (kept for reference)

Usage:
    python3 scripts/make_fig1_fig3_variants.py
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

# ── helpers ────────────────────────────────────────────────────────────────────

def load(path, crop=None):
    """Load image; crop = (left, top, right, bottom) in pixels."""
    img = Image.open(path).convert("RGB")
    if crop:
        img = img.crop(crop)
    return np.array(img)


def panel_letter(ax, letter, fontsize=18, color="black"):
    ax.text(
        0.025, 0.97, letter,
        transform=ax.transAxes,
        fontsize=fontsize, fontweight="bold", color=color,
        va="top", ha="left",
    )


def clean_ax(ax):
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])


def add_arrow(ax, xy_tip, xy_base, color="red", lw=2.5, head_width=20,
              outline=True):
    """Add an annotation arrow; coordinates are in image pixel space (data coords)."""
    props = dict(
        arrowstyle="-|>",
        color=color,
        lw=lw,
        mutation_scale=head_width,
    )
    if outline:
        # White outline underneath for contrast
        outline_props = dict(
            arrowstyle="-|>",
            color="white",
            lw=lw + 2.5,
            mutation_scale=head_width + 4,
        )
        ax.annotate("", xy=xy_tip, xytext=xy_base,
                    xycoords="data", textcoords="data",
                    arrowprops=outline_props)
    ax.annotate("", xy=xy_tip, xytext=xy_base,
                xycoords="data", textcoords="data",
                arrowprops=props)


# ── Figure 1 (setup) ──────────────────────────────────────────────────────────
#
# Source images (768 × 1024 each):
#   tip_box_growth_medium.jpeg      — tip box with etiolated seedlings
#   windmill_filter_paper.jpeg      — empty jar with windmill filter paper
#   windmill_filter_paper_inserted  — seedling threaded through aperture
#   windmill_filter_paper_tucked    — seedling in final position
#
# Crops (left, top, right, bottom) in pixels:
TIP_CROP_MOD   = (20,  80, 748, 960)   # moderate — removes counter at bottom
TIP_CROP_WIDE  = (0,   30, 768, 994)   # near-full frame
WIND_CROP_MOD  = (40, 120, 728, 960)   # moderate — tightens on jar
WIND_CROP_TIGHT = (80, 180, 690, 880)  # tighter — removes more surround

SETUP_IMGS = [
    FIG_DIR / "tip_box_growth_medium.jpeg",
    FIG_DIR / "windmill_filter_paper_inserted.jpeg",
    FIG_DIR / "windmill_filter_paper_tucked.jpeg",
]
LETTERS = "ABC"


def make_setup_v1():
    """1 × 3 row — shows sequence left to right."""
    fig, axes = plt.subplots(1, 3, figsize=(12, 5.2), facecolor="white")
    fig.subplots_adjust(wspace=0.025, hspace=0,
                        left=0.01, right=0.99, top=0.99, bottom=0.01)
    for idx, (ax, path) in enumerate(zip(axes, SETUP_IMGS)):
        ax.imshow(load(path))
        panel_letter(ax, LETTERS[idx], color="white")
        clean_ax(ax)
    return fig


def make_setup_v2():
    """2 × 2 grid, tip box top-left spanning two columns on left."""
    fig, axes = plt.subplots(1, 3, figsize=(12, 5.2), facecolor="white")
    fig.subplots_adjust(wspace=0.025, hspace=0,
                        left=0.01, right=0.99, top=0.99, bottom=0.01)
    for idx, (ax, path) in enumerate(zip(axes, SETUP_IMGS)):
        ax.imshow(load(path))
        panel_letter(ax, LETTERS[idx], color="white")
        clean_ax(ax)
    return fig


def make_setup_v3():
    """Tip box (A) on left, windmill stages (B–C) stacked on right."""
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(10, 6), facecolor="white")
    gs = GridSpec(2, 2, figure=fig, wspace=0.025, hspace=0.025,
                  left=0.01, right=0.99, top=0.99, bottom=0.01)
    ax_a = fig.add_subplot(gs[:, 0])
    ax_a.imshow(load(SETUP_IMGS[0]))
    panel_letter(ax_a, "A", color="white")
    clean_ax(ax_a)
    for i, (path, letter) in enumerate(zip(SETUP_IMGS[1:], "BC")):
        ax = fig.add_subplot(gs[i, 1])
        ax.imshow(load(path))
        panel_letter(ax, letter, color="white", fontsize=11)
        clean_ax(ax)
    return fig


# ── Figure 3 (curly root) ─────────────────────────────────────────────────────
#
# Source: figures/curly_lsh2oe.jpg (3024 × 4032)
#
# Arrow tip (arrowhead) targets the prominent S-curve hairy root in the
# centre of the image; arrow base is offset toward the upper-right.
# All coordinates are in pixels of the CROPPED image.
#
# Crop definitions (left, top, right, bottom):
CURLY_SRC = FIG_DIR / "curly_lsh2oe.jpg"

# v1 — wide: full context, minimal crop
CURLY_CROP_WIDE   = (150,  250, 2900, 3800)
# v2 — medium: exclude leaves (top-right) and label (bottom-right)
CURLY_CROP_MED    = (0,    600, 2400, 3350)
# v3 — tight: zoom into the root mass only
CURLY_CROP_TIGHT  = (0,   1000, 2100, 3100)

# Arrow tip in full-image pixels: approximately (880, 2200)
# (the S-curve hairy root in the centre-left of the image)
# Arrow base: pulled toward upper-right for a clean angle
ARROW_TIP_FULL  = (880,  2200)
ARROW_BASE_FULL = (1550, 1580)


def _adjust_arrow(tip_full, base_full, crop):
    """Translate full-image pixel coords to cropped-image pixel coords."""
    l, t, r, b = crop
    tip  = (tip_full[0]  - l, tip_full[1]  - t)
    base = (base_full[0] - l, base_full[1] - t)
    return tip, base


def make_curly(crop, title=None):
    """Single-panel curly root figure with red arrow."""
    img = load(CURLY_SRC, crop)
    h, w = img.shape[:2]

    fig, ax = plt.subplots(figsize=(6, 6 * h / w), facecolor="white")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    ax.imshow(img)
    clean_ax(ax)

    tip, base = _adjust_arrow(ARROW_TIP_FULL, ARROW_BASE_FULL, crop)
    add_arrow(ax, tip, base, color="#ff2222", lw=2.8, head_width=22)

    return fig


def make_curly_v1():
    """Wide crop — full context, shows whole plant + root."""
    return make_curly(CURLY_CROP_WIDE)


def make_curly_v2():
    """Medium crop — excludes shoot/leaves and label, focuses on root."""
    return make_curly(CURLY_CROP_MED)


def make_curly_v3():
    """Tight crop — zoomed into root mass, maximum detail on the curl."""
    return make_curly(CURLY_CROP_TIGHT)


# ── Figure 1 composed (final) ─────────────────────────────────────────────────
#
# Layout: tip box (A) portrait on left, windmill inserted (B) + tucked (C)
# stacked on right.  Width ratio is derived from the images' natural aspect
# ratios so that the right column height matches the left column height.
#
# tip_box:           w=739, h=867  → aspect w/h = 0.852
# windmill_inserted: w=722, h=719  → aspect w/h = 1.004
# windmill_tucked:   w=759, h=688  → aspect w/h = 1.103
#
# For equal total heights:
#   left_width  = H * 0.852
#   right_width = H/2 * mean(1.004, 1.103) ≈ H * 0.527
#   ratio left:right ≈ 1.617 : 1

def make_setup_composed():
    from matplotlib.gridspec import GridSpec

    tip   = load(SETUP_IMGS[0])   # 739 × 867
    wind1 = load(SETUP_IMGS[1])   # 722 × 719
    wind2 = load(SETUP_IMGS[2])   # 759 × 688

    h, w = tip.shape[:2]
    tip_aspect   = w / h           # 0.852
    w1_h, w1_w   = wind1.shape[:2]
    w2_h, w2_w   = wind2.shape[:2]
    right_aspect = ((w1_w / w1_h) + (w2_w / w2_h)) / 2  # ≈ 1.053, per image
    # right column width for equal column heights:
    right_col = right_aspect / 2   # each image gets half the total height

    fig_h = 8.0
    fig_w = fig_h * (tip_aspect + right_col) + 0.1  # +0.1 for gap

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")
    gs  = GridSpec(
        2, 2, figure=fig,
        width_ratios=[tip_aspect, right_col],
        height_ratios=[w1_h, w2_h],
        wspace=0.02, hspace=0.02,
        left=0.005, right=0.995, top=0.995, bottom=0.005,
    )

    ax_a = fig.add_subplot(gs[:, 0])   # tip box, spans both rows
    ax_b = fig.add_subplot(gs[0, 1])   # windmill inserted
    ax_c = fig.add_subplot(gs[1, 1])   # windmill tucked

    for ax, img, lbl in [(ax_a, tip, "A"), (ax_b, wind1, "B"), (ax_c, wind2, "C")]:
        ax.imshow(img)
        panel_letter(ax, lbl)
        clean_ax(ax)

    return fig


# ── main ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Final composed figure
    fig = make_setup_composed()
    out = FIG_DIR / "setup_composed.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved  {out.relative_to(REPO)}")

    # Draft variants (kept for reference)
    for tag, func in [
        ("setup_v1", make_setup_v1),
        ("setup_v2", make_setup_v2),
        ("setup_v3", make_setup_v3),
        ("curly_v1", make_curly_v1),
        ("curly_v2", make_curly_v2),
        ("curly_v3", make_curly_v3),
    ]:
        fig = func()
        out = FIG_DIR / f"{tag}.png"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved  {out.relative_to(REPO)}")
