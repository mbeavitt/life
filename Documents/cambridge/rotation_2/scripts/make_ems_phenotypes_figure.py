#!/usr/bin/env python3
"""
make_ems_phenotypes_figure.py

3×2 collage of representative atypical nodule morphologies from the M3 EMS
controlled-conditions experiment. Each panel shows a GUS-stained nodule image
with the line number below.

Output: figures/ems_phenotype_collage.png
"""

from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

REPO    = Path(__file__).resolve().parent.parent
SRC     = REPO / "source_images" / "cool random keyence images"
FIG_DIR = REPO / "figures"

CELL_H     = 700    # image height within each cell (px)
LABEL_H    = 80     # height of text strip below each image
GAP        = 18     # gap between cells
MARGIN     = 30
LABEL_SIZE = 28
PANEL_SIZE = 52
BG         = (255, 255, 255)

PANELS = [
    # (path, panel_letter, line_id)
    (
        SRC / "1077-1-2" / "2026_03_10_1077-1-2000001_cluster.jpg",
        "A", "1077-1 (2)",
    ),
    (
        SRC / "1132-2-1_all multilobed and wonky with odd splitting:branching phenotype"
            / "2026_03_10_1132-2-1000001_multilobed.jpg",
        "B", "1132-2 (1)",
    ),
    (
        SRC / "1132-2-1_all multilobed and wonky with odd splitting:branching phenotype"
            / "2026_03_10_1132-2-1000011_multilobed.jpg",
        "C", "1132-2 (1)",
    ),
    (
        SRC / "2026_03_10_1237-1-1000006_fanshaped.jpg",
        "D", "1237-1 (1)",
    ),
    (
        SRC / "2026_03_10_1310-1-2000001_hypernod_elongated.jpg",
        "E", "1310-1 (2)",
    ),
    (
        SRC / "2026_03_10_1310-1-3000009_clubshaped.jpg",
        "F", "1310-1 (3)",
    ),
]


def get_font(size):
    for path in ["/System/Library/Fonts/Helvetica.ttc", "/Library/Fonts/Arial.ttf"]:
        try:
            return ImageFont.truetype(path, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def make_cell(img_path, letter, line_id):
    """Return a PIL image of one collage cell: photo + line number strip below."""
    img = Image.open(img_path).convert("RGB")
    w, h = img.size
    new_w = int(w * CELL_H / h)
    img   = img.resize((new_w, CELL_H), Image.LANCZOS)

    cell = Image.new("RGB", (new_w, CELL_H + LABEL_H), BG)
    cell.paste(img, (0, 0))

    draw       = ImageDraw.Draw(cell)
    font_panel = get_font(PANEL_SIZE)
    font_label = get_font(LABEL_SIZE)

    # Panel letter — plain black text, top-left
    draw.text((14, 10), letter, fill=(0, 0, 0), font=font_panel)

    # Divider line between image and label strip
    draw.line([(0, CELL_H), (new_w, CELL_H)], fill=(200, 200, 200), width=1)

    # Label strip: line number centred
    bbox = draw.textbbox((0, 0), line_id, font=font_label)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    tx = (new_w - tw) // 2
    ty = CELL_H + (LABEL_H - th) // 2
    draw.text((tx, ty), line_id, fill=(40, 40, 40), font=font_label)

    return cell


def hstack(cells, gap=GAP):
    total_w = sum(c.width for c in cells) + gap * (len(cells) - 1)
    max_h   = max(c.height for c in cells)
    canvas  = Image.new("RGB", (total_w, max_h), BG)
    x = 0
    for c in cells:
        canvas.paste(c, (x, (max_h - c.height) // 2))
        x += c.width + gap
    return canvas


def vstack(rows, gap=GAP):
    max_w   = max(r.width for r in rows)
    total_h = sum(r.height for r in rows) + gap * (len(rows) - 1)
    canvas  = Image.new("RGB", (max_w, total_h), BG)
    y = 0
    for r in rows:
        canvas.paste(r, ((max_w - r.width) // 2, y))
        y += r.height + gap
    return canvas


def main():
    cells = []
    for path, letter, line_id in PANELS:
        print(f"  {letter}: {path.name}")
        cells.append(make_cell(path, letter, line_id))

    row1 = hstack(cells[:3])
    row2 = hstack(cells[3:])
    fig  = vstack([row1, row2], gap=GAP * 2)

    # Outer margin
    final = Image.new("RGB", (fig.width + 2 * MARGIN, fig.height + 2 * MARGIN), BG)
    final.paste(fig, (MARGIN, MARGIN))

    out = FIG_DIR / "ems_phenotype_collage.png"
    final.save(out, dpi=(300, 300))
    print(f"\nSaved: {out}  ({final.width} × {final.height} px)")


if __name__ == "__main__":
    main()
