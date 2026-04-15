#!/usr/bin/env python3
"""
make_614_figure.py
------------------
Assembles the controlled-conditions figure for EMS accession 614-2.

Layout
------
Row 1  (M2 JIC screen phenotype, no calibration):
  A  614-2_clusters_M2     B  614-2_elongated_M2

Row 2  (M3 controlled conditions — whole root system):
  C  614-2 (1) FAR          D  614-2 (2) FAR          E  JI2822 FAR

Row 3  (M3 controlled conditions — nodule close-up):
  F  614-2 (1) CLOSE        G  614-2 (2) CLOSE         H  JI2822 CLOSE

Scale calibration is read from source filenames:
  *_700px1cm*  → 700 px = 1 cm   → 5 mm scale bar  (350 px)
  *_1000px10cm* → 100 px/cm      → 5 cm scale bar  (500 px)
  JI2822 FAR: 1100px = 10 cm     → 5 cm scale bar  (550 px)

Output: figures/614-2_controlled_conditions.png
"""

from PIL import Image, ImageDraw, ImageFont
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC  = os.path.join(BASE, "source_images")
FIG  = os.path.join(BASE, "figures")

# ── Source images ─────────────────────────────────────────────────────────────

SOURCES = {
    "M2_clusters":   f"{SRC}/614-2/614-2_clusters_M2.jpg",
    "M2_elongated":  f"{SRC}/614-2/614-2_elongated_M2.png",
    "c1_far":        f"{SRC}/614-2/614-2-1_FAR_M3_1000px10cm.jpeg",
    "c2_far":        f"{SRC}/614-2/614-2-2_FAR_M3_1000px10cm.jpeg",
    "c1_close":      f"{SRC}/614-2/614-2-1_CLOSE_M3_700px1cm.JPG",
    "c2_close":      f"{SRC}/614-2/614-2-2_CLOSE_M3_700px1cm.JPG",
    "ji_far":        f"{SRC}/JI2822/EMS_JI2822-5_FAR_1100px10cm.jpeg",
    "ji_close":      f"{SRC}/JI2822/JI2822_700px1cm.JPG",
}

# (px_per_cm, bar_cm, label_text)
SCALE_BARS = {
    "c1_far":   (100, 5,   "5 cm"),
    "c2_far":   (100, 5,   "5 cm"),
    "ji_far":   (110, 5,   "5 cm"),
    "c1_close": (700, 0.5, "5 mm"),
    "c2_close": (700, 0.5, "5 mm"),
    "ji_close": (700, 0.5, "5 mm"),
}

# ── Target normalised heights (px) ────────────────────────────────────────────

M2_H    = 450    # M2 images displayed at this height
WIDE_H  = 1000   # All wide shots normalised to this height
CLOSE_H = 750    # All close-ups normalised to this height

# ── Helpers ───────────────────────────────────────────────────────────────────

def get_font(size, bold=False):
    candidates = [
        "/System/Library/Fonts/Helvetica.ttc",
        "/Library/Fonts/Arial.ttf",
    ]
    for path in candidates:
        try:
            return ImageFont.truetype(path, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def add_scale_bar(img: Image.Image, px_per_cm: float,
                  bar_cm: float, label: str) -> Image.Image:
    """Overlay a white scale bar with label onto img."""
    img = img.copy()
    draw = ImageDraw.Draw(img)
    w, h = img.size
    bar_px = int(px_per_cm * bar_cm)
    bar_h  = max(10, h // 60)
    margin = max(20, h // 40)
    x1, x2 = w - margin - bar_px, w - margin
    y2 = h - margin
    y1 = y2 - bar_h
    draw.rectangle([x1, y1, x2, y2], fill=(255, 255, 255))
    font = get_font(max(24, h // 25))
    bbox = draw.textbbox((0, 0), label, font=font)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    draw.text((x1 + (bar_px - tw) // 2, y1 - th - 6),
              label, fill=(255, 255, 255), font=font)
    return img


def add_panel_label(img: Image.Image, letter: str,
                    font_size: int = 32) -> Image.Image:
    """Black-box white-letter panel label in top-left corner."""
    img = img.copy()
    draw = ImageDraw.Draw(img)
    font = get_font(font_size)
    pad = 10
    bbox = draw.textbbox((0, 0), letter, font=font)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    x0, y0 = 14, 14
    draw.rectangle([x0 - pad, y0 - pad, x0 + tw + pad, y0 + th + pad],
                   fill=(0, 0, 0))
    draw.text((x0, y0), letter, fill=(255, 255, 255), font=font)
    return img


def resize_h(img: Image.Image, target_h: int) -> Image.Image:
    w, h = img.size
    return img.resize((int(w * target_h / h), target_h), Image.LANCZOS)


def hstack(images: list, gap: int = 15,
           bg: tuple = (255, 255, 255)) -> Image.Image:
    total_w = sum(i.width for i in images) + gap * (len(images) - 1)
    max_h   = max(i.height for i in images)
    canvas  = Image.new("RGB", (total_w, max_h), bg)
    x = 0
    for img in images:
        canvas.paste(img, (x, (max_h - img.height) // 2))
        x += img.width + gap
    return canvas


def vstack(rows: list, gaps: list, bg: tuple = (255, 255, 255)) -> Image.Image:
    """Stack rows vertically; gaps[i] is the space after row i."""
    max_w   = max(r.width for r in rows)
    total_h = sum(r.height for r in rows) + sum(gaps)
    canvas  = Image.new("RGB", (max_w, total_h), bg)
    y = 0
    for i, row in enumerate(rows):
        canvas.paste(row, ((max_w - row.width) // 2, y))
        y += row.height + (gaps[i] if i < len(gaps) else 0)
    return canvas


# ── Build figure ──────────────────────────────────────────────────────────────

def main():
    print("Loading images...")
    imgs = {k: Image.open(v).convert("RGB") for k, v in SOURCES.items()}
    for k, v in imgs.items():
        print(f"  {k:15s}  {v.width} x {v.height}")

    print("\nApplying scale bars...")
    for key, params in SCALE_BARS.items():
        imgs[key] = add_scale_bar(imgs[key], *params)

    print("Normalising and labelling panels...")
    # Use a single consistent font size across all panels so labels are the
    # same absolute size regardless of the row's normalised height.
    label_font = max(24, CLOSE_H // 22)   # ~34px; anchored to smallest row

    row1_panels = [
        add_panel_label(resize_h(imgs["M2_clusters"],  M2_H),  "A (614-2)",          label_font),
        add_panel_label(resize_h(imgs["M2_elongated"], M2_H),  "B (614-2)",          label_font),
    ]
    row2_panels = [
        add_panel_label(resize_h(imgs["c1_far"],  WIDE_H),  "C (614-2 (1))",         label_font),
        add_panel_label(resize_h(imgs["c2_far"],  WIDE_H),  "D (614-2 (2))",         label_font),
        add_panel_label(resize_h(imgs["ji_far"],  WIDE_H),  "E (JI2822 control)",    label_font),
    ]
    row3_panels = [
        add_panel_label(resize_h(imgs["c1_close"], CLOSE_H), "F (614-2 (1))",        label_font),
        add_panel_label(resize_h(imgs["c2_close"], CLOSE_H), "G (614-2 (2))",        label_font),
        add_panel_label(resize_h(imgs["ji_close"], CLOSE_H), "H (JI2822 control)",   label_font),
    ]

    row1 = hstack(row1_panels, gap=20)
    row2 = hstack(row2_panels, gap=15)
    row3 = hstack(row3_panels, gap=15)

    print(f"  Row 1 (M2):    {row1.width} x {row1.height}")
    print(f"  Row 2 (wide):  {row2.width} x {row2.height}")
    print(f"  Row 3 (close): {row3.width} x {row3.height}")

    # Larger gap between M2 row and M3 rows to visually separate them
    fig = vstack([row1, row2, row3], gaps=[60, 20])

    # Outer margin
    margin = 30
    final = Image.new("RGB",
                      (fig.width + 2 * margin, fig.height + 2 * margin),
                      (255, 255, 255))
    final.paste(fig, (margin, margin))

    # ── M2 / M3 annotation ───────────────────────────────────────────────────
    # Divider sits in the middle of the 60px gap between row1 and row2
    divider_y = margin + row1.height + 30

    draw = ImageDraw.Draw(final)

    # Full-width black rule
    draw.rectangle([0, divider_y - 2, final.width, divider_y + 2], fill=(0, 0, 0))

    # x-centre of annotation: middle of the blank strip to the right of row2
    row2_right = margin + (fig.width + row2.width) // 2
    ann_cx = (row2_right + final.width) // 2

    ann_font = get_font(80)
    offset  = 110   # distance from divider to centre of each annotation
    tri_w   = 28    # half-width of arrowhead
    tri_h   = 32    # height of arrowhead
    gap     = 8     # gap between arrowhead and text

    for label, above in [("M2", True), ("M3", False)]:
        bbox = draw.textbbox((0, 0), label, font=ann_font)
        tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
        total_h = tri_h + gap + th
        cy = divider_y - offset if above else divider_y + offset

        if above:
            # Triangle pointing up, text below
            tri_tip_y = cy - total_h // 2
            tri_base_y = tri_tip_y + tri_h
            txt_y = tri_base_y + gap
            tri_pts = [(ann_cx, tri_tip_y),
                       (ann_cx - tri_w, tri_base_y),
                       (ann_cx + tri_w, tri_base_y)]
        else:
            # Text above, triangle pointing down
            txt_y = cy - total_h // 2
            tri_base_y = txt_y + th + gap + 5
            tri_tip_y = tri_base_y + tri_h
            tri_pts = [(ann_cx, tri_tip_y),
                       (ann_cx - tri_w, tri_base_y),
                       (ann_cx + tri_w, tri_base_y)]

        draw.polygon(tri_pts, fill=(0, 0, 0))
        draw.text((ann_cx - tw // 2, txt_y), label, fill=(0, 0, 0), font=ann_font)

    out = os.path.join(FIG, "614-2_controlled_conditions.png")
    final.save(out, dpi=(300, 300))
    print(f"\nSaved: {out}  ({final.width} x {final.height} px)")


if __name__ == "__main__":
    main()
