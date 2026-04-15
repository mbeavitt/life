"""
add_scale_bars.py
-----------------
Adds white scale bars to calibrated microscopy/photography images for
use in report figures. Scale calibration is encoded in the source
filenames (e.g. "700px1cm" means 700 pixels = 1 cm).

Outputs are written to figures/ and named for use in report.typ.
M2 screen images (no calibration available) are copied as-is.

Run from the rotation_2 project root:
    python3 scripts/add_scale_bars.py
"""

from PIL import Image, ImageDraw, ImageFont
import shutil
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC  = os.path.join(BASE, "source_images")
FIG  = os.path.join(BASE, "figures")


def add_scale_bar(img_path, out_path, px_per_cm, bar_cm, label, dark_bg=True):
    """
    Overlay a scale bar on an image and save to out_path.

    Parameters
    ----------
    img_path   : str   Path to source image.
    out_path   : str   Destination path.
    px_per_cm  : float Pixels per centimetre (from filename calibration).
    bar_cm     : float Desired scale bar length in centimetres.
    label      : str   Text label rendered above the bar (e.g. "5 mm").
    dark_bg    : bool  If True, draw white bar/text; else draw black.
    """
    img = Image.open(img_path).convert("RGB")
    w, h = img.size
    draw = ImageDraw.Draw(img)

    bar_px = int(px_per_cm * bar_cm)
    bar_h  = max(10, h // 60)
    margin = max(20, h // 40)

    x2 = w - margin;  x1 = x2 - bar_px
    y2 = h - margin;  y1 = y2 - bar_h

    color = (255, 255, 255) if dark_bg else (0, 0, 0)
    draw.rectangle([x1, y1, x2, y2], fill=color)

    font = None
    for font_path in [
        "/System/Library/Fonts/Helvetica.ttc",
        "/Library/Fonts/Arial.ttf",
    ]:
        try:
            font = ImageFont.truetype(font_path, size=max(24, h // 25))
            break
        except OSError:
            continue
    if font is None:
        font = ImageFont.load_default()

    bbox = draw.textbbox((0, 0), label, font=font)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    draw.text((x1 + (bar_px - tw) // 2, y1 - th - 6), label, fill=color, font=font)

    img.save(out_path, quality=95)
    print(f"  {os.path.basename(img_path):50s} → {os.path.basename(out_path)}")


# ---------------------------------------------------------------------------
# Images with scale calibration
# (input_path, output_name, px_per_cm, bar_cm, bar_label, dark_background)
# ---------------------------------------------------------------------------
CALIBRATED = [
    # 614-2 M3 close-ups: 700 px = 1 cm  →  5 mm bar
    (f"{SRC}/614-2/614-2-1_CLOSE_M3_700px1cm.JPG",     "614-2-1_close_m3.jpg",  700,  0.5, "5 mm", True),
    (f"{SRC}/614-2/614-2-2_CLOSE_M3_700px1cm.JPG",     "614-2-2_close_m3.jpg",  700,  0.5, "5 mm", True),
    # 614-2 M3 far shots: 1000 px = 10 cm  →  5 cm bar
    (f"{SRC}/614-2/614-2-1_FAR_M3_1000px10cm.jpeg",    "614-2-1_far_m3.jpg",    100,  5,   "5 cm", True),
    (f"{SRC}/614-2/614-2-2_FAR_M3_1000px10cm.jpeg",    "614-2-2_far_m3.jpg",    100,  5,   "5 cm", True),
    # JI2822 close: 700 px = 1 cm  →  5 mm bar
    (f"{SRC}/JI2822/JI2822_700px1cm.JPG",              "ji2822_close.jpg",      700,  0.5, "5 mm", True),
    # JI2822 far: 1100 px = 10 cm  →  5 cm bar (white bar on any background)
    (f"{SRC}/JI2822/EMS_JI2822-5_FAR_1100px10cm.jpeg", "ji2822_far.jpg",        110,  5,   "5 cm", True),
    # Curly LSH2oe hairy root: 19 px/mm = 190 px/cm (cross-calibrated from fig 8)  →  5 cm bar
    (f"{FIG}/curly_lsh2oe.jpg",                         "curly_lsh2oe.jpg",      190,  5,   "5 cm", False),
]

# ---------------------------------------------------------------------------
# M2 screen images — no calibration, copy as-is
# ---------------------------------------------------------------------------
UNCALIBRATED = [
    (f"{SRC}/614-2/614-2_clusters_M2.jpg",   "614-2_clusters_m2.jpg"),
    (f"{SRC}/614-2/614-2_elongated2_M2.jpg", "614-2_elongated2_m2.jpg"),
]


if __name__ == "__main__":
    print("Adding scale bars to calibrated images...")
    for args in CALIBRATED:
        img_path, out_name = args[0], args[1]
        add_scale_bar(img_path, os.path.join(FIG, out_name), *args[2:])

    print("\nCopying uncalibrated M2 images...")
    for src_name, out_name in UNCALIBRATED:
        shutil.copy(src_name, os.path.join(FIG, out_name))
        print(f"  Copied {out_name}")

    print("\nDone.")
