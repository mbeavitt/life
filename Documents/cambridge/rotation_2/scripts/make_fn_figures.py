#!/usr/bin/env python3
"""
make_fn_figures.py

Generates two figures for the rotation report from the FN indel dataset:

  1. fn_priority_composite.png  — clean 3-panel locus plot for LSH1, NOOT1,
                                   NOOT2 with shared legend
  2. fn_dataset_overview.png    — 2×2 faceted overview of the full FN dataset

Outputs go to figures/.
"""

import gzip
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path

REPO         = Path(__file__).resolve().parent.parent
INDELS_GFF   = REPO / "fn_indels.gff3.gz"
ORTHOLOGS_GFF = REPO / "ortholog_genes.gff3.gz"
OUT          = REPO / "figures"

OZ_TO_CHR = {
    "OZ075428.1": "chr1", "OZ075429.1": "chr2", "OZ075430.1": "chr3",
    "OZ075431.1": "chr4", "OZ075432.1": "chr5", "OZ075433.1": "chr6",
    "OZ075434.1": "chr7",
}
CHR_SIZES_BP = {
    "chr1": 480540775, "chr2": 507322098, "chr3": 538651075,
    "chr4": 511007960, "chr5": 647395306, "chr6": 535013209, "chr7": 563422397,
}
CHROMS = [f"chr{i}" for i in range(1, 8)]

COL_DEL  = "#c0392b"
COL_DUP  = "#2471a3"
COL_NONE = "#bdc3c7"
COL_BOTH = "#8e44ad"

# Core priority loci with their gene names
PRIORITY_LOCI = [
    ("PSAT_LOCUS25697", "LSH1"),
    ("PSAT_LOCUS13168", "NOOT1"),
    ("PSAT_LOCUS26318", "NOOT2"),
]

sns.set_theme(style="whitegrid", font_scale=1.0)


# ── Helpers ───────────────────────────────────────────────────────────────────

def parse_attrs(attr_str):
    d = {}
    for part in attr_str.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
    return d


def parse_gff3(path, chrom_map=None):
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom = parts[0]
            if chrom_map:
                chrom = chrom_map.get(chrom, chrom)
            rows.append({
                "chrom": chrom,
                "type":  parts[2],
                "start": int(parts[3]),
                "end":   int(parts[4]),
                "attrs": parse_attrs(parts[8]),
            })
    return pd.DataFrame(rows)


# ── Load data ─────────────────────────────────────────────────────────────────

print("Loading data...")
indels = parse_gff3(INDELS_GFF)
indels["fn_line"]  = indels["attrs"].apply(lambda a: a.get("fn_line", ""))
indels["fn_sub"]   = indels["attrs"].apply(lambda a: a.get("fn_subline", ""))
indels["zygosity"] = indels["attrs"].apply(lambda a: a.get("zygosity", ""))
indels["size_kb"]  = (indels["end"] - indels["start"] + 1) / 1000
dels = indels[indels["type"] == "deletion"].copy()
dups = indels[indels["type"] == "duplication"].copy()

orthologs = parse_gff3(ORTHOLOGS_GFF, chrom_map=OZ_TO_CHR)
orthologs["symbol"] = orthologs["attrs"].apply(
    lambda a: next((v for v in [a.get("medicago_symbol", ""),
                                a.get("Name", ""), a.get("gene_id", "")] if v), "")
)
orthologs["locus"] = orthologs["attrs"].apply(lambda a: a.get("gene_id", ""))

# Intersect
hits = []
for _, gene in orthologs.iterrows():
    ov = indels[
        (indels["chrom"] == gene["chrom"]) &
        (indels["start"] <= gene["end"]) &
        (indels["end"]   >= gene["start"])
    ]
    for _, indel in ov.iterrows():
        hits.append({
            "locus":         gene["locus"],
            "gene_chrom":    gene["chrom"],
            "gene_start":    gene["start"],
            "gene_end":      gene["end"],
            "indel_type":    indel["type"],
            "zygosity":      indel["zygosity"],
            "fn_line":       indel["fn_line"],
            "fn_sub":        indel["fn_sub"],
            "indel_start":   indel["start"],
            "indel_end":     indel["end"],
            "indel_size_kb": indel["size_kb"],
        })
hits_df = pd.DataFrame(hits)

hit_loci_del = set(hits_df[hits_df["indel_type"] == "deletion"]["locus"])
hit_loci_dup = set(hits_df[hits_df["indel_type"] == "duplication"]["locus"])

def hit_status(locus):
    d, u = locus in hit_loci_del, locus in hit_loci_dup
    if d and u: return "both"
    if d:       return "deletion"
    if u:       return "duplication"
    return "none"

orthologs["hit_status"] = orthologs["locus"].apply(hit_status)
print(f"  {len(indels)} indels, {len(orthologs)} orthologs, {len(hits_df)} overlaps")


# ── Figure 1: priority composite ──────────────────────────────────────────────

print("Making priority composite...")

fig, axes = plt.subplots(3, 1, figsize=(10, 7),
                          gridspec_kw={"hspace": 0.55})

for ax, (locus, gene_name) in zip(axes, PRIORITY_LOCI):
    gene      = orthologs[orthologs["locus"] == locus].iloc[0]
    # Collapse to homozygous calls only
    gene_hits = (hits_df[(hits_df["locus"] == locus) &
                          (hits_df["zygosity"] == "homozygous")]
                 .sort_values("indel_type").reset_index(drop=True))

    x_min = min(gene_hits["indel_start"].min(), gene["start"])
    x_max = max(gene_hits["indel_end"].max(),   gene["end"])
    span  = x_max - x_min
    pad   = span * 0.06

    bar_h = 0.38   # narrowed by ~30%

    gene_display_w = max(gene["end"] - gene["start"], span * 0.012)
    gene_display_s = (gene["start"] + gene["end"]) / 2 - gene_display_w / 2
    gene_right_mb  = (gene_display_s + gene_display_w) / 1e6

    # Deletion bar
    for i, row in gene_hits.iterrows():
        ax.barh(i, (row["indel_end"] - row["indel_start"]) / 1e6,
                left=row["indel_start"] / 1e6,
                height=bar_h, color=COL_DEL, alpha=0.88, zorder=2)
        bar_mid = (row["indel_start"] + row["indel_end"]) / 2 / 1e6
        ax.text(bar_mid, i, f"{row['indel_size_kb']:.0f} kb",
                ha="center", va="center", fontsize=7.5,
                color="white", fontweight="bold", zorder=3)

    # Gene body — solid black bar, same height as deletion bar, on top
    ax.barh(0, gene_display_w / 1e6, left=gene_display_s / 1e6,
            height=bar_h, color="black", zorder=4)

    ax.set_xlim((x_min - pad) / 1e6, (x_max + pad) / 1e6)
    ax.set_ylim(-0.5, 0.5)
    ax.set_yticks([])
    ax.set_xlabel(f"{gene['chrom']} position (Mb)", fontsize=8.5)
    sns.despine(ax=ax, left=True)

    # Gene name to the right of the gene body with a short slanted pointer
    x_label = gene_right_mb + span * 0.09 / 1e6
    ax.annotate(
        gene_name,
        xy=(gene_right_mb, 0),
        xytext=(x_label, 0.33),
        fontsize=11, fontweight="bold",
        ha="left", va="center",
        arrowprops=dict(
            arrowstyle="-",
            connectionstyle="arc3,rad=-0.2",
            color="black", lw=1.0,
        ),
        annotation_clip=False,
    )

# Shared legend
legend_handles = [
    mpatches.Patch(color=COL_DEL, alpha=0.88, label="Deletion"),
    mpatches.Patch(color="black", label="Gene body"),
]
fig.legend(handles=legend_handles, loc="upper right",
           bbox_to_anchor=(0.98, 0.98), fontsize=9, framealpha=0.9)

out1 = OUT / "fn_priority_composite.png"
fig.savefig(out1, dpi=180, bbox_inches="tight")
plt.close()
print(f"  → {out1.name}")


# ── Figure 2: dataset overview ────────────────────────────────────────────────

print("Making dataset overview...")

fig, axes = plt.subplots(2, 2, figsize=(12, 8),
                          gridspec_kw={"hspace": 0.42, "wspace": 0.32})
ax_size, ax_line, ax_chrom, ax_map = axes.flat

# ── A: size distributions ─────────────────────────────────────────────────────
log_bins = np.logspace(
    np.log10(max(indels["size_kb"].min(), 1)),
    np.log10(indels["size_kb"].max()), 55
)
ax_size.hist(dels["size_kb"], bins=log_bins, alpha=0.65, color=COL_DEL,
             label=f"Deletions (n={len(dels):,})")
ax_size.hist(dups["size_kb"], bins=log_bins, alpha=0.65, color=COL_DUP,
             label=f"Duplications (n={len(dups):,})")
ax_size.set_xscale("log")
ax_size.set_xlabel("Indel size (kb, log scale)")
ax_size.set_ylabel("Count")
ax_size.set_title("A   Size distributions", loc="left", fontweight="bold")
ax_size.legend(fontsize=8, framealpha=0.8)
sns.despine(ax=ax_size)

# ── B: indels per line distribution ──────────────────────────────────────────
per_line = indels.groupby("fn_line").size()
ax_line.hist(per_line, bins=30, color="#7f8c8d", edgecolor="white", linewidth=0.4)
ax_line.axvline(per_line.median(), color="black", linestyle="--", linewidth=1,
                label=f"Median = {per_line.median():.0f}")
ax_line.set_xlabel("Indels per FN line")
ax_line.set_ylabel("Number of lines")
ax_line.set_title("B   Indels per line", loc="left", fontweight="bold")
ax_line.legend(fontsize=8, framealpha=0.8)
sns.despine(ax=ax_line)

# ── C: per-chromosome Mb ──────────────────────────────────────────────────────
chrom_mb = (
    indels.assign(size_mb=lambda d: d["size_kb"] / 1000)
    .groupby(["chrom", "type"])["size_mb"].sum()
    .unstack(fill_value=0).reindex(CHROMS)
)
x = np.arange(len(CHROMS))
w = 0.38
ax_chrom.bar(x - w/2, chrom_mb.get("deletion", 0), width=w,
             color=COL_DEL, alpha=0.85, label="Deletions")
ax_chrom.bar(x + w/2, chrom_mb.get("duplication", 0), width=w,
             color=COL_DUP, alpha=0.85, label="Duplications")
ax_chrom.set_xticks(x)
ax_chrom.set_xticklabels(CHROMS, fontsize=8.5)
ax_chrom.set_xlabel("Chromosome")
ax_chrom.set_ylabel("Total sequence affected (Mb)")
ax_chrom.set_title("C   FN indel coverage per chromosome", loc="left", fontweight="bold")
ax_chrom.legend(fontsize=8, framealpha=0.8)
sns.despine(ax=ax_chrom)

# ── D: genome-wide ortholog hit map ───────────────────────────────────────────
STATUS_COL = {"none": COL_NONE, "deletion": COL_DEL,
              "duplication": COL_DUP, "both": COL_BOTH}
STATUS_Z   = {"none": 1, "deletion": 3, "duplication": 3, "both": 4}
chr_y = {c: i for i, c in enumerate(reversed(CHROMS))}

for _, gene in orthologs.iterrows():
    if gene["chrom"] not in chr_y:
        continue
    y   = chr_y[gene["chrom"]]
    mid = (gene["start"] + gene["end"]) / 2 / 1e6
    ax_map.plot(mid, y, "|",
                color=STATUS_COL[gene["hit_status"]],
                markersize=9, markeredgewidth=1.6,
                zorder=STATUS_Z[gene["hit_status"]])
for chrom, y in chr_y.items():
    ax_map.hlines(y, 0, CHR_SIZES_BP[chrom] / 1e6,
                  color="#aab7b8", linewidth=0.8, zorder=0)

ax_map.set_yticks(list(chr_y.values()))
ax_map.set_yticklabels(list(reversed(CHROMS)), fontsize=8.5)
ax_map.set_xlabel("Position (Mb)")
ax_map.set_xlim(0, max(CHR_SIZES_BP.values()) / 1e6 + 5)
ax_map.set_title("D   Ortholog gene hit map", loc="left", fontweight="bold")
ax_map.legend(
    handles=[mpatches.Patch(color=STATUS_COL[s],
             label={"none": "No hit", "deletion": "Deletion",
                    "duplication": "Duplication", "both": "Both"}[s])
             for s in ["deletion", "duplication", "both", "none"]],
    fontsize=8, framealpha=0.9, loc="lower right",
)
sns.despine(ax=ax_map, left=True)

out2 = OUT / "fn_dataset_overview.png"
fig.savefig(out2, dpi=180, bbox_inches="tight")
plt.close()
print(f"  → {out2.name}")

print("Done.")
