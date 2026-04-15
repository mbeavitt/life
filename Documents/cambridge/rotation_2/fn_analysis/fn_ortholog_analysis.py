"""
FN Indel × Nodulation Ortholog Intersection Analysis
=====================================================
Intersects fast neutron deletion/duplication calls with nodulation ortholog
gene models in Pisum sativum JI2822. Generates per-gene locus plots and
compiles a PDF report via typst.

Usage:
    python fn_ortholog_analysis.py
"""

import gzip
import os
import subprocess
from datetime import date

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns

# ── Config ────────────────────────────────────────────────────────────────────

BROWSER_DIR   = "/home/mbeavitt/Documents/ks_analysis/browser/ji2822"
ANALYSIS_DIR  = "/home/mbeavitt/Documents/ks_analysis"
INDELS_GFF    = f"{BROWSER_DIR}/fn_indels.gff3.gz"
ORTHOLOGS_GFF = f"{BROWSER_DIR}/ortholog_genes.gff3.gz"
PLOT_DIR      = f"{ANALYSIS_DIR}/plots"
TYP_PATH      = f"{ANALYSIS_DIR}/fn_indel_ortholog_report.typ"
PDF_PATH      = f"{ANALYSIS_DIR}/fn_indel_ortholog_report.pdf"

os.makedirs(PLOT_DIR, exist_ok=True)

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

PRIORITY_GENES = ["LSH1", "LSH2", "NOOT1", "NOOT2", "NF-YA1"]

COL_DEL  = "#c0392b"
COL_DUP  = "#2471a3"
COL_BOTH = "#8e44ad"
COL_NONE = "#bdc3c7"

sns.set_theme(style="whitegrid", font_scale=1.05)

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
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
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


def typ_escape(s):
    """Escape special characters for typst content blocks."""
    return (str(s)
            .replace("\\", "\\\\")
            .replace("#", "\\#")
            .replace("@", "\\@")
            .replace("]", "\\]")
            .replace("[", "\\[")
            .replace("<", "\\<")
            .replace(">", "\\>"))


def typ_cell(s):
    return f"[{typ_escape(s)}]"


# ── Load data ─────────────────────────────────────────────────────────────────

print("Loading FN indels...")
indels = parse_gff3(INDELS_GFF)
indels["fn_line"]  = indels["attrs"].apply(lambda a: a.get("fn_line", ""))
indels["fn_sub"]   = indels["attrs"].apply(lambda a: a.get("fn_subline", ""))
indels["zygosity"] = indels["attrs"].apply(lambda a: a.get("zygosity", ""))
indels["size_kb"]  = (indels["end"] - indels["start"] + 1) / 1000
dels = indels[indels["type"] == "deletion"].copy()
dups = indels[indels["type"] == "duplication"].copy()

print("Loading ortholog genes...")
orthologs = parse_gff3(ORTHOLOGS_GFF, chrom_map=OZ_TO_CHR)
orthologs["symbol"] = orthologs["attrs"].apply(
    lambda a: next((v for v in [a.get("medicago_symbol", ""), a.get("Name", ""), a.get("gene_id", "")] if v), "")
)
orthologs["description"] = orthologs["attrs"].apply(lambda a: a.get("description", ""))
orthologs["locus"]       = orthologs["attrs"].apply(lambda a: a.get("gene_id", ""))

print(f"  {len(indels)} indels ({len(dels)} deletions, {len(dups)} duplications)")
print(f"  {len(orthologs)} ortholog genes across {indels['fn_line'].nunique()} FN lines")

# ── Intersection ──────────────────────────────────────────────────────────────

print("Intersecting...")
hits = []
for _, gene in orthologs.iterrows():
    ov = indels[
        (indels["chrom"] == gene["chrom"]) &
        (indels["start"] <= gene["end"]) &
        (indels["end"]   >= gene["start"])
    ]
    for _, indel in ov.iterrows():
        hits.append({
            "gene_name":     gene["symbol"],
            "locus":         gene["locus"],
            "gene_chrom":    gene["chrom"],
            "gene_start":    gene["start"],
            "gene_end":      gene["end"],
            "description":   gene["description"],
            "indel_type":    indel["type"],
            "zygosity":      indel["zygosity"],
            "fn_line":       indel["fn_line"],
            "fn_sub":        indel["fn_sub"],
            "indel_start":   indel["start"],
            "indel_end":     indel["end"],
            "indel_size_kb": indel["size_kb"],
        })

hits_df = pd.DataFrame(hits)
print(f"  {len(hits_df)} overlaps across {hits_df['gene_name'].nunique()} genes "
      f"and {hits_df['fn_line'].nunique()} FN lines")

gene_summary = (
    hits_df
    .groupby(["gene_name", "locus", "gene_chrom", "gene_start", "gene_end", "description"])
    .apply(lambda g: pd.Series({
        "n_lines_del":   g[g["indel_type"] == "deletion"]["fn_line"].nunique(),
        "n_lines_dup":   g[g["indel_type"] == "duplication"]["fn_line"].nunique(),
        "n_lines_total": g["fn_line"].nunique(),
    }), include_groups=False)
    .reset_index()
    .sort_values("n_lines_total", ascending=False)
)

hit_loci_del = set(hits_df[hits_df["indel_type"] == "deletion"]["locus"])
hit_loci_dup = set(hits_df[hits_df["indel_type"] == "duplication"]["locus"])

def hit_status(locus):
    d, u = locus in hit_loci_del, locus in hit_loci_dup
    if d and u: return "both"
    if d:       return "deletion"
    if u:       return "duplication"
    return "none"

orthologs["hit_status"] = orthologs["locus"].apply(hit_status)

# ── Per-gene locus plots ──────────────────────────────────────────────────────

print("Generating per-gene locus plots...")

def plot_gene_locus(gene_name, locus, gene_chrom, gene_start, gene_end,
                    gene_hits, outpath):
    """
    Horizontal bar chart showing each overlapping FN indel to scale,
    with the gene body highlighted underneath.
    """
    # Deduplicate: group hemi+homo pairs for the same subline+type into one row
    # (keep both as separate bars so the size difference is visible)
    rows = gene_hits.sort_values(["indel_type", "fn_sub", "zygosity"]).reset_index(drop=True)
    n = len(rows)

    fig_h = max(2.0, 0.45 * n + 1.0)
    fig, ax = plt.subplots(figsize=(9, fig_h))

    x_min = min(rows["indel_start"].min(), gene_start)
    x_max = max(rows["indel_end"].max(),   gene_end)
    span  = x_max - x_min
    pad   = span * 0.05
    x_min -= pad
    x_max += pad

    # Minimum visible gene width: at least 1% of the total x span
    gene_display_width = max(gene_end - gene_start, span * 0.01)
    gene_display_start = (gene_start + gene_end) / 2 - gene_display_width / 2

    # Gene body: shaded rectangle spanning full y range as background guide
    ax.axvspan(gene_display_start / 1e6, (gene_display_start + gene_display_width) / 1e6,
               color="black", alpha=0.12, zorder=1, label="_nolegend_")
    # Gene body bar at y = -0.6 (below indel rows)
    gene_bar_y = -0.6
    ax.barh(gene_bar_y, gene_display_width / 1e6,
            left=gene_display_start / 1e6, height=0.35,
            color="black", zorder=3)
    ax.text((gene_display_start + gene_display_width / 2) / 1e6, gene_bar_y,
            gene_name, ha="center", va="center",
            fontsize=7, color="white", fontweight="bold", zorder=4)

    # Indel bars (no bar labels — see colour legend and y-tick labels)
    for i, row in rows.iterrows():
        col   = COL_DEL if row["indel_type"] == "deletion" else COL_DUP
        alpha = 0.45 if row["zygosity"] == "hemizygous" else 0.85
        ax.barh(i, (row["indel_end"] - row["indel_start"]) / 1e6,
                left=row["indel_start"] / 1e6,
                height=0.6, color=col, alpha=alpha, zorder=2)

    ax.set_xlim(x_min / 1e6, x_max / 1e6)
    ax.set_ylim(-1.1, n)
    ax.set_xlabel(f"{gene_chrom} position (Mb)")
    # Y-tick labels: subline + zygosity + size
    ytick_labels = [
        f"{r['fn_sub']}  {'hemi' if r['zygosity']=='hemizygous' else r['zygosity']}  {r['indel_size_kb']:.0f} kb"
        for _, r in rows.iterrows()
    ]
    ax.set_yticks(range(n))
    ax.set_yticklabels(ytick_labels, fontsize=7.5)
    gene_kb = (gene_end - gene_start + 1) / 1000
    ax.set_title(f"{gene_name}  ·  {locus}  ·  {gene_chrom}:{gene_start/1e6:.3f}–{gene_end/1e6:.3f} Mb  ({gene_kb:.1f} kb gene)",
                 fontsize=9)

    legend_handles = [
        mpatches.Patch(color=COL_DEL, alpha=0.7, label="Deletion"),
        mpatches.Patch(color=COL_DUP, alpha=0.7, label="Duplication"),
        mpatches.Patch(color="black", alpha=0.5, label="Gene body"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=8, framealpha=0.9)

    sns.despine(left=True)
    plt.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()


# Generate a plot for every hit gene
gene_plot_paths = {}  # locus → relative path from ANALYSIS_DIR
for _, row in gene_summary.iterrows():
    locus     = row["locus"]
    safe_name = locus.replace(":", "_").replace("/", "_")
    fname     = f"plots/gene_{safe_name}.png"
    fpath     = f"{ANALYSIS_DIR}/{fname}"
    gene_hits = hits_df[hits_df["locus"] == locus]
    plot_gene_locus(
        row["gene_name"], locus, row["gene_chrom"],
        int(row["gene_start"]), int(row["gene_end"]),
        gene_hits, fpath
    )
    gene_plot_paths[locus] = fname

print(f"  {len(gene_plot_paths)} locus plots saved.")

# ── Summary stat plots ────────────────────────────────────────────────────────

print("Generating summary plots...")

# (a) Dels + dups per FN line (top 30)
line_counts = (
    indels.groupby(["fn_line", "type"])
    .size().unstack(fill_value=0)
    .assign(total=lambda d: d.sum(axis=1))
    .sort_values("total", ascending=False)
    .head(30)
    .drop(columns="total")
)
fig, ax = plt.subplots(figsize=(14, 5))
line_counts.plot(kind="bar", stacked=True, ax=ax,
                 color=[COL_DEL, COL_DUP], width=0.8)
ax.set_xlabel("FN line"); ax.set_ylabel("Number of calls")
ax.set_title("Deletions and duplications per FN line (top 30 by total count)")
ax.legend(["Deletion", "Duplication"])
plt.xticks(rotation=45, ha="right", fontsize=8)
plt.tight_layout()
fig.savefig(f"{PLOT_DIR}/fn_stats_dels_dups_per_line.png", dpi=150); plt.close()

# (b) Deletion size distribution (log x)
fig, ax = plt.subplots(figsize=(9, 5))
log_bins = np.logspace(np.log10(max(dels["size_kb"].min(), 1)), np.log10(dels["size_kb"].max()), 60)
for zyg, col, alpha in [("hemizygous", COL_DEL, 0.45), ("homozygous", "#922b21", 0.55)]:
    subset = dels[dels["zygosity"] == zyg]["size_kb"]
    ax.hist(subset, bins=log_bins, alpha=alpha, color=col, label=f"{zyg} (n={len(subset):,})")
ax.set_xscale("log"); ax.set_xlabel("Deletion size (kb, log scale)"); ax.set_ylabel("Count")
ax.set_title("Deletion size distribution by zygosity"); ax.legend()
plt.tight_layout(); fig.savefig(f"{PLOT_DIR}/fn_stats_del_size_dist.png", dpi=150); plt.close()

# (c) Duplication size distribution (log x)
fig, ax = plt.subplots(figsize=(9, 5))
log_bins_dup = np.logspace(np.log10(max(dups["size_kb"].min(), 1)), np.log10(dups["size_kb"].max()), 60)
ax.hist(dups["size_kb"], bins=log_bins_dup, alpha=0.6, color=COL_DUP, label=f"n={len(dups):,}")
ax.set_xscale("log"); ax.set_xlabel("Duplication size (kb, log scale)"); ax.set_ylabel("Count")
ax.set_title("Duplication size distribution"); ax.legend()
plt.tight_layout(); fig.savefig(f"{PLOT_DIR}/fn_stats_dup_size_dist.png", dpi=150); plt.close()

# (d) Per-chromosome Mb deleted / duplicated
chrom_mb = (
    indels.assign(size_mb=indels["size_kb"] / 1000)
    .groupby(["chrom", "type"])["size_mb"].sum()
    .unstack(fill_value=0).reindex(CHROMS)
)
fig, ax = plt.subplots(figsize=(9, 5))
chrom_mb.plot(kind="bar", ax=ax, color=[COL_DEL, COL_DUP], width=0.7)
ax.set_xlabel("Chromosome"); ax.set_ylabel("Total Mb")
ax.set_title("Total Mb deleted and duplicated per chromosome")
ax.legend(["Deletion", "Duplication"]); plt.xticks(rotation=0)
plt.tight_layout(); fig.savefig(f"{PLOT_DIR}/fn_stats_genome_coverage.png", dpi=150); plt.close()

# (e) Deletion size violin by chromosome
fig, ax = plt.subplots(figsize=(10, 5))
data_by_chr = [dels[dels["chrom"] == c]["size_kb"].values for c in CHROMS]
vp = ax.violinplot(data_by_chr, positions=range(7), showmedians=True, showextrema=True)
for body in vp["bodies"]:
    body.set_facecolor(COL_DEL); body.set_alpha(0.6)
ax.set_xticks(range(7)); ax.set_xticklabels(CHROMS)
ax.set_xlabel("Chromosome"); ax.set_ylabel("Deletion size (kb)")
ax.set_title("Deletion size distribution per chromosome")
plt.tight_layout(); fig.savefig(f"{PLOT_DIR}/fn_stats_del_size_by_chrom.png", dpi=150); plt.close()

# (f) Genome-wide ortholog coverage map
STATUS_COLOUR = {"none": COL_NONE, "deletion": COL_DEL, "duplication": COL_DUP, "both": COL_BOTH}
STATUS_LABEL  = {"none": "No hit", "deletion": "Deletion", "duplication": "Duplication", "both": "Both"}
STATUS_Z      = {"none": 1, "deletion": 3, "duplication": 3, "both": 4}

fig, ax = plt.subplots(figsize=(14, 5))
chr_y = {c: i for i, c in enumerate(reversed(CHROMS))}
for _, gene in orthologs.iterrows():
    if gene["chrom"] not in chr_y: continue
    y   = chr_y[gene["chrom"]]
    mid = (gene["start"] + gene["end"]) / 2 / 1e6
    ax.plot(mid, y, "|", color=STATUS_COLOUR[gene["hit_status"]],
            markersize=10, markeredgewidth=1.8, zorder=STATUS_Z[gene["hit_status"]])
for chrom, y in chr_y.items():
    ax.hlines(y, 0, CHR_SIZES_BP[chrom] / 1e6, color="#aab7b8", linewidth=1, zorder=0)
ax.set_yticks(list(chr_y.values())); ax.set_yticklabels(list(reversed(CHROMS)))
ax.set_xlabel("Position (Mb)")
ax.set_title("Nodulation ortholog genes coloured by FN indel hit status")
ax.set_xlim(0, max(CHR_SIZES_BP.values()) / 1e6 + 5)
ax.legend(handles=[mpatches.Patch(color=STATUS_COLOUR[s], label=STATUS_LABEL[s])
                   for s in ["none", "deletion", "duplication", "both"]],
          loc="lower right", framealpha=0.9)
sns.despine(left=True)
plt.tight_layout(); fig.savefig(f"{PLOT_DIR}/fn_stats_ortholog_map.png", dpi=150); plt.close()

print("  All summary plots saved.")

# ── Summary statistics ────────────────────────────────────────────────────────

n_lines      = indels["fn_line"].nunique()
total_del_mb = dels["size_kb"].sum() / 1000
total_dup_mb = dups["size_kb"].sum() / 1000
n_genes_del  = hits_df[hits_df["indel_type"] == "deletion"]["locus"].nunique()
n_genes_dup  = hits_df[hits_df["indel_type"] == "duplication"]["locus"].nunique()
n_genes_any  = hits_df["locus"].nunique()

# ── Build typst document ──────────────────────────────────────────────────────

print("Writing typst source...")

T = []
def t(*args): T.append(" ".join(str(a) for a in args))

t(f"""\
#set document(title: "FN Indel × Nodulation Ortholog Intersection Report")
#set page(paper: "a4", margin: (x: 2.5cm, y: 2.5cm))
#set text(font: "Libertinus Serif", size: 10pt)
#set heading(numbering: "1.")
#set par(justify: true)
#show table: set text(size: 9pt)
#show figure.caption: set text(size: 9pt, style: "italic")

#align(center)[
  #v(2em)
  #text(size: 20pt, weight: "bold")[FN Indel × Nodulation Ortholog \\ Intersection Report]
  #v(0.6em)
  #text(size: 13pt)[_Pisum sativum_ JI2822]
  #v(0.4em)
  #text(size: 10pt)[{date.today()}]
  #v(2em)
]

#outline(depth: 2)
#pagebreak()
""")

# Methods
t("= Methods\n")
t(f"Fast neutron (FN) deletion and duplication calls ({len(dels)} deletions, {len(dups)} duplications "
  f"across {n_lines} FN lines) were intersected with {len(orthologs)} nodulation ortholog gene models "
  f"mapped to the _Pisum sativum_ JI2822 v1.3 assembly. FN indel coordinates are from skim-sequencing "
  f"coverage analysis (Supplementary Tables 5--6). Ortholog coordinates were mapped from "
  f"_Medicago truncatula_ A17 via synteny. An overlap is defined as any FN interval sharing ≥1 bp "
  f"with a gene model. Hemizygous and homozygous deletion calls are treated independently (a single "
  f"subline may contribute both, representing the outer and inner coverage thresholds of the same "
  f"deletion event). Chromosome aliases (chr1--7 ↔ OZ075428--34.1) were resolved via chromAliases.txt.\n")

# Summary stats
t("= Summary Statistics\n")
t(f"""#figure(
  table(
    columns: (2fr, 1fr),
    table.header([*Metric*], [*Value*]),
    [FN lines], [{n_lines}],
    [Total deletions], [{len(dels)}],
    [Total duplications], [{len(dups)}],
    [Mean deletion size], [{dels['size_kb'].mean():.0f} kb],
    [Median deletion size], [{dels['size_kb'].median():.0f} kb],
    [Mean duplication size], [{dups['size_kb'].mean():.0f} kb],
    [Median duplication size], [{dups['size_kb'].median():.0f} kb],
    [Total Mb deleted (genome-wide)], [{total_del_mb:.1f} Mb],
    [Total Mb duplicated (genome-wide)], [{total_dup_mb:.1f} Mb],
    [Ortholog genes hit by ≥1 deletion], [{n_genes_del}],
    [Ortholog genes hit by ≥1 duplication], [{n_genes_dup}],
    [Ortholog genes hit by any indel], [{n_genes_any} / {len(orthologs)}],
    [Ortholog genes with no hit], [{len(orthologs) - n_genes_any}],
  ),
  caption: [Summary of FN indel calls and ortholog intersections.],
)
""")

t("""#figure(
  image("plots/fn_stats_ortholog_map.png", width: 100%),
  caption: [Genome-wide distribution of the 488 nodulation ortholog genes coloured by FN indel hit status. Grey: no hit; red: ≥1 deletion; blue: ≥1 duplication; purple: both.],
)
#pagebreak()
""")

# Priority genes
t("= Priority Genes\n")
for pg in PRIORITY_GENES:
    gene_rows = orthologs[orthologs["symbol"].str.contains(
        pg.replace("Mt", ""), case=False, na=False)]
    t(f"== {typ_escape(pg)}\n")
    if gene_rows.empty:
        t("_Not found in the ortholog gene set._\n"); continue
    for _, gr in gene_rows.iterrows():
        gene_hits = hits_df[hits_df["locus"] == gr["locus"]]
        gene_kb   = (gr["end"] - gr["start"] + 1) / 1000
        t(f"*Locus:* {typ_escape(gr['locus'])} · "
          f"*Position:* {gr['chrom']}:{gr['start']:,}--{gr['end']:,} ({gene_kb:.1f} kb) · "
          f"*Description:* _{typ_escape(gr['description'] or '—')}_\n")
        if gene_hits.empty:
            t("No FN indels overlap this gene.\n")
        else:
            n_l = gene_hits["fn_line"].nunique()
            t(f"{len(gene_hits)} hit(s) across {n_l} FN line(s):\n")
            t("""#table(
  columns: (1.5fr, 1fr, 1.2fr, 2.5fr, 1fr),
  table.header([*FN subline*], [*Type*], [*Zygosity*], [*Interval*], [*Size (kb)*]),""")
            for _, h in gene_hits.sort_values(["indel_type", "fn_sub"]).iterrows():
                interval = f"{h['gene_chrom']}:{h['indel_start']:,}--{h['indel_end']:,}"
                size_str = f"{h['indel_size_kb']:.0f}"
                t(f"  {typ_cell(h['fn_sub'])}, {typ_cell(h['indel_type'])}, "
                  f"{typ_cell(h['zygosity'])}, {typ_cell(interval)}, {typ_cell(size_str)},")
            t(")\n")
            plot_rel = gene_plot_paths[gr["locus"]]
            t(f"""#figure(
  image("{plot_rel}", width: 95%),
  caption: [Locus plot for {typ_escape(gr['symbol'])}. Each bar is one FN indel drawn to genomic scale; y-axis labels give subline, zygosity, and size. Opacity distinguishes hemizygous (lighter) from homozygous (darker). Gene body shown as black bar (minimum display width applied where gene is too small to see at this scale).],
)
""")

t("#pagebreak()")

# All candidate genes
t("= All Candidate Genes\n")
t(f"All {len(gene_summary)} ortholog genes overlapped by at least one FN indel, "
  f"ranked by number of distinct FN lines. Priority genes marked ★.\n")

t("""#table(
  columns: (0.35fr, 1.3fr, 1.3fr, 0.55fr, 2.2fr, 0.65fr, 0.65fr, 0.65fr),
  table.header(
    [*\\#*], [*Gene*], [*Locus*], [*Chr*], [*Description*],
    [*Del*], [*Dup*], [*Total*]
  ),""")
for rank, (_, row) in enumerate(gene_summary.iterrows(), 1):
    short = row["gene_name"].replace("Mt", "")
    star  = " ★" if any(p.replace("Mt", "") in short for p in PRIORITY_GENES) else ""
    desc  = (row["description"] or "—")[:55]
    t(f"  [{rank}], [{typ_escape(row['gene_name'])}{star}], "
      f"[{typ_escape(row['locus'])}], [{row['gene_chrom']}], "
      f"[{typ_escape(desc)}], [{int(row['n_lines_del'])}], "
      f"[{int(row['n_lines_dup'])}], [{int(row['n_lines_total'])}],")
t(")\n")

# Locus plots for all hit genes
t("= Locus Plots\n")
t("One plot per hit gene, drawn to genomic scale. Each bar is one FN indel call; "
  "y-axis labels give subline, zygosity, and indel size. The gene body is shown as a black bar "
  "beneath (minimum display width applied where gene is too small to see at the indel scale).\n")

for _, row in gene_summary.iterrows():
    locus    = row["locus"]
    plot_rel = gene_plot_paths[locus]
    short    = row["gene_name"].replace("Mt", "")
    star     = " ★" if any(p.replace("Mt", "") in short for p in PRIORITY_GENES) else ""
    caption  = (f"{typ_escape(row['gene_name'])}{star} · {typ_escape(locus)} · "
                f"{row['gene_chrom']} · {int(row['n_lines_total'])} FN line(s)")
    t(f"""#figure(
  image("{plot_rel}", width: 100%),
  caption: [{caption}],
)
""")

# Summary stat figures appendix
t("= Supplementary Figures\n")
for fname, caption in [
    ("fn_stats_dels_dups_per_line.png",
     "Deletions and duplications per FN line (top 30 by total count)."),
    ("fn_stats_del_size_dist.png",
     "Deletion size distribution (log scale) by zygosity."),
    ("fn_stats_dup_size_dist.png",
     "Duplication size distribution (log scale)."),
    ("fn_stats_del_size_by_chrom.png",
     "Deletion size distribution per chromosome (violin)."),
    ("fn_stats_genome_coverage.png",
     "Total sequence deleted and duplicated per chromosome (Mb)."),
]:
    t(f"""#figure(
  image("plots/{fname}", width: 100%),
  caption: [{typ_escape(caption)}],
)
""")

with open(TYP_PATH, "w") as fh:
    fh.write("\n".join(T) + "\n")
print(f"  Typst source written: {TYP_PATH}")

# ── Compile PDF ───────────────────────────────────────────────────────────────

print("Compiling PDF with typst...")
result = subprocess.run(
    ["typst", "compile", TYP_PATH, PDF_PATH],
    capture_output=True, text=True,
    cwd=ANALYSIS_DIR,  # so relative image paths resolve correctly
)
if result.returncode == 0:
    print(f"  PDF written: {PDF_PATH}")
else:
    print("  typst compilation failed:")
    print(result.stderr)
