#!/usr/bin/env python3
"""
FN deletion/duplication × nodulation ortholog intersection analysis
Pisum sativum JI2822
"""

import gzip
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path
from collections import defaultdict

# ── paths ──────────────────────────────────────────────────────────────────
BASE   = Path("/home/mbeavitt/Documents/ks_analysis/browser/ji2822")
FN_GFF = BASE / "fn_indels.gff3.gz"
OG_GFF = BASE / "ortholog_genes.gff3.gz"
ALIAS  = BASE / "chromAliases.txt"
OUTDIR = BASE          # plots saved here
REPORT = Path("/home/mbeavitt/Documents/ks_analysis/fn_indel_ortholog_report.md")

# ── priority genes (search these names) ────────────────────────────────────
PRIORITY_NAMES = ["LSH1", "LSH2", "NOOT1", "NOOT2",
                  "NFYa1", "NFYalpha1", "NF-YA1", "NFYA1"]
PRIORITY_DISPLAY = {
    "LSH1":  "LSH1",
    "LSH2":  "LSH2",
    "NOOT1": "NOOT1",
    "NOOT2": "NOOT2",
    "NFYa1": "NFYa1",
}

# ── seaborn style ───────────────────────────────────────────────────────────
sns.set_theme(style="whitegrid", font_scale=1.1)
PALETTE = sns.color_palette("colorblind")

# ═══════════════════════════════════════════════════════════════════════════
# 1.  HELPER: parse GFF3 attributes
# ═══════════════════════════════════════════════════════════════════════════
def parse_attrs(attr_str):
    d = {}
    for item in attr_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def read_gff3(path):
    rows = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            rows.append(parts)
    cols = ["seqid", "source", "type", "start", "end",
            "score", "strand", "phase", "attributes"]
    df = pd.DataFrame(rows, columns=cols)
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)
    return df


# ═══════════════════════════════════════════════════════════════════════════
# 2.  LOAD DATA
# ═══════════════════════════════════════════════════════════════════════════
print("Loading FN indels …")
fn_df = read_gff3(FN_GFF)

# Parse FN attributes
def parse_fn_attrs(row):
    d = parse_attrs(row["attributes"])
    return pd.Series({
        "indel_id":  d.get("ID", ""),
        "fn_line":   d.get("fn_line", ""),
        "fn_subline": d.get("fn_subline", ""),
        "zygosity":  d.get("zygosity", ""),
        "indel_type": row["type"],      # deletion / duplication
        "chrom":     row["seqid"],
        "start":     row["start"],
        "end":       row["end"],
        "size":      row["end"] - row["start"] + 1,
    })

fn = fn_df.apply(parse_fn_attrs, axis=1)
print(f"  FN records: {len(fn)} ({fn['indel_type'].value_counts().to_dict()})")

# ── chromosome alias map ─────────────────────────────────────────────────
alias_df = pd.read_csv(ALIAS, sep="\t", header=None,
                       names=["oz_id", "chr_name", "num"])
oz_to_chr = dict(zip(alias_df["oz_id"], alias_df["chr_name"]))
print("  Alias map:", oz_to_chr)

print("Loading ortholog genes …")
og_df = read_gff3(OG_GFF)

def parse_og_attrs(row):
    d = parse_attrs(row["attributes"])
    oz = row["seqid"]
    return pd.Series({
        "gene_id":     d.get("gene_id", d.get("ID", "").replace("gene:", "")),
        "name":        d.get("Name", ""),
        "alias":       d.get("Alias", ""),
        "medicago_symbol": d.get("medicago_symbol", ""),
        "pisum_symbol":    d.get("pisum_symbol", ""),
        "description": d.get("description", ""),
        "orthogroup":  d.get("orthogroup", ""),
        "oz_id":       oz,
        "chrom":       oz_to_chr.get(oz, oz),
        "start":       row["start"],
        "end":         row["end"],
        "strand":      row["strand"],
    })

og = og_df.apply(parse_og_attrs, axis=1)
print(f"  Ortholog genes: {len(og)}")


# ═══════════════════════════════════════════════════════════════════════════
# 3.  INTERSECTION
# ═══════════════════════════════════════════════════════════════════════════
print("Intersecting …")
hits = []

# Index orthologs by chromosome for speed
og_by_chr = defaultdict(list)
for _, g in og.iterrows():
    og_by_chr[g["chrom"]].append(g)

for _, indel in fn.iterrows():
    chrom = indel["chrom"]
    i_start, i_end = indel["start"], indel["end"]
    for gene in og_by_chr.get(chrom, []):
        # overlap: indel_start <= gene_end AND indel_end >= gene_start
        if i_start <= gene["end"] and i_end >= gene["start"]:
            hits.append({
                "gene_id":     gene["gene_id"],
                "gene_name":   gene["name"],
                "gene_alias":  gene["alias"],
                "medicago_symbol": gene["medicago_symbol"],
                "pisum_symbol":    gene["pisum_symbol"],
                "description": gene["description"],
                "gene_chrom":  gene["chrom"],
                "gene_start":  gene["start"],
                "gene_end":    gene["end"],
                "indel_id":    indel["indel_id"],
                "indel_type":  indel["indel_type"],
                "fn_line":     indel["fn_line"],
                "fn_subline":  indel["fn_subline"],
                "zygosity":    indel["zygosity"],
                "indel_start": i_start,
                "indel_end":   i_end,
                "indel_size_kb": (i_end - i_start + 1) / 1000,
            })

hits_df = pd.DataFrame(hits)
print(f"  Total hit records: {len(hits_df)}")
if len(hits_df) > 0:
    print(f"  Unique genes hit: {hits_df['gene_id'].nunique()}")
    print(f"  Unique FN lines with at least one hit: {hits_df['fn_line'].nunique()}")


# ═══════════════════════════════════════════════════════════════════════════
# 4.  PLOTS
# ═══════════════════════════════════════════════════════════════════════════
print("Generating plots …")

# ── 4a. Deletions and duplications per FN line (top 30) ─────────────────
per_line = fn.groupby(["fn_line", "indel_type"]).size().unstack(fill_value=0)
for col in ["deletion", "duplication"]:
    if col not in per_line.columns:
        per_line[col] = 0
per_line["total"] = per_line["deletion"] + per_line["duplication"]
top30_lines = per_line.nlargest(30, "total")

fig, ax = plt.subplots(figsize=(14, 6))
x = np.arange(len(top30_lines))
w = 0.45
ax.bar(x - w/2, top30_lines["deletion"],   width=w, label="Deletion",
       color=PALETTE[0], edgecolor="white", linewidth=0.5)
ax.bar(x + w/2, top30_lines["duplication"], width=w, label="Duplication",
       color=PALETTE[1], edgecolor="white", linewidth=0.5)
ax.set_xticks(x)
ax.set_xticklabels(top30_lines.index, rotation=45, ha="right", fontsize=9)
ax.set_xlabel("FN line")
ax.set_ylabel("Number of calls")
ax.set_title("Deletions and duplications per FN line (top 30 by total count)")
ax.legend(frameon=True)
fig.tight_layout()
fig.savefig(OUTDIR / "fn_stats_dels_dups_per_line.png", dpi=150)
plt.close(fig)
print("  Saved fn_stats_dels_dups_per_line.png")

# ── 4b. Deletion size distribution (hemizygous vs homozygous) ───────────
dels = fn[fn["indel_type"] == "deletion"].copy()
dels["size_kb"] = dels["size"] / 1000

fig, ax = plt.subplots(figsize=(10, 5))
for zyg, colour, label in [
    ("hemizygous",  PALETTE[0], "Hemizygous"),
    ("homozygous",  PALETTE[2], "Homozygous"),
    ("unknown",     PALETTE[3], "Unknown"),
]:
    sub = dels[dels["zygosity"] == zyg]["size_kb"]
    if len(sub) == 0:
        continue
    sns.histplot(sub, bins=60, kde=True, ax=ax, color=colour,
                 label=f"{label} (n={len(sub)})", alpha=0.55, stat="density",
                 edgecolor="none")
ax.set_xlabel("Deletion size (kb)")
ax.set_ylabel("Density")
ax.set_title("Distribution of deletion sizes by zygosity")
ax.legend(frameon=True)
# clip x-axis at 99th percentile for readability
clip = np.percentile(dels["size_kb"], 99)
ax.set_xlim(0, clip)
fig.tight_layout()
fig.savefig(OUTDIR / "fn_stats_del_size_dist.png", dpi=150)
plt.close(fig)
print("  Saved fn_stats_del_size_dist.png")

# ── 4c. Duplication size distribution ────────────────────────────────────
dups = fn[fn["indel_type"] == "duplication"].copy()
dups["size_kb"] = dups["size"] / 1000

fig, ax = plt.subplots(figsize=(10, 5))
sns.histplot(dups["size_kb"], bins=50, kde=True, ax=ax,
             color=PALETTE[1], alpha=0.7, stat="density", edgecolor="none",
             label=f"Duplications (n={len(dups)})")
ax.set_xlabel("Duplication size (kb)")
ax.set_ylabel("Density")
ax.set_title("Distribution of duplication sizes")
ax.legend(frameon=True)
clip_d = np.percentile(dups["size_kb"], 99)
ax.set_xlim(0, clip_d)
fig.tight_layout()
fig.savefig(OUTDIR / "fn_stats_dup_size_dist.png", dpi=150)
plt.close(fig)
print("  Saved fn_stats_dup_size_dist.png")

# ── 4d. Genome coverage per chromosome ────────────────────────────────────
fn["size_mb"] = fn["size"] / 1_000_000
cov = fn.groupby(["chrom", "indel_type"])["size_mb"].sum().unstack(fill_value=0)
for col in ["deletion", "duplication"]:
    if col not in cov.columns:
        cov[col] = 0
cov = cov.sort_index()

fig, ax = plt.subplots(figsize=(9, 5))
x = np.arange(len(cov))
w = 0.4
ax.bar(x - w/2, cov["deletion"],   width=w, label="Deleted",
       color=PALETTE[0], edgecolor="white")
ax.bar(x + w/2, cov["duplication"], width=w, label="Duplicated",
       color=PALETTE[1], edgecolor="white")
ax.set_xticks(x)
ax.set_xticklabels(cov.index, rotation=0)
ax.set_xlabel("Chromosome")
ax.set_ylabel("Total affected (Mb)")
ax.set_title("Genome coverage by FN deletions and duplications per chromosome")
ax.legend(frameon=True)
fig.tight_layout()
fig.savefig(OUTDIR / "fn_stats_genome_coverage.png", dpi=150)
plt.close(fig)
print("  Saved fn_stats_genome_coverage.png")

# ── 4e. Hits per gene (top 30) ────────────────────────────────────────────
if len(hits_df) > 0:
    gene_hits = (hits_df[hits_df["indel_type"] == "deletion"]
                 .groupby("gene_id")["fn_line"].nunique()
                 .rename("n_del_lines"))
    gene_hits_dup = (hits_df[hits_df["indel_type"] == "duplication"]
                     .groupby("gene_id")["fn_line"].nunique()
                     .rename("n_dup_lines"))
    gene_summary = (hits_df.groupby(["gene_id", "gene_name"])["fn_line"]
                    .nunique().reset_index()
                    .rename(columns={"fn_line": "n_total_lines"}))
    gene_summary = gene_summary.merge(gene_hits.reset_index(), on="gene_id", how="left")
    gene_summary = gene_summary.merge(gene_hits_dup.reset_index(), on="gene_id", how="left")
    gene_summary[["n_del_lines", "n_dup_lines"]] = gene_summary[
        ["n_del_lines", "n_dup_lines"]].fillna(0).astype(int)
    gene_summary = gene_summary.sort_values("n_total_lines", ascending=False)

    top30_genes = gene_summary.head(30)
    # use gene_name if available, otherwise gene_id
    labels = [n if n and n != gid else gid
              for n, gid in zip(top30_genes["gene_name"], top30_genes["gene_id"])]

    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(top30_genes))
    ax.bar(x - w/2, top30_genes["n_del_lines"], width=w,
           label="Deletion lines", color=PALETTE[0], edgecolor="white")
    ax.bar(x + w/2, top30_genes["n_dup_lines"], width=w,
           label="Duplication lines", color=PALETTE[1], edgecolor="white")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_xlabel("Ortholog gene")
    ax.set_ylabel("Number of distinct FN lines")
    ax.set_title("Top 30 ortholog genes by number of distinct FN lines overlapping")
    ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fn_stats_hits_per_gene.png", dpi=150)
    plt.close(fig)
    print("  Saved fn_stats_hits_per_gene.png")


# ═══════════════════════════════════════════════════════════════════════════
# 5.  SUMMARY STATISTICS
# ═══════════════════════════════════════════════════════════════════════════
dels_kb  = dels["size_kb"]
dups_kb  = dups["size_kb"]
total_dels  = len(dels)
total_dups  = len(dups)
total_lines = fn["fn_line"].nunique()
total_del_mb = dels["size"].sum() / 1_000_000
total_dup_mb = dups["size"].sum() / 1_000_000

n_genes_hit_del = hits_df[hits_df["indel_type"] == "deletion"]["gene_id"].nunique() if len(hits_df) else 0
n_genes_hit_dup = hits_df[hits_df["indel_type"] == "duplication"]["gene_id"].nunique() if len(hits_df) else 0

print(f"\nSummary statistics:")
print(f"  Total FN lines:        {total_lines}")
print(f"  Total deletions:       {total_dels}")
print(f"  Total duplications:    {total_dups}")
print(f"  Mean deletion size:    {dels_kb.mean():.1f} kb")
print(f"  Median deletion size:  {dels_kb.median():.1f} kb")
print(f"  Mean duplication size: {dups_kb.mean():.1f} kb")
print(f"  Median dup size:       {dups_kb.median():.1f} kb")
print(f"  Total Mb deleted:      {total_del_mb:.1f}")
print(f"  Genes hit by deletion: {n_genes_hit_del}")
print(f"  Genes hit by dup:      {n_genes_hit_dup}")


# ═══════════════════════════════════════════════════════════════════════════
# 6.  PRIORITY GENE LOOKUP
# ═══════════════════════════════════════════════════════════════════════════
print("\nLooking up priority genes …")

def find_priority_gene(og, search_names):
    """Return rows in ortholog df where Name or Alias matches any search_name."""
    results = []
    for sn in search_names:
        # case-insensitive search in Name and Alias columns
        mask = (og["name"].str.contains(sn, case=False, na=False) |
                og["alias"].str.contains(sn, case=False, na=False) |
                og["medicago_symbol"].str.contains(sn, case=False, na=False) |
                og["pisum_symbol"].str.contains(sn, case=False, na=False))
        results.append(og[mask])
    combined = pd.concat(results).drop_duplicates(subset=["gene_id"])
    return combined


# Search for each priority gene individually
priority_groups = {
    "LSH1":  ["LSH1"],
    "LSH2":  ["LSH2"],
    "NOOT1": ["NOOT1"],
    "NOOT2": ["NOOT2"],
    "NFYa1": ["NFYa1", "NFYalpha1", "NF-YA1", "NFYA1", "MtNF-YA1"],
}

priority_results = {}
for display, names in priority_groups.items():
    found = find_priority_gene(og, names)
    priority_results[display] = found
    print(f"  {display}: {len(found)} gene(s) found")
    if len(found) > 0:
        print(f"    -> {found[['gene_id','name','chrom','start','end']].to_string(index=False)}")


def get_hits_for_gene(gene_id, hits_df):
    if len(hits_df) == 0:
        return pd.DataFrame()
    return hits_df[hits_df["gene_id"] == gene_id].sort_values(["fn_line", "indel_type"])


# ═══════════════════════════════════════════════════════════════════════════
# 7.  WRITE MARKDOWN REPORT
# ═══════════════════════════════════════════════════════════════════════════
print("\nWriting report …")

def fmt_kb(n):
    return f"{n/1000:.1f} kb"

def md_table(df, columns=None, rename=None):
    """Return a markdown table string from a DataFrame."""
    if columns:
        df = df[columns]
    if rename:
        df = df.rename(columns=rename)
    lines = []
    lines.append("| " + " | ".join(str(c) for c in df.columns) + " |")
    lines.append("| " + " | ".join(["---"] * len(df.columns)) + " |")
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(v) for v in row.values) + " |")
    return "\n".join(lines)

lines = []

lines.append("# FN Deletion/Duplication × Nodulation Ortholog Intersection Report")
lines.append(f"\n*Generated: 2026-03-30 | Pisum sativum JI2822*\n")

lines.append("## Methods\n")
lines.append(
    "Fast neutron (FN) deletion and duplication calls were obtained from a genome-wide "
    "structural variant dataset for *Pisum sativum* JI2822 (`fn_indels.gff3.gz`; "
    f"{total_dels} deletions and {total_dups} duplications across {total_lines} FN lines). "
    "Nodulation ortholog genes were identified by synteny and orthology to *Medicago truncatula* "
    f"(`ortholog_genes.gff3.gz`; {len(og)} genes). "
    "Chromosome names in the ortholog file (OZ075428.1–OZ075434.1) were mapped to chr1–chr7 "
    "using a chromosome alias table. "
    "Interval intersection was performed using the criterion: "
    "*indel_start* ≤ *gene_end* AND *indel_end* ≥ *gene_start* (i.e., any overlap). "
    "For each gene, the number of distinct FN lines carrying a deletion or duplication overlapping "
    "that gene was counted separately; zygosity (hemizygous / homozygous / unknown) was recorded "
    "from the FN call attributes. All analyses were performed in Python (pandas, matplotlib, seaborn)."
)

# ── Priority genes ─────────────────────────────────────────────────────────
lines.append("\n---\n")
lines.append("## Priority Genes\n")
lines.append(
    "The following genes are of primary biological interest for nodulation in pea. "
    "Each is highlighted regardless of whether FN hits were found.\n"
)

for display, names in priority_groups.items():
    found_genes = priority_results[display]
    lines.append(f"### {display}\n")

    if len(found_genes) == 0:
        # NFYa1 special note
        if display == "NFYa1":
            lines.append(
                f"> **Note:** Searched under names: {', '.join(names)}. "
                "No gene matching any of these names was found in the ortholog gene set. "
                "MtNF-YA1 may not be represented in this ortholog collection, or may be "
                "annotated under a different identifier.\n"
            )
        else:
            lines.append(f"> **Note:** No gene matching '{display}' was found in the ortholog gene set.\n")
        continue

    for _, gene in found_genes.iterrows():
        size_kb = (gene["end"] - gene["start"] + 1) / 1000
        lines.append(f"**PSAT locus:** `{gene['gene_id']}`  ")
        lines.append(f"**Name:** {gene['name']}  ")
        lines.append(f"**Coordinates:** {gene['chrom']}:{gene['start']:,}–{gene['end']:,} ({size_kb:.1f} kb)  ")
        if gene["description"]:
            lines.append(f"**Description:** {gene['description']}  ")
        if gene["pisum_symbol"]:
            lines.append(f"**Pea symbol:** {gene['pisum_symbol']}  ")
        lines.append("")

        gene_hits = get_hits_for_gene(gene["gene_id"], hits_df)
        if len(gene_hits) == 0:
            lines.append("*No FN deletions or duplications overlap this gene.*\n")
        else:
            n_lines = gene_hits["fn_line"].nunique()
            n_del   = gene_hits[gene_hits["indel_type"] == "deletion"]["fn_line"].nunique()
            n_dup   = gene_hits[gene_hits["indel_type"] == "duplication"]["fn_line"].nunique()
            lines.append(
                f"**Overlapping FN lines:** {n_lines} total ({n_del} with deletion, {n_dup} with duplication)\n"
            )
            # one row per indel record
            tbl = gene_hits[["fn_subline", "indel_type", "zygosity",
                              "indel_start", "indel_end", "indel_size_kb"]].copy()
            tbl["indel_size_kb"] = tbl["indel_size_kb"].map(lambda v: f"{v:.1f}")
            tbl["interval"] = tbl.apply(
                lambda r: f"{r['indel_start']:,}–{r['indel_end']:,}", axis=1)
            tbl = tbl[["fn_subline", "indel_type", "zygosity", "interval", "indel_size_kb"]]
            tbl.columns = ["FN subline", "Type", "Zygosity", "Interval (bp)", "Size (kb)"]
            lines.append(md_table(tbl))
            lines.append("")

# ── Top candidates ─────────────────────────────────────────────────────────
lines.append("\n---\n")
lines.append("## Top Candidate Genes\n")

if len(hits_df) == 0:
    lines.append("*No ortholog genes were overlapped by any FN indel.*\n")
else:
    # rebuild full gene summary with descriptions
    gene_desc = og.set_index("gene_id")[["name", "chrom", "description",
                                          "pisum_symbol"]].drop_duplicates()
    gs = gene_summary.merge(gene_desc, on="gene_id", how="left")
    gs = gs.sort_values("n_total_lines", ascending=False).head(20)

    lines.append(
        "Ranked by number of distinct FN lines with any overlapping deletion or duplication "
        "(top 20 shown).\n"
    )
    tbl20 = gs[["gene_name", "gene_id", "chrom", "description",
                 "n_del_lines", "n_dup_lines", "n_total_lines"]].copy()
    tbl20.columns = ["Gene name", "PSAT locus", "Chr",
                     "Description", "Lines (del)", "Lines (dup)", "Lines (total)"]
    lines.append(md_table(tbl20))
    lines.append("")

# ── Summary statistics ─────────────────────────────────────────────────────
lines.append("\n---\n")
lines.append("## Summary Statistics\n")

lines.append("### FN library overview\n")
lines.append(f"| Metric | Value |")
lines.append(f"| --- | --- |")
lines.append(f"| Total FN lines | {total_lines} |")
lines.append(f"| Total deletion calls | {total_dels:,} |")
lines.append(f"| Total duplication calls | {total_dups:,} |")
lines.append(f"| Mean deletion size | {dels_kb.mean():.1f} kb |")
lines.append(f"| Median deletion size | {dels_kb.median():.1f} kb |")
lines.append(f"| Max deletion size | {dels_kb.max():.1f} kb |")
lines.append(f"| Mean duplication size | {dups_kb.mean():.1f} kb |")
lines.append(f"| Median duplication size | {dups_kb.median():.1f} kb |")
lines.append(f"| Max duplication size | {dups_kb.max():.1f} kb |")
lines.append(f"| Total genome deleted | {total_del_mb:.1f} Mb |")
lines.append(f"| Total genome duplicated | {total_dup_mb:.1f} Mb |")
lines.append("")

lines.append("### Ortholog gene intersection\n")
lines.append(f"| Metric | Value |")
lines.append(f"| --- | --- |")
lines.append(f"| Total ortholog genes analysed | {len(og)} |")
lines.append(f"| Ortholog genes hit by ≥1 deletion | {n_genes_hit_del} |")
lines.append(f"| Ortholog genes hit by ≥1 duplication | {n_genes_hit_dup} |")
lines.append(f"| Total ortholog genes hit (any) | {hits_df['gene_id'].nunique() if len(hits_df) else 0} |")
lines.append(f"| Total intersection records | {len(hits_df)} |")
lines.append("")

# Per-chromosome coverage table
cov_tbl = cov.copy()
cov_tbl.index.name = "Chromosome"
cov_tbl.columns = ["Deleted (Mb)", "Duplicated (Mb)"]
cov_tbl = cov_tbl.reset_index()
cov_tbl["Deleted (Mb)"]    = cov_tbl["Deleted (Mb)"].map(lambda v: f"{v:.2f}")
cov_tbl["Duplicated (Mb)"] = cov_tbl["Duplicated (Mb)"].map(lambda v: f"{v:.2f}")
lines.append("### Per-chromosome coverage\n")
lines.append(md_table(cov_tbl))
lines.append("")

# ── Plots reference ────────────────────────────────────────────────────────
lines.append("\n---\n")
lines.append("## Figures\n")
lines.append(
    "All figures are saved as PNG files in "
    f"`{OUTDIR}/`.\n"
)
plot_descriptions = [
    ("fn_stats_dels_dups_per_line.png",
     "Deletions and duplications per FN line (top 30 by total count)."),
    ("fn_stats_del_size_dist.png",
     "Distribution of deletion sizes (histogram + KDE), stratified by zygosity."),
    ("fn_stats_dup_size_dist.png",
     "Distribution of duplication sizes (histogram + KDE)."),
    ("fn_stats_genome_coverage.png",
     "Total Mb deleted and duplicated per chromosome."),
    ("fn_stats_hits_per_gene.png",
     "Top 30 nodulation ortholog genes ranked by number of distinct FN lines overlapping."),
]
for fname, desc in plot_descriptions:
    lines.append(f"- **{fname}** — {desc}")
lines.append("")

REPORT.write_text("\n".join(lines) + "\n")
print(f"Report written to {REPORT}")
print("Done.")
