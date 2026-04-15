#!/usr/bin/env python3
"""
kmer_clustering.py

Hierarchical clustering of pea accessions by k-mer presence/absence at
nodulation gene loci.

Modes:
  --genes [NAME ...]   Per-gene clustermaps for specific genes (default: key
                       nodulation genes LSH1, LSH2, NIN, NOOT1, NOOT2, NFYa1)
  --global             Global clustermap: accessions clustered on ALL variable
                       k-mer positions across every gene in the VCF, displayed
                       as a genes × accessions summary heatmap

Usage:
    python3 scripts/kmer_clustering.py                      # per-gene, default genes
    python3 scripts/kmer_clustering.py --genes LSH1 NIN     # per-gene, specific
    python3 scripts/kmer_clustering.py --global             # global clustering

Outputs:
    figures/kmer_clustering/{GeneName}.png
    figures/kmer_clustering/global_clustermap.png
"""

import argparse
import gzip
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import fcluster, linkage as sp_linkage

REPO          = Path(__file__).resolve().parent.parent
VCF           = Path("/Users/michaelbeavitt/Documents/motif_analysis/pea-browser/kmer_tracks/kmer_presence_sorted.vcf.gz")
GFF           = Path("/Users/michaelbeavitt/Documents/motif_analysis/pea-browser/ortholog_genes.gff3.gz")
OUT           = REPO / "figures" / "kmer_clustering"
SPECIES_CSV   = REPO / "scripts" / "accession_species.csv"

SPECIES_PALETTE = {
    "Pisum sativum":                "#4e79a7",
    "Pisum sativum subsp. elatius": "#f28e2b",
    "Pisum fulvum":                 "#e15759",
    "Pisum abyssinicum":            "#76b7b2",
}
SPECIES_LABELS = {
    "Pisum sativum":                r"$\it{P.\ sativum}$",
    "Pisum sativum subsp. elatius": r"$\it{P.\ sativum}$ subsp. $\it{elatius}$",
    "Pisum fulvum":                 r"$\it{P.\ fulvum}$",
    "Pisum abyssinicum":            r"$\it{P.\ abyssinicum}$",
}

# Feature values in the VCF are numbered (exon_1, intron_2 etc.)
# Normalise by stripping the trailing _N before colour lookup
FEATURE_COLOURS = {
    "promoter":   "#4e79a7",
    "exon":       "#f28e2b",
    "intron":     "#bab0ac",
    "downstream": "#59a14f",
    "UTR":        "#76b7b2",
    "CDS":        "#f28e2b",
}
DEFAULT_FEAT_COL = "#e0e0e0"

def feat_colour(raw_feat):
    base = re.sub(r"_\d+$", "", raw_feat)   # exon_1 → exon, intron_2 → intron
    return FEATURE_COLOURS.get(base, DEFAULT_FEAT_COL)

DEFAULT_GENES = ["LSH1", "LSH2", "NIN", "NOOT1", "NOOT2", "NF-YA1"]


# ── Gene name lookup ─────────────────────────────────────────────────────────

def load_gene_map(gff_path):
    """Returns ({locus_id: name}, {name/alias: locus_id})."""
    locus_to_name, name_to_locus = {}, {}
    with gzip.open(gff_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            attrs = parts[8]
            gid = re.search(r"gene_id=([^;]+)", attrs)
            nam = re.search(r"Name=([^;]+)", attrs)
            if gid and nam:
                locus = gid.group(1).strip()
                name  = nam.group(1).strip()
                locus_to_name[locus] = name
                name_to_locus[name] = locus
                short = re.sub(r"^Mt", "", name)
                if short not in name_to_locus:
                    name_to_locus[short] = locus
    return locus_to_name, name_to_locus


def load_species(csv_path):
    """Returns Series: sample_name → species string."""
    df = pd.read_csv(csv_path)
    return df.set_index("sample")["species"]


def resolve_genes(gene_list, name_to_locus):
    """Resolve gene name strings to locus IDs."""
    target_loci = []
    for g in gene_list:
        locus = name_to_locus.get(g) or name_to_locus.get("Mt" + g)
        if locus:
            target_loci.append(locus)
        else:
            print(f"  WARNING: '{g}' not found in GFF — skipping", file=sys.stderr)
    return target_loci


# ── VCF streaming ────────────────────────────────────────────────────────────

def _parse_line(fields, samples):
    """Parse one VCF data line. Returns (locus, relpos, feature, values_array)."""
    info   = fields[7]
    gm     = re.search(r"GENE=([^;]+)", info)
    rm     = re.search(r"RELPOS=(\d+)", info)
    fm     = re.search(r"FEATURE=([^;]+)", info)
    locus  = gm.group(1) if gm else None
    relpos = int(rm.group(1)) if rm else 0
    feat   = fm.group(1) if fm else "unknown"
    gts    = fields[9:9 + len(samples)]
    values = np.array([1 if g.startswith("1") else 0 for g in gts], dtype=np.uint8)
    return locus, relpos, feat, values


def stream_vcf(vcf_path, target_loci=None):
    """
    Yield (locus_id, relpos_index, feat_series, matrix) per gene.
    matrix shape: (n_positions, n_samples), dtype uint8.
    If target_loci is None, yield all genes.
    """
    target_set = set(target_loci) if target_loci is not None else None
    samples = None
    current_locus = None
    rows, relpos_list, feat_list = [], [], []

    def emit():
        if not rows:
            return None
        return np.stack(rows), relpos_list[:], feat_list[:]

    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]
                continue

            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue

            locus, relpos, feat, values = _parse_line(fields, samples)
            if locus is None:
                continue
            if target_set is not None and locus not in target_set:
                if current_locus != locus:
                    # Flush if we were accumulating a different gene
                    if current_locus in (target_set or set()) and rows:
                        yield current_locus, *emit()
                    current_locus = locus
                    rows, relpos_list, feat_list = [], [], []
                continue

            if locus != current_locus:
                if current_locus is not None and rows:
                    if target_set is None or current_locus in target_set:
                        yield current_locus, *emit()
                current_locus = locus
                rows, relpos_list, feat_list = [], [], []

            rows.append(values)
            relpos_list.append(relpos)
            feat_list.append(feat)

    if current_locus is not None and rows:
        if target_set is None or current_locus in target_set:
            yield current_locus, *emit()

    return samples   # not actually returned via yield, accessed via closure below


def get_samples(vcf_path):
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                return line.strip().split("\t")[9:]
    return []


# ── Per-gene clustermap ───────────────────────────────────────────────────────

def make_gene_clustermap(locus_id, gene_name, matrix, relpos_list, feat_list,
                         samples, out_path, n_clusters=None):
    """
    matrix: (n_positions, n_samples) uint8
    Clusters accessions; keeps positions in RELPOS order.
    Returns cluster-label Series (accession → cluster int) or None.
    """
    import matplotlib.patches as mpatches

    mat_df = pd.DataFrame(matrix, index=relpos_list, columns=samples)

    # Variable positions only
    variable = mat_df[mat_df.std(axis=1) > 0]
    feat_var  = pd.Series(feat_list, index=relpos_list)[variable.index]

    n_var, n_acc = len(variable), len(samples)
    print(f"  {gene_name}: {n_var} variable / {len(mat_df)} positions, {n_acc} accessions")

    if n_var < 10:
        print(f"    → skipping (< 10 variable positions)")
        return None

    heatmap = variable.T   # accessions × positions

    # Pass feature colours to reserve space; we redraw the axes as a gene model
    feat_colours_series = feat_var.map(feat_colour)
    feat_df = feat_colours_series.to_frame(name=gene_name)

    fig_h = max(8, n_acc * 0.012 + 5)

    g = sns.clustermap(
        heatmap,
        method="ward",
        metric="euclidean",
        col_cluster=False,
        cmap="Blues", vmin=0, vmax=1,
        xticklabels=False, yticklabels=False,
        figsize=(min(20, n_var * 0.015 + 10), min(18, fig_h)),
        col_colors=feat_df,
        cbar_pos=None,
        dendrogram_ratio=(0.15, 0.0),
        colors_ratio=0.14,   # enlarged to accommodate gene model
        linewidths=0,
    )
    g.ax_heatmap.set_xlabel(
        f"k-mer position in gene window  (n = {n_var} variable positions)",
        labelpad=8, fontsize=15,
    )
    g.ax_heatmap.set_ylabel(f"Accessions  (n = {n_acc})", labelpad=8, fontsize=15)
    g.figure.suptitle(f"{gene_name}  ·  {locus_id}",
                      fontsize=18, fontweight="bold", y=1.01)

    # ── Redraw col_colors axis as a genome-browser gene model ────────────────────
    # Heights (ymin, height) per feature type — exons tallest, introns thin
    FEATURE_HEIGHTS = {
        "promoter":   (0.22, 0.56),
        "exon":       (0.05, 0.90),
        "CDS":        (0.05, 0.90),
        "UTR":        (0.20, 0.60),
        "intron":     (0.46, 0.08),
        "downstream": (0.22, 0.56),
    }
    DEFAULT_H = (0.25, 0.50)

    ax_model = g.ax_col_colors
    ax_model.cla()

    xlim   = g.ax_heatmap.get_xlim()
    n_cols = heatmap.shape[1]
    x_scale = (xlim[1] - xlim[0]) / n_cols

    ax_model.set_xlim(*xlim)
    # bottom padding (-0.55 → 0) acts as breathing room before the heatmap below
    ax_model.set_ylim(-0.55, 1.05)
    ax_model.set_yticks([])
    ax_model.set_xticks([])
    for sp in ax_model.spines.values():
        sp.set_visible(False)

    # Thin backbone through entire window
    ax_model.axhline(0.50, color="#aaaaaa", lw=0.8, zorder=1)

    # Group consecutive positions by base feature and draw blocks of varying height
    base_feats = [re.sub(r"_\d+$", "", f) for f in feat_var.values]
    i = 0
    while i < len(base_feats):
        j = i + 1
        while j < len(base_feats) and base_feats[j] == base_feats[i]:
            j += 1
        feat  = base_feats[i]
        ymin, ht = FEATURE_HEIGHTS.get(feat, DEFAULT_H)
        color = FEATURE_COLOURS.get(feat, DEFAULT_FEAT_COL)
        ax_model.add_patch(
            mpatches.Rectangle(
                (xlim[0] + i * x_scale, ymin), (j - i) * x_scale, ht,
                fc=color, ec="none", zorder=2, alpha=0.9,
            )
        )
        i = j

    # ── Legends: both placed outside right, anchored to heatmap ─────────────────
    present_feats = {re.sub(r"_\d+$", "", f) for f in feat_list}
    feat_order = ["promoter", "exon", "CDS", "UTR", "intron", "downstream"]
    feat_handles = [
        mpatches.Patch(fc=FEATURE_COLOURS[f], label=f.capitalize())
        for f in feat_order if f in present_feats
    ]
    extra_artists = []
    if feat_handles:
        feat_legend = g.ax_heatmap.legend(
            handles=feat_handles, title="Feature", title_fontsize=12,
            loc="lower left", bbox_to_anchor=(1.02, 1.02),
            fontsize=12, frameon=True,
            edgecolor="#cccccc", facecolor="white", framealpha=0.9,
        )
        g.ax_heatmap.add_artist(feat_legend)
        extra_artists.append(feat_legend)

    pa_handles = [
        mpatches.Patch(fc="#084594", label="Present"),
        mpatches.Patch(fc="#f7fbff", ec="#cccccc", lw=0.5, label="Absent"),
    ]
    pa_legend = g.ax_heatmap.legend(
        handles=pa_handles, title="k-mer", title_fontsize=12,
        loc="upper left", bbox_to_anchor=(1.02, 1.0),
        fontsize=12, frameon=True,
        edgecolor="#cccccc", facecolor="white", framealpha=0.9,
    )
    extra_artists.append(pa_legend)

    g.figure.savefig(out_path, dpi=150, bbox_inches="tight",
                     bbox_extra_artists=extra_artists)
    plt.close(g.figure)
    print(f"    → {out_path.name}")

    # Extract cluster labels from row linkage
    if n_clusters and hasattr(g, "dendrogram_row") and g.dendrogram_row is not None:
        lnk = g.dendrogram_row.linkage
        labels = fcluster(lnk, n_clusters, criterion="maxclust")
        return pd.Series(labels, index=heatmap.index, name=gene_name)
    return None


# ── Global clustermap ─────────────────────────────────────────────────────────

def make_global_clustermap(gene_summaries, locus_to_name, samples,
                           core_loci, out_path, dot, n_pos, n_clusters=10,
                           species_series=None):
    """
    Cluster accessions using exact pairwise Euclidean distances computed from
    the accumulated dot product matrix (no compression).
    Display heatmap uses mean kmer presence per gene (for visualisation only).
    """
    from scipy.spatial.distance import squareform

    print(f"  Global: computing pairwise distances from {n_pos:,} variable positions "
          f"× {len(samples)} accessions")

    # Exact pairwise squared Euclidean from dot product:
    # ||x_i - x_j||^2 = ||x_i||^2 + ||x_j||^2 - 2 x_i·x_j
    norms   = np.diag(dot)
    D_sq    = norms[:, None] + norms[None, :] - 2 * dot
    D_sq    = np.clip(D_sq, 0, None)   # guard against tiny negatives from float32
    D       = np.sqrt(D_sq)
    lnk     = sp_linkage(squareform(D, checks=False), method="ward")
    cluster_labels = fcluster(lnk, n_clusters, criterion="maxclust")
    cluster_series = pd.Series(cluster_labels, index=samples, name="global_cluster")

    # Save cluster assignments
    csv_path = out_path.parent / "clusters_global.csv"
    cluster_series.to_csv(csv_path, header=True)
    print(f"    → {csv_path.name}  ({n_clusters} clusters)")

    from scipy.cluster.hierarchy import leaves_list, dendrogram as sp_dendrogram
    from matplotlib.gridspec import GridSpec

    col_order       = leaves_list(lnk)
    palette         = sns.color_palette("tab10", n_clusters)
    cluster_ordered = cluster_series.iloc[col_order].values

    # Display matrix: mean kmer presence per gene (heatmap only, not used for clustering)
    gene_labels = [locus_to_name.get(l, l) for l in gene_summaries]
    all_data    = np.stack(list(gene_summaries.values()))
    all_df      = pd.DataFrame(all_data, index=gene_labels, columns=samples)
    all_df      = all_df[all_df.std(axis=1) > 0.01]

    display_names = [locus_to_name.get(l, l) for l in core_loci if l in gene_summaries]
    display_df    = all_df.loc[all_df.index.isin(display_names)]
    if display_df.empty:
        display_df = all_df

    # Sort rows by mean kmer presence across all accessions (ascending → most variable at top)
    display_cols    = display_df.iloc[:, col_order]
    row_means       = display_cols.mean(axis=1)
    display_ordered = display_cols.iloc[row_means.argsort()]

    # Labels: always only the DEFAULT_GENES core set, regardless of what is displayed
    core_name_set = set(DEFAULT_GENES) | {f"Mt{g}" for g in DEFAULT_GENES}

    # ── Layout: 4 rows × 2 cols; col 1 is narrow colorbar column ──────────────
    n_display = len(display_ordered)
    heatmap_h = max(4, min(12, n_display * 0.18))
    fig = plt.figure(figsize=(18, 4 + heatmap_h))
    gs  = GridSpec(4, 2, figure=fig,
                   height_ratios=[2.5, 0.25, 0.25, heatmap_h],
                   width_ratios=[1, 0.025],
                   hspace=0.04, wspace=0.02)

    ax_dendro  = fig.add_subplot(gs[0, 0])
    ax_bar     = fig.add_subplot(gs[1, 0])
    ax_species = fig.add_subplot(gs[2, 0])
    ax_heat    = fig.add_subplot(gs[3, 0])
    ax_cbar    = fig.add_subplot(gs[3, 1])
    # placeholders keep widths consistent
    for ax in (fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1]),
               fig.add_subplot(gs[2, 1])):
        ax.axis("off")

    # ── Dendrogram ─────────────────────────────────────────────────────────────
    sp_dendrogram(lnk, ax=ax_dendro, no_labels=True, color_threshold=0,
                  link_color_func=lambda _: "#888888")
    # scipy leaf positions: 5, 15, 25, … → range (0, n*10)
    x_min, x_max = ax_dendro.get_xlim()
    ax_dendro.set_xlim(x_min, x_max)
    ax_dendro.axis("off")
    ax_dendro.set_title(
        "Global k-mer diversity",
        fontsize=15, fontweight="bold",
    )

    # ── Cluster colour bar (same x-range as dendrogram) ────────────────────────
    for i, clust in enumerate(cluster_ordered):
        x_lo = x_min + (x_max - x_min) * i       / len(samples)
        x_hi = x_min + (x_max - x_min) * (i + 1) / len(samples)
        ax_bar.axvspan(x_lo, x_hi, color=palette[clust - 1], lw=0)
    ax_bar.set_xlim(x_min, x_max)
    ax_bar.set_yticks([])
    ax_bar.set_xticks([])
    ax_bar.set_ylabel("Cluster", fontsize=12, rotation=0, labelpad=38, va="center")

    # ── Species bar ────────────────────────────────────────────────────────────
    if species_series is not None:
        species_ordered = species_series.reindex(samples).iloc[col_order]
        for i, sp in enumerate(species_ordered):
            x_lo = x_min + (x_max - x_min) * i       / len(samples)
            x_hi = x_min + (x_max - x_min) * (i + 1) / len(samples)
            color = SPECIES_PALETTE.get(sp, "#cccccc") if pd.notna(sp) else "#cccccc"
            ax_species.axvspan(x_lo, x_hi, color=color, lw=0)
    ax_species.set_xlim(x_min, x_max)
    ax_species.set_yticks([])
    ax_species.set_xticks([])
    ax_species.set_ylabel("Species", fontsize=12, rotation=0, labelpad=38, va="center")

    # ── Heatmap (extent matches dendrogram x-range) ────────────────────────────
    n_genes = len(display_ordered)
    im = ax_heat.imshow(
        display_ordered.values, aspect="auto", cmap="Blues",
        vmin=0, vmax=1, interpolation="nearest",
        extent=[x_min, x_max, n_genes - 0.5, -0.5],
    )
    ax_heat.set_xlim(x_min, x_max)
    ax_heat.set_ylim(n_genes - 0.5, -0.5)
    ax_heat.set_yticks([])

    # Annotate core gene rows with pointer lines; spread collisions with curved arcs
    labeled  = [(i, name) for i, name in enumerate(display_ordered.index)
                if name in core_name_set]
    if labeled:
        min_sep = 1.5   # minimum label separation in data (row) units
        data_ys = [float(y) for y, _ in labeled]
        text_ys = list(data_ys)
        for _ in range(200):
            moved = False
            for i in range(len(text_ys) - 1):
                if text_ys[i + 1] - text_ys[i] < min_sep:
                    mid = (text_ys[i] + text_ys[i + 1]) / 2
                    text_ys[i]     = mid - min_sep / 2
                    text_ys[i + 1] = mid + min_sep / 2
                    moved = True
            if not moved:
                break
        for (data_y, name), text_y in zip(labeled, text_ys):
            offset = text_y - data_y
            rad    = -0.35 * np.sign(offset) if abs(offset) > 0.1 else 0.0
            ax_heat.annotate(
                name,
                xy=(0.0, data_y),
                xycoords=("axes fraction", "data"),
                xytext=(-0.01, text_y),
                textcoords=("axes fraction", "data"),
                fontsize=13, fontweight="bold",
                ha="right", va="center",
                clip_on=False, annotation_clip=False,
                arrowprops=dict(
                    arrowstyle="-",
                    connectionstyle=f"arc3,rad={rad}",
                    color="#555555", lw=0.8, clip_on=False,
                ),
            )
    ax_heat.set_xticks([])
    ax_heat.set_xlabel(f"Accessions  (n = {len(samples)}, ordered by clustering)",
                       labelpad=6, fontsize=12)

    plt.colorbar(im, cax=ax_cbar, label="Mean k-mer presence")
    ax_cbar.yaxis.set_tick_params(labelsize=11)
    ax_cbar.set_ylabel("Mean k-mer presence", fontsize=11)

    # ── Cluster legend (upper left) ────────────────────────────────────────────
    cluster_handles = [plt.Rectangle((0, 0), 1, 1, fc=palette[i], label=f"Cluster {i+1}")
                       for i in range(n_clusters)]
    cluster_legend = ax_dendro.legend(
        handles=cluster_handles, loc="upper left", bbox_to_anchor=(0.0, 1.0),
        fontsize=11, ncol=2, frameon=True,
        title="Cluster", title_fontsize=11,
        edgecolor="#cccccc", facecolor="white", framealpha=0.9)
    ax_dendro.add_artist(cluster_legend)

    # ── Species legend (upper right) ───────────────────────────────────────────
    if species_series is not None:
        present = [sp for sp in SPECIES_PALETTE if sp in species_series.values]
        sp_handles = [plt.Rectangle((0, 0), 1, 1,
                                    fc=SPECIES_PALETTE[sp],
                                    label=SPECIES_LABELS.get(sp, sp))
                      for sp in present]
        ax_dendro.legend(
            handles=sp_handles, loc="upper right", bbox_to_anchor=(1.0, 1.0),
            fontsize=11, ncol=1, frameon=True,
            title="Species", title_fontsize=11,
            edgecolor="#cccccc", facecolor="white", framealpha=0.9)

    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"    → {out_path.name}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--genes", nargs="*", metavar="GENE",
                        help="Gene names for per-gene plots (default: key nodulation genes)")
    parser.add_argument("--global", dest="do_global", action="store_true",
                        help="Also produce global clustermap across all orthologs")
    parser.add_argument("--clusters", type=int, default=10, metavar="N",
                        help="Number of clusters to extract (default: 10)")
    parser.add_argument("--show-all-genes", dest="show_all_genes", action="store_true",
                        help="Show all ortholog genes in the global heatmap (not just core 6)")
    parser.add_argument("--core-only", dest="core_only", action="store_true",
                        help="Cluster on core genes only (not all orthologs)")
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)

    print("Loading gene name map …")
    locus_to_name, name_to_locus = load_gene_map(GFF)

    samples = get_samples(VCF)
    print(f"Samples in VCF: {len(samples)}")

    species_series = None
    if SPECIES_CSV.exists():
        species_series = load_species(SPECIES_CSV)
        print(f"Loaded species for {len(species_series)} accessions")

    # Resolve target gene loci for per-gene plots
    gene_list   = args.genes if args.genes is not None else DEFAULT_GENES
    target_loci = resolve_genes(gene_list, name_to_locus)

    if args.do_global and not args.core_only:
        all_loci = None   # stream all genes
        print(f"\nMode: per-gene ({len(target_loci)} genes) + global (all genes)")
    elif args.do_global and args.core_only:
        all_loci = target_loci   # stream core genes only
        print(f"\nMode: per-gene ({len(target_loci)} genes) + global (core genes only)")
    else:
        all_loci = target_loci
        print(f"\nMode: per-gene only ({len(target_loci)} genes)")

    per_gene_set = set(target_loci)
    gene_summaries = {}   # for display heatmap: {locus_id: mean_presence array}
    n_acc = len(samples)
    global_dot   = np.zeros((n_acc, n_acc), dtype=np.float64)
    global_n_pos = 0

    print(f"Streaming VCF …")
    for locus_id, matrix, relpos_list, feat_list in stream_vcf(VCF, target_loci=all_loci):
        gene_name = locus_to_name.get(locus_id, locus_id)

        # Per-gene plot
        if locus_id in per_gene_set:
            out_path = OUT / f"{gene_name.replace('/', '_').replace(' ', '_')}.png"
            make_gene_clustermap(locus_id, gene_name, matrix, relpos_list,
                                 feat_list, samples, out_path,
                                 n_clusters=args.clusters)

        # Accumulate for global clustering
        if args.do_global:
            mat_df   = pd.DataFrame(matrix, columns=samples)
            variable = mat_df[mat_df.std(axis=1) > 0]
            if len(variable) >= 1:
                # Display summary: mean presence per accession per gene
                gene_summaries[locus_id] = variable.mean(axis=0).values
                # Clustering: accumulate dot product incrementally (no compression)
                X = variable.values.T.astype(np.float32)   # (n_acc, n_var_pos)
                global_dot   += X @ X.T
                global_n_pos += X.shape[1]

    if args.do_global and gene_summaries:
        print(f"\nBuilding global clustermap ({len(gene_summaries)} genes) …")
        if args.show_all_genes:
            display_loci = list(gene_summaries.keys())
            out_name = "global_clustermap_full.png"
        else:
            display_loci = target_loci
            out_name = "global_clustermap.png"
        make_global_clustermap(
            gene_summaries, locus_to_name, samples,
            core_loci=display_loci,
            out_path=OUT / out_name,
            dot=global_dot,
            n_pos=global_n_pos,
            n_clusters=args.clusters,
            species_series=species_series,
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
