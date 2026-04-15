# FN Indel × Nodulation Ortholog Intersection Report

## Methods

Fast neutron (FN) deletion and duplication calls (1603 deletions, 451 duplications across 197 FN lines) were intersected with 488 nodulation ortholog gene models mapped to the *Pisum sativum* JI2822 v1.3 assembly. FN indel coordinates are from skim-sequencing coverage analysis (Supplementary Tables 5–6). Ortholog coordinates were mapped from *Medicago truncatula* A17 via synteny. An overlap is defined as any FN interval sharing ≥1 bp with a gene model. Hemizygous and homozygous deletion calls are treated independently (a single subline may contribute both for the same locus, representing different calling stringencies over the same coverage drop). Chromosome aliases (chr1–7 ↔ OZ075428–34.1) were resolved using the JBrowse chromAliases.txt.

## Summary Statistics

| Metric | Value |
|--------|-------|
| FN lines | 197 |
| Total deletions | 1603 |
| Total duplications | 451 |
| Mean deletion size | 592 kb |
| Median deletion size | 230 kb |
| Mean duplication size | 446 kb |
| Median duplication size | 250 kb |
| Total Mb deleted (genome-wide) | 949.0 Mb |
| Total Mb duplicated (genome-wide) | 201.4 Mb |
| Ortholog genes hit by ≥1 deletion | 40 |
| Ortholog genes hit by ≥1 duplication | 15 |
| Ortholog genes hit by any indel | 55 / 488 |
| Ortholog genes with no hit | 433 |

## Priority Genes

### LSH1

**Locus:** PSAT_LOCUS25697  
**Position:** chr6:69,313,200–69,314,098 (0.9 kb)  
**Description:** protein LIGHT-DEPENDENT SHORT HYPOCOTYLS 4-like

**2 hit(s) across 1 FN line(s):**

| FN subline | Type | Zygosity | Interval | Size (kb) |
|------------|------|----------|----------|-----------|
| FN2026/7 | deletion | hemizygous | chr6:68,280,000–70,230,000 | 1950 |
| FN2026/7 | deletion | homozygous | chr6:68,300,000–70,230,000 | 1930 |

### LSH2

**Locus:** PSAT_LOCUS13554  
**Position:** chr3:370,690,287–370,691,234 (0.9 kb)  
**Description:** protein LIGHT-DEPENDENT SHORT HYPOCOTYLS 4-like

*No FN indels overlap this gene.*

### NOOT1

**Locus:** PSAT_LOCUS13168  
**Position:** chr3:317,691,479–317,694,007 (2.5 kb)  
**Description:** NODULE ROOT 1

**2 hit(s) across 1 FN line(s):**

| FN subline | Type | Zygosity | Interval | Size (kb) |
|------------|------|----------|----------|-----------|
| FN3185/1325 | deletion | hemizygous | chr3:317,250,000–320,440,000 | 3190 |
| FN3185/1325 | deletion | homozygous | chr3:317,250,000–320,440,000 | 3190 |

### NOOT2

**Locus:** PSAT_LOCUS26318  
**Position:** chr6:184,006,759–184,009,064 (2.3 kb)  
**Description:** NODULE ROOT 2

**2 hit(s) across 1 FN line(s):**

| FN subline | Type | Zygosity | Interval | Size (kb) |
|------------|------|----------|----------|-----------|
| FN2087/2 | deletion | hemizygous | chr6:183,530,000–185,600,000 | 2070 |
| FN2087/2 | deletion | homozygous | chr6:183,530,000–185,600,000 | 2070 |

### NF-YA1

**Locus:** PSAT_LOCUS26253  
**Position:** chr6:161,422,771–161,424,650 (1.9 kb)  
**Description:** nuclear factor Y subunit A1

*No FN indels overlap this gene.*

## Top Candidate Genes

Ranked by number of distinct FN lines with an overlapping deletion or duplication. Priority genes marked ★.

| Rank | Gene | Locus | Chr | Description | Lines w/ del | Lines w/ dup | Total lines |
|------|------|-------|-----|-------------|-------------|-------------|-------------|
| 1 | **PSAT_LOCUS12094** | PSAT_LOCUS12094 | chr3 | — | 4 | 0 | 4 |
| 2 | **MtFLS1** | PSAT_LOCUS9418 | chr2 | flavonol synthase | 2 | 0 | 2 |
| 3 | **MtGH3** | PSAT_LOCUS9284 | chr2 | Gretchen Hagen3 (auxin-responsive  GH3 family) | 2 | 0 | 2 |
| 4 | **MtCKX4** | PSAT_LOCUS13232 | chr3 | cytokinin oxydase 4 | 2 | 0 | 2 |
| 5 | **MtCASTOR** | PSAT_LOCUS14811 | chr3 | Lotus japonicus CASTOR ion  channel-like | 2 | 0 | 2 |
| 6 | **FUN** | PSAT_LOCUS9316 | chr2 | Putative transcription factor bZIP family | 1 | 0 | 1 |
| 7 | **IAMT1** | PSAT_LOCUS21425 | chr5 | Putative indole-3-acetate O-methyltransferase | 1 | 0 | 1 |
| 8 | **LOG-like|LOG1, cytokinin riboside 5'-monophosphate phosphoribohydrolase** | PSAT_LOCUS13817 | chr3 | — | 1 | 0 | 1 |
| 9 | **MEX1Â ** | PSAT_LOCUS30328 | chr7 | — | 0 | 1 | 1 |
| 10 | **MtBHLHB** | PSAT_LOCUS31164 | chr7 | bHLH Transcription factor | 1 | 0 | 1 |
| 11 | **MtARF4a** | PSAT_LOCUS29732 | chr7 | — | 1 | 0 | 1 |
| 12 | **Mt3CHA** | PSAT_LOCUS33237 | chr7 | CCCH-type zinc finger | 0 | 1 | 1 |
| 13 | **Mt12b** | PSAT_LOCUS12620 | chr3 | — | 1 | 0 | 1 |
| 14 | **MtCKX3** | PSAT_LOCUS5090 | chr1 | cytokinin oxidase 3 | 1 | 0 | 1 |
| 15 | **MtCKX6** | PSAT_LOCUS20805 | chr5 | cytokinin oxydase 6 | 1 | 0 | 1 |
| 16 | **MtCAMTA1** | PSAT_LOCUS3780 | chr1 | Calmodulin-binding transcription activator 1 | 0 | 1 | 1 |
| 17 | **MtBRI** | PSAT_LOCUS23282 | chr5 | Brassinosteroid Insensitive 1 (BR receptor) | 1 | 0 | 1 |
| 18 | **MtDGK** | PSAT_LOCUS16447 | chr4 | — | 0 | 1 | 1 |
| 19 | **MtCYCA1** | PSAT_LOCUS9381 | chr2 | — | 1 | 0 | 1 |
| 20 | **MtEXP7** | PSAT_LOCUS7810 | chr2 | — | 1 | 0 | 1 |
| 21 | **MtDMI1** | PSAT_LOCUS3652 | chr1 | does not make infection 1|ion  channel | 0 | 1 | 1 |
| 22 | **MtFPN2** | PSAT_LOCUS29222 | chr7 | ferroportin 2 | 1 | 0 | 1 |
| 23 | **MtGRF1** | PSAT_LOCUS30323 | chr7 | growth-regulating factor 1 (plant-specific transcription fac | 0 | 1 | 1 |
| 24 | **MtHOM_TALE_A ** | PSAT_LOCUS17159 | chr4 | Putative transcription factor Homeodomain-TALE-BEL family | 0 | 1 | 1 |
| 25 | **Epidermal Patterning Factor-Like** | PSAT_LOCUS24205 | chr5 | signalling peptide | 1 | 0 | 1 |
| 26 | **MtHPT4** | PSAT_LOCUS19295 | chr5 | Histidine PhosphoTransfer protein 4 | 1 | 0 | 1 |
| 27 | **MtIPT1** | PSAT_LOCUS28589 | chr6 | isopentenyl transferase | 1 | 0 | 1 |
| 28 | **MtLSH1** ★ | PSAT_LOCUS25697 | chr6 | protein LIGHT-DEPENDENT SHORT HYPOCOTYLS 4-like | 1 | 0 | 1 |
| 29 | **MtLINL** | PSAT_LOCUS33304 | chr7 | — | 0 | 1 | 1 |
| 30 | **MtMATE58** | PSAT_LOCUS28149 | chr6 | — | 1 | 0 | 1 |

## Figures

![Stacked bar: deletions and duplications per FN line (top 30 by total)](fn_stats_dels_dups_per_line.png)
*Stacked bar: deletions and duplications per FN line (top 30 by total)*

![Deletion size distribution, hemizygous vs homozygous](fn_stats_del_size_dist.png)
*Deletion size distribution, hemizygous vs homozygous*

![Duplication size distribution](fn_stats_dup_size_dist.png)
*Duplication size distribution*

![Total Mb deleted / duplicated per chromosome](fn_stats_genome_coverage.png)
*Total Mb deleted / duplicated per chromosome*

![Top 30 ortholog genes by number of overlapping FN lines](fn_stats_hits_per_gene.png)
*Top 30 ortholog genes by number of overlapping FN lines*

![Deletion size distribution per chromosome (violin)](fn_stats_del_size_by_chrom.png)
*Deletion size distribution per chromosome (violin)*

