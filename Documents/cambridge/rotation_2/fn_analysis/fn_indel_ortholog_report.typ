#set document(title: "FN Indel × Nodulation Ortholog Intersection Report")
#set page(paper: "a4", margin: (x: 2.5cm, y: 2.5cm))
#set text(font: "Libertinus Serif", size: 10pt)
#set heading(numbering: "1.")
#set par(justify: true)
#show table: set text(size: 9pt)
#show figure.caption: set text(size: 9pt, style: "italic")

#align(center)[
  #v(2em)
  #text(size: 20pt, weight: "bold")[FN Indel × Nodulation Ortholog \ Intersection Report]
  #v(0.6em)
  #text(size: 13pt)[_Pisum sativum_ JI2822]
  #v(0.4em)
  #text(size: 10pt)[2026-03-31]
  #v(2em)
]

#outline(depth: 2)
#pagebreak()

= Methods

Fast neutron (FN) deletion and duplication calls (1603 deletions, 451 duplications across 197 FN lines) were intersected with 488 nodulation ortholog gene models mapped to the _Pisum sativum_ JI2822 v1.3 assembly. FN indel coordinates are from skim-sequencing coverage analysis (Supplementary Tables 5--6). Ortholog coordinates were mapped from _Medicago truncatula_ A17 via synteny. An overlap is defined as any FN interval sharing ≥1 bp with a gene model. Hemizygous and homozygous deletion calls are treated independently (a single subline may contribute both, representing the outer and inner coverage thresholds of the same deletion event). Chromosome aliases (chr1--7 ↔ OZ075428--34.1) were resolved via chromAliases.txt.

= Summary Statistics

#figure(
  table(
    columns: (2fr, 1fr),
    table.header([*Metric*], [*Value*]),
    [FN lines], [197],
    [Total deletions], [1603],
    [Total duplications], [451],
    [Mean deletion size], [592 kb],
    [Median deletion size], [230 kb],
    [Mean duplication size], [446 kb],
    [Median duplication size], [250 kb],
    [Total Mb deleted (genome-wide)], [949.0 Mb],
    [Total Mb duplicated (genome-wide)], [201.4 Mb],
    [Ortholog genes hit by ≥1 deletion], [42],
    [Ortholog genes hit by ≥1 duplication], [15],
    [Ortholog genes hit by any indel], [57 / 488],
    [Ortholog genes with no hit], [431],
  ),
  caption: [Summary of FN indel calls and ortholog intersections.],
)

#figure(
  image("plots/fn_stats_ortholog_map.png", width: 100%),
  caption: [Genome-wide distribution of the 488 nodulation ortholog genes coloured by FN indel hit status. Grey: no hit; red: ≥1 deletion; blue: ≥1 duplication; purple: both.],
)
#pagebreak()

= Priority Genes

== LSH1

*Locus:* PSAT_LOCUS25697 · *Position:* chr6:69,313,200--69,314,098 (0.9 kb) · *Description:* _protein LIGHT-DEPENDENT SHORT HYPOCOTYLS 4-like_

2 hit(s) across 1 FN line(s):

#table(
  columns: (1.5fr, 1fr, 1.2fr, 2.5fr, 1fr),
  table.header([*FN subline*], [*Type*], [*Zygosity*], [*Interval*], [*Size (kb)*]),
  [FN2026/7], [deletion], [hemizygous], [chr6:68,280,000--70,230,000], [1950],
  [FN2026/7], [deletion], [homozygous], [chr6:68,300,000--70,230,000], [1930],
)

#figure(
  image("plots/gene_PSAT_LOCUS25697.png", width: 95%),
  caption: [Locus plot for MtLSH1. Each bar is one FN indel drawn to genomic scale; y-axis labels give subline, zygosity, and size. Opacity distinguishes hemizygous (lighter) from homozygous (darker). Gene body shown as black bar (minimum display width applied where gene is too small to see at this scale).],
)

== LSH2

*Locus:* PSAT_LOCUS13554 · *Position:* chr3:370,690,287--370,691,234 (0.9 kb) · *Description:* _protein LIGHT-DEPENDENT SHORT HYPOCOTYLS 4-like_

No FN indels overlap this gene.

== NOOT1

*Locus:* PSAT_LOCUS13168 · *Position:* chr3:317,691,479--317,694,007 (2.5 kb) · *Description:* _NODULE ROOT 1_

2 hit(s) across 1 FN line(s):

#table(
  columns: (1.5fr, 1fr, 1.2fr, 2.5fr, 1fr),
  table.header([*FN subline*], [*Type*], [*Zygosity*], [*Interval*], [*Size (kb)*]),
  [FN3185/1325], [deletion], [hemizygous], [chr3:317,250,000--320,440,000], [3190],
  [FN3185/1325], [deletion], [homozygous], [chr3:317,250,000--320,440,000], [3190],
)

#figure(
  image("plots/gene_PSAT_LOCUS13168.png", width: 95%),
  caption: [Locus plot for MtNOOT1. Each bar is one FN indel drawn to genomic scale; y-axis labels give subline, zygosity, and size. Opacity distinguishes hemizygous (lighter) from homozygous (darker). Gene body shown as black bar (minimum display width applied where gene is too small to see at this scale).],
)

== NOOT2

*Locus:* PSAT_LOCUS26318 · *Position:* chr6:184,006,759--184,009,064 (2.3 kb) · *Description:* _NODULE ROOT 2_

2 hit(s) across 1 FN line(s):

#table(
  columns: (1.5fr, 1fr, 1.2fr, 2.5fr, 1fr),
  table.header([*FN subline*], [*Type*], [*Zygosity*], [*Interval*], [*Size (kb)*]),
  [FN2087/2], [deletion], [hemizygous], [chr6:183,530,000--185,600,000], [2070],
  [FN2087/2], [deletion], [homozygous], [chr6:183,530,000--185,600,000], [2070],
)

#figure(
  image("plots/gene_PSAT_LOCUS26318.png", width: 95%),
  caption: [Locus plot for MtNOOT2. Each bar is one FN indel drawn to genomic scale; y-axis labels give subline, zygosity, and size. Opacity distinguishes hemizygous (lighter) from homozygous (darker). Gene body shown as black bar (minimum display width applied where gene is too small to see at this scale).],
)

== NF-YA1

*Locus:* PSAT_LOCUS26253 · *Position:* chr6:161,422,771--161,424,650 (1.9 kb) · *Description:* _nuclear factor Y subunit A1_

No FN indels overlap this gene.

#pagebreak()
= All Candidate Genes

All 57 ortholog genes overlapped by at least one FN indel, ranked by number of distinct FN lines. Priority genes marked ★.

#table(
  columns: (0.35fr, 1.3fr, 1.3fr, 0.55fr, 2.2fr, 0.65fr, 0.65fr, 0.65fr),
  table.header(
    [*\#*], [*Gene*], [*Locus*], [*Chr*], [*Description*],
    [*Del*], [*Dup*], [*Total*]
  ),
  [1], [PSAT_LOCUS12094], [PSAT_LOCUS12094], [chr3], [—], [4], [0], [4],
  [2], [MtFLS1], [PSAT_LOCUS9418], [chr2], [flavonol synthase], [2], [0], [2],
  [3], [MtGH3], [PSAT_LOCUS9284], [chr2], [Gretchen Hagen3 (auxin-responsive  GH3 family)], [2], [0], [2],
  [4], [MtCKX4], [PSAT_LOCUS13232], [chr3], [cytokinin oxydase 4], [2], [0], [2],
  [5], [MtCASTOR], [PSAT_LOCUS14811], [chr3], [Lotus japonicus CASTOR ion  channel-like], [2], [0], [2],
  [6], [FUN], [PSAT_LOCUS9316], [chr2], [Putative transcription factor bZIP family], [1], [0], [1],
  [7], [IAMT1], [PSAT_LOCUS21425], [chr5], [Putative indole-3-acetate O-methyltransferase], [1], [0], [1],
  [8], [LOG-like|LOG1, cytokinin riboside 5'-monophosphate phosphoribohydrolase], [PSAT_LOCUS13817], [chr3], [—], [1], [0], [1],
  [9], [MEX1Â ], [PSAT_LOCUS30328], [chr7], [—], [0], [1], [1],
  [10], [MtBHLHB], [PSAT_LOCUS31164], [chr7], [bHLH Transcription factor], [1], [0], [1],
  [11], [MtARF4a], [PSAT_LOCUS29732], [chr7], [—], [1], [0], [1],
  [12], [Mt3CHA], [PSAT_LOCUS33237], [chr7], [CCCH-type zinc finger], [0], [1], [1],
  [13], [Mt12b], [PSAT_LOCUS12620], [chr3], [—], [1], [0], [1],
  [14], [MtCKX3], [PSAT_LOCUS5090], [chr1], [cytokinin oxidase 3], [1], [0], [1],
  [15], [MtCKX6], [PSAT_LOCUS20805], [chr5], [cytokinin oxydase 6], [1], [0], [1],
  [16], [MtCAMTA1], [PSAT_LOCUS3780], [chr1], [Calmodulin-binding transcription activator 1], [0], [1], [1],
  [17], [MtBRI], [PSAT_LOCUS23282], [chr5], [Brassinosteroid Insensitive 1 (BR receptor)], [1], [0], [1],
  [18], [MtDGK], [PSAT_LOCUS16447], [chr4], [—], [0], [1], [1],
  [19], [MtCYCA1], [PSAT_LOCUS9381], [chr2], [—], [1], [0], [1],
  [20], [MtEXP7], [PSAT_LOCUS7810], [chr2], [—], [1], [0], [1],
  [21], [MtDMI1], [PSAT_LOCUS3652], [chr1], [does not make infection 1|ion  channel], [0], [1], [1],
  [22], [MtFPN2], [PSAT_LOCUS29222], [chr7], [ferroportin 2], [1], [0], [1],
  [23], [MtGRF1], [PSAT_LOCUS30323], [chr7], [growth-regulating factor 1 (plant-specific transcriptio], [0], [1], [1],
  [24], [MtHOM_TALE_A ], [PSAT_LOCUS17159], [chr4], [Putative transcription factor Homeodomain-TALE-BEL fami], [0], [1], [1],
  [25], [Epidermal Patterning Factor-Like], [PSAT_LOCUS24205], [chr5], [signalling peptide], [1], [0], [1],
  [26], [MtHPT4], [PSAT_LOCUS19295], [chr5], [Histidine PhosphoTransfer protein 4], [1], [0], [1],
  [27], [MtIPT1], [PSAT_LOCUS28589], [chr6], [isopentenyl transferase], [1], [0], [1],
  [28], [MtLSH1 ★], [PSAT_LOCUS25697], [chr6], [protein LIGHT-DEPENDENT SHORT HYPOCOTYLS 4-like], [1], [0], [1],
  [29], [MtLINL], [PSAT_LOCUS33304], [chr7], [—], [0], [1], [1],
  [30], [MtMATE58], [PSAT_LOCUS28149], [chr6], [—], [1], [0], [1],
  [31], [MtMSD1|MtGA20ox|MtGA20ox3], [PSAT_LOCUS28173], [chr6], [MAIN STEM DWARF1 (Gibberellin 20 oxidase|gibberellin bi], [1], [0], [1],
  [32], [MtN21], [PSAT_LOCUS20597], [chr5], [Mtnodulin 21], [1], [0], [1],
  [33], [MtLYR6], [PSAT_LOCUS12623], [chr3], [LYK (LysM domain-RLK)-related 6], [1], [0], [1],
  [34], [MtNCR564], [PSAT_LOCUS7420], [chr1], [nodule-specific cysteine-richÂ (NCR) peptide 564], [1], [0], [1],
  [35], [MtNCR564], [PSAT_LOCUS7423], [chr1], [nodule-specific cysteine-richÂ (NCR) peptide 564], [1], [0], [1],
  [36], [MtNOOT1 ★], [PSAT_LOCUS13168], [chr3], [NODULE ROOT 1], [1], [0], [1],
  [37], [MtNOOT2 ★], [PSAT_LOCUS26318], [chr6], [NODULE ROOT 2], [1], [0], [1],
  [38], [MtNYE1], [PSAT_LOCUS7808], [chr2], [non-yellowing 1 (=stay green 1, SG1)], [1], [0], [1],
  [39], [MtOFP_A], [PSAT_LOCUS33317], [chr7], [Putative transcription factor OFP family], [0], [1], [1],
  [40], [MtPLT3], [PSAT_LOCUS9232], [chr2], [plethora ERF], [1], [0], [1],
  [41], [MtNCR564], [PSAT_LOCUS13843], [chr3], [nodule-specific cysteine-richÂ (NCR) peptide 564], [1], [0], [1],
  [42], [MtSTF|MtWOX1a], [PSAT_LOCUS33504], [chr7], [—], [1], [0], [1],
  [43], [MtSTY2|MtSTYL2], [PSAT_LOCUS17064], [chr4], [SHORT INTERNODES/STYLISH 2|STY-like 2], [0], [1], [1],
  [44], [MtSYMREM1|MtREM2.2], [PSAT_LOCUS32965], [chr7], [symbiotic remorin 1|remorin 2.2], [0], [1], [1],
  [45], [MtSUC2], [PSAT_LOCUS5889], [chr1], [—], [0], [1], [1],
  [46], [MtbZIPB], [PSAT_LOCUS24545], [chr5], [bZIP Transcription factor], [1], [0], [1],
  [47], [PID], [PSAT_LOCUS32565], [chr7], [PINOID], [0], [1], [1],
  [48], [PSAT_LOCUS12014], [PSAT_LOCUS12014], [chr3], [—], [1], [0], [1],
  [49], [PSAT_LOCUS14356], [PSAT_LOCUS14356], [chr3], [—], [1], [0], [1],
  [50], [PSAT_LOCUS32838], [PSAT_LOCUS32838], [chr7], [—], [0], [1], [1],
  [51], [PSAT_LOCUS8846], [PSAT_LOCUS8846], [chr2], [—], [1], [0], [1],
  [52], [PSAT_LOCUS9272], [PSAT_LOCUS9272], [chr2], [—], [1], [0], [1],
  [53], [Putative endo-polygalacturonase], [PSAT_LOCUS9139], [chr2], [—], [0], [1], [1],
  [54], [Putative transcription factor AP2-EREBP family, LEP], [PSAT_LOCUS12889], [chr3], [—], [1], [0], [1],
  [55], [SMR-like ], [PSAT_LOCUS6606], [chr1], [Putative cyclin-dependent protein kinase inhibitor SMR], [1], [0], [1],
  [56], [TRX], [PSAT_LOCUS8847], [chr2], [—], [1], [0], [1],
  [57], [trihelix transcription factor (GT) 4], [PSAT_LOCUS26122], [chr6], [Myb/SANT-like DNA-binding domain protein, PETAL LOSS li], [1], [0], [1],
)

= Locus Plots

One plot per hit gene, drawn to genomic scale. Each bar is one FN indel call; y-axis labels give subline, zygosity, and indel size. The gene body is shown as a black bar beneath (minimum display width applied where gene is too small to see at the indel scale).

#figure(
  image("plots/gene_PSAT_LOCUS12094.png", width: 100%),
  caption: [PSAT_LOCUS12094 · PSAT_LOCUS12094 · chr3 · 4 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS9418.png", width: 100%),
  caption: [MtFLS1 · PSAT_LOCUS9418 · chr2 · 2 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS9284.png", width: 100%),
  caption: [MtGH3 · PSAT_LOCUS9284 · chr2 · 2 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS13232.png", width: 100%),
  caption: [MtCKX4 · PSAT_LOCUS13232 · chr3 · 2 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS14811.png", width: 100%),
  caption: [MtCASTOR · PSAT_LOCUS14811 · chr3 · 2 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS9316.png", width: 100%),
  caption: [FUN · PSAT_LOCUS9316 · chr2 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS21425.png", width: 100%),
  caption: [IAMT1 · PSAT_LOCUS21425 · chr5 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS13817.png", width: 100%),
  caption: [LOG-like|LOG1, cytokinin riboside 5'-monophosphate phosphoribohydrolase · PSAT_LOCUS13817 · chr3 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS30328.png", width: 100%),
  caption: [MEX1Â  · PSAT_LOCUS30328 · chr7 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS31164.png", width: 100%),
  caption: [MtBHLHB · PSAT_LOCUS31164 · chr7 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS29732.png", width: 100%),
  caption: [MtARF4a · PSAT_LOCUS29732 · chr7 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS33237.png", width: 100%),
  caption: [Mt3CHA · PSAT_LOCUS33237 · chr7 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS12620.png", width: 100%),
  caption: [Mt12b · PSAT_LOCUS12620 · chr3 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS5090.png", width: 100%),
  caption: [MtCKX3 · PSAT_LOCUS5090 · chr1 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS20805.png", width: 100%),
  caption: [MtCKX6 · PSAT_LOCUS20805 · chr5 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS3780.png", width: 100%),
  caption: [MtCAMTA1 · PSAT_LOCUS3780 · chr1 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS23282.png", width: 100%),
  caption: [MtBRI · PSAT_LOCUS23282 · chr5 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS16447.png", width: 100%),
  caption: [MtDGK · PSAT_LOCUS16447 · chr4 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS9381.png", width: 100%),
  caption: [MtCYCA1 · PSAT_LOCUS9381 · chr2 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS7810.png", width: 100%),
  caption: [MtEXP7 · PSAT_LOCUS7810 · chr2 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS3652.png", width: 100%),
  caption: [MtDMI1 · PSAT_LOCUS3652 · chr1 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS29222.png", width: 100%),
  caption: [MtFPN2 · PSAT_LOCUS29222 · chr7 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS30323.png", width: 100%),
  caption: [MtGRF1 · PSAT_LOCUS30323 · chr7 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS17159.png", width: 100%),
  caption: [MtHOM_TALE_A  · PSAT_LOCUS17159 · chr4 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS24205.png", width: 100%),
  caption: [Epidermal Patterning Factor-Like · PSAT_LOCUS24205 · chr5 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS19295.png", width: 100%),
  caption: [MtHPT4 · PSAT_LOCUS19295 · chr5 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS28589.png", width: 100%),
  caption: [MtIPT1 · PSAT_LOCUS28589 · chr6 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS25697.png", width: 100%),
  caption: [MtLSH1 ★ · PSAT_LOCUS25697 · chr6 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS33304.png", width: 100%),
  caption: [MtLINL · PSAT_LOCUS33304 · chr7 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS28149.png", width: 100%),
  caption: [MtMATE58 · PSAT_LOCUS28149 · chr6 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS28173.png", width: 100%),
  caption: [MtMSD1|MtGA20ox|MtGA20ox3 · PSAT_LOCUS28173 · chr6 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS20597.png", width: 100%),
  caption: [MtN21 · PSAT_LOCUS20597 · chr5 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS12623.png", width: 100%),
  caption: [MtLYR6 · PSAT_LOCUS12623 · chr3 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS7420.png", width: 100%),
  caption: [MtNCR564 · PSAT_LOCUS7420 · chr1 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS7423.png", width: 100%),
  caption: [MtNCR564 · PSAT_LOCUS7423 · chr1 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS13168.png", width: 100%),
  caption: [MtNOOT1 ★ · PSAT_LOCUS13168 · chr3 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS26318.png", width: 100%),
  caption: [MtNOOT2 ★ · PSAT_LOCUS26318 · chr6 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS7808.png", width: 100%),
  caption: [MtNYE1 · PSAT_LOCUS7808 · chr2 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS33317.png", width: 100%),
  caption: [MtOFP_A · PSAT_LOCUS33317 · chr7 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS9232.png", width: 100%),
  caption: [MtPLT3 · PSAT_LOCUS9232 · chr2 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS13843.png", width: 100%),
  caption: [MtNCR564 · PSAT_LOCUS13843 · chr3 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS33504.png", width: 100%),
  caption: [MtSTF|MtWOX1a · PSAT_LOCUS33504 · chr7 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS17064.png", width: 100%),
  caption: [MtSTY2|MtSTYL2 · PSAT_LOCUS17064 · chr4 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS32965.png", width: 100%),
  caption: [MtSYMREM1|MtREM2.2 · PSAT_LOCUS32965 · chr7 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS5889.png", width: 100%),
  caption: [MtSUC2 · PSAT_LOCUS5889 · chr1 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS24545.png", width: 100%),
  caption: [MtbZIPB · PSAT_LOCUS24545 · chr5 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS32565.png", width: 100%),
  caption: [PID · PSAT_LOCUS32565 · chr7 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS12014.png", width: 100%),
  caption: [PSAT_LOCUS12014 · PSAT_LOCUS12014 · chr3 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS14356.png", width: 100%),
  caption: [PSAT_LOCUS14356 · PSAT_LOCUS14356 · chr3 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS32838.png", width: 100%),
  caption: [PSAT_LOCUS32838 · PSAT_LOCUS32838 · chr7 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS8846.png", width: 100%),
  caption: [PSAT_LOCUS8846 · PSAT_LOCUS8846 · chr2 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS9272.png", width: 100%),
  caption: [PSAT_LOCUS9272 · PSAT_LOCUS9272 · chr2 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS9139.png", width: 100%),
  caption: [Putative endo-polygalacturonase · PSAT_LOCUS9139 · chr2 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS12889.png", width: 100%),
  caption: [Putative transcription factor AP2-EREBP family, LEP · PSAT_LOCUS12889 · chr3 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS6606.png", width: 100%),
  caption: [SMR-like  · PSAT_LOCUS6606 · chr1 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS8847.png", width: 100%),
  caption: [TRX · PSAT_LOCUS8847 · chr2 · 1 FN line(s)],
)

#figure(
  image("plots/gene_PSAT_LOCUS26122.png", width: 100%),
  caption: [trihelix transcription factor (GT) 4 · PSAT_LOCUS26122 · chr6 · 1 FN line(s)],
)

= Supplementary Figures

#figure(
  image("plots/fn_stats_dels_dups_per_line.png", width: 100%),
  caption: [Deletions and duplications per FN line (top 30 by total count).],
)

#figure(
  image("plots/fn_stats_del_size_dist.png", width: 100%),
  caption: [Deletion size distribution (log scale) by zygosity.],
)

#figure(
  image("plots/fn_stats_dup_size_dist.png", width: 100%),
  caption: [Duplication size distribution (log scale).],
)

#figure(
  image("plots/fn_stats_del_size_by_chrom.png", width: 100%),
  caption: [Deletion size distribution per chromosome (violin).],
)

#figure(
  image("plots/fn_stats_genome_coverage.png", width: 100%),
  caption: [Total sequence deleted and duplicated per chromosome (Mb).],
)

