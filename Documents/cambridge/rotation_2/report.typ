// report.typ — Rotation 2 Scientific Report
// Compile with: typst compile report.typ

#import "@preview/wordometer:0.1.4": word-count, total-words
#show: word-count

// ── Page & document metadata ────────────────────────────────────────────────

#set document(
  title: "Investigating the genetic basis of nodulation via chemical mutagenesis screening and natural diversity panel mining in Pisum sativum",
  author: "Michael Beavitt",
  date: datetime(year: 2026, month: 4, day: 2),
)

// Body: Palatino — humanist serif, journal-like; headings: Helvetica Neue sans-serif
#let sans = ("Helvetica Neue", "Helvetica", "Arial")
#let serif = ("Palatino", "New Computer Modern")

#set page(
  paper: "a4",
  margin: (top: 2.5cm, bottom: 2.8cm, left: 2.7cm, right: 2.5cm),
  numbering: "1",
)

// ── Typography ──────────────────────────────────────────────────────────────

#set text(font: serif, size: 10.5pt, lang: "en")
#set par(justify: true, leading: 0.68em, spacing: 1.15em)
#set heading(numbering: none)
#set figure(supplement: [Figure], numbering: "1")
#set figure.caption(separator: [. ], position: bottom)
#set bibliography(style: "apa")

// Level 1: bold sans-serif, sentence case, no rule — like a Nature section heading
#show heading.where(level: 1): it => {
  v(1.8em, weak: true)
  block(below: 0.5em)[
    #text(font: sans, size: 10.5pt, weight: "bold")[#it.body]
  ]
}

// Level 2: bold italic sans-serif
#show heading.where(level: 2): it => {
  v(1.1em, weak: true)
  block(below: 0.65em)[
    #text(font: sans, size: 10.5pt, weight: "bold", style: "italic")[#it.body]
  ]
}

// Captions in sans-serif (Nature style)
#show figure.caption: it => text(font: sans, size: 9pt)[#it]
#show figure: set block(breakable: false)

// ── Utilities ───────────────────────────────────────────────────────────────

#let todo(body) = block(
  fill: rgb("#fffff0"),
  stroke: (left: 2.5pt + rgb("#cccc55")),
  inset: (left: 10pt, right: 8pt, top: 5pt, bottom: 5pt),
  width: 100%,
)[
  #text(size: 9.5pt, fill: rgb("#666600"), style: "italic")[*TODO:* #body]
]

#let ph(body) = text(fill: luma(160), style: "italic")[\<#body\>]

#let panel(lbl, path, width: 100%) = align(left)[
  #text(weight: "bold", size: 9.5pt)[#lbl] \
  #image(path, width: width)
]

// ── Title block ─────────────────────────────────────────────────────────────

#v(2cm)
#text(font: serif, size: 19pt, weight: "bold", tracking: -0.4pt)[
  Investigating the genetic basis of nodulation via
  chemical mutagenesis screening and natural
  diversity panel mining in _Pisum sativum_
]
#v(0.9cm)
#text(font: sans, size: 10pt, weight: "bold")[Michael Beavitt]
#v(0.15cm)
#text(font: sans, size: 9pt, fill: luma(110), style: "italic")[
  Sainsbury Laboratory, University of Cambridge · April 2026
] \
#v(0.1cm)
#text(font: sans, size: 9pt, fill: luma(110), style: "italic")[
  Supervisor: Dr Katharina Schiessl · 4968 words
]

#v(0.8cm)
#text(font: sans, size: 8.5pt, fill: luma(60), style: "italic")[
  I confirm that the material in this report is not copied from any published
  material, nor is it a paraphrase or abstract of any published material, nor
  is it copied or paraphrased from any AI-generated output unless it is
  identified as such and a full source reference is given. I confirm that,
  other than where indicated as above, this document is my own work.
]

#v(0.5cm)
#block(
  stroke: (top: 1pt + luma(180), bottom: 0.5pt + luma(220)),
  inset: (top: 10pt, bottom: 10pt),
  width: 100%,
)[
  #text(font: sans, size: 8pt, weight: "bold", fill: luma(80))[SUMMARY]
  #v(0.4em)
  #set par(justify: true)
  _Pisum sativum_ (common pea) was the original genetics model, prized for its
  well-segregating traits and ease of crossing. It also represents one of the
  most agriculturally relevant nitrogen-fixing plants --- but is its identity as
  a model plant still relevant in the genomics era? Here we
  investigate the potential of pea by exploiting natural and induced diversity
  datasets: an EMS screen of over 1,400 plants identified 11 heritable aberrant
  nodulation phenotypes, a fast neutron panel yielded 57 reverse-genetics
  candidates including predicted knockouts of LSH1, NOOT1, and NOOT2, and
  _k_-mer clustering of 656 _Pisum_ accessions partitioned the diversity panel
  into ten haplotypic groups for natural variant mining. Additionally, we
  demonstrate transgenic overexpression of nodulation-related genes in
  _P. sativum_ cv. Avola via hairy root transformation, providing a functional
  validation system for candidates in an otherwise transformation-recalcitrant
  species.
]
#v(0.5cm)

#pagebreak()

// ── Introduction ─────────────────────────────────────────────────────────────

= Introduction

While all plants require reactive nitrogen to grow, no eukaryotic organisms have
the capability of converting atmospheric N#sub[2] into a reactive form suitable
for use in growth and development. Many prokaryotes are able to fix nitrogen
under anaerobic conditions, and a clade of land plants termed the
nitrogen-fixing clade evolved mechanisms to cultivate nitrogen fixing bacteria
approximately 100 million years ago, with specialised organs (nodules) likely
evolving multiple times within the nitrogen fixing clade since then @doyle2011phylogenetic.
Nodulation mechanisms and pathways are well established in _Medicago
truncatula_, one of the most deeply studied plants in the nitrogen-fixing clade
from a modern reverse genetics perspective. However, shoot phenotypes are
difficult to categorise, and _Medicago_ is far from being put to use in an
agricultural setting. To this end, it is of significant interest to elucidate
nodulation pathways in _Pisum sativum_ (garden pea). Six genes of particular
interest form the core of this study. NIN (NODULE INCEPTION) is a master
transcriptional regulator of nodule organogenesis,
required for both infection thread formation and cortical cell division
@schauser1999plant. NOOT1 and NOOT2 (NODULE ROOT) are BLADE-ON-PETIOLE (BOP)-class transcription
factors that maintain nodule identity by repressing root developmental programmes
within the nodule meristem @magne2018mtnodule. NF-YA1 (Nuclear Factor Y subunit
A1) is a CCAAT-binding transcription factor acting downstream of the symbiosis
signalling cascade to promote nodule organogenesis @laporte2014ccaat. LSH1 and
LSH2 (LIGHT SENSITIVE SHORT HYPOCOTYL) — factors typically implicated as shoot phenotype regulators — have recently
been identified as positive regulators of nodulation @lee2024light, making them
particularly interesting candidates for understanding how nodulation-relevant
programmes may have been co-opted from more general developmental pathways.
Many genes which coordinate nodulation are highly conserved plant genes
involved in general root growth/arbuscular mycorrhiza symbiosis but are likely
co-opted by plants in the nitrogen fixing clade either through modification,
duplication or nodule-specific regulation @yu2024conserved. One working
evolutionary theory is that more general regulatory programmes have been
co-opted by root nodulation systems as opposed to the evolution of a
fundamentally new signalling architecture @griesmann2018phylogenomics.  As a diploid crop species with a long history as a model for plant genetics,
_P. sativum_ offers well-segregating traits and an established nodulation
phenotype, making it a promising system for disentangling the mechanistic
basis of this co-option.

Three resources were made available through collaborators: a natural pea diversity
germplasm dataset along with ~30× short-read paired-end genome sequencing data via
the John Innes Centre (JIC) @feng2025genomic; an ethyl methanesulfonate (EMS)
chemical mutagenesis panel of over 2,000 lines in the "JI2822" pea background,
also via the JIC (unpublished); and a fast neutron (FN) mutagenesis panel, also in the JI2822
background @ellis2026mapped. EMS alkylates guanine residues, causing G/C-to-A/T
transition single nucleotide polymorphisms (SNPs) genome-wide, while fast neutron irradiation induces double-strand
breaks that are repaired imprecisely, resulting predominantly in deletions and
other structural insertion/deletion (indel) variants. Our challenge was to maximise the value of each of these
resources in pursuit of establishing _P. sativum_ as a tractable genetic model for
dissecting the molecular basis of nodulation.

The EMS panel was the product of the PCGIN consortium and our phenotyping and
screening was undertaken as part of a shared effort in two stages: 908 M2 plants were
phenotyped in 2025, with seeds from 44 flagged lines (126 siblings) taken
forward to controlled conditions in the M3 generation. A further 524 M2 plants
were screened in 2026 where 86 additional
candidates were identified that are yet to be advanced to M3. In both stages, plants were flagged for
unusual nodule morphology, hypernodulation, or correlated root and shoot
phenotypes. The FN panel, with its
propensity to generate large deletions, was treated as a reverse genetics
resource. Insertion/duplication calls were intersected with a set of 488
nodulation-relevant gene models manually curated from the _Medicago truncatula_
annotation and mapped to _P. sativum_ JI2822 coordinates via OrthoFinder,
identifying lines likely to carry complete gene knockouts. For the natural
diversity panel, the goal was to interpret the 656 whole genome sequencing accessions (short reads, 30X average
coverage) — identifying genes and accessions that carry the most variation, and
grouping accessions into clusters that may reflect distinct nodule-related haplotypic
backgrounds and thus be informative for candidate selection further on.

A key question is how to validate candidate mutants or natural variant effects
functionally in _P. sativum_.
Pea has a well-deserved reputation for being recalcitrant to transformation
@dalmais2008utilldb @pandey2022optimized,
and our collaborators at the JIC had reported that the JI2822 background was
particularly intractable. We therefore set out to establish hairy root
transformation in the cv. Avola background as a parallel workstream, on the
basis that a tractable transformation system would be essential for any
functional validation of candidates emerging from the screens.

// ── Results ──────────────────────────────────────────────────────────────────

= Results

== M2 greenhouse screen reveals a spectrum of aberrant nodulation phenotypes

To assess the induced genetic diversity of nodulation phenotypes, we screened
an EMS-mutagenised population of second generation (M2) JI2822 pea lines
across two stages in 2025 and 2026. As these lines were grown in the greenhouse
for shoot phenotyping, we relied on nodulation phenotypes induced by natural
ambient rhizobial strains, prioritising unusual nodule morphologies or
hypernodulation, particularly when combined with unusual root or shoot
phenotypes (extremely low-complexity roots, unusually short internodes,
serrated leaves, or large features). Plants without nodules were not flagged,
as the greenhouse conditions made it impossible to distinguish a genuine
non-nodulation phenotype from decayed nodules, lack of nodulation due to sufficient nitrogen, or
failure of the roots to encounter a compatible rhizobial strain.

In total, 1,432 plants were assessed across both stages, of which 439 were
flagged as showing notable phenotypes: 48 were hypernodulators (supernumerary
nodules across the root system), 279 displayed aberrant nodule morphology
(multilobed, clustered, brain-like, or supercluster forms), 40 showed
markedly low root complexity, and the remainder showed other notable features
not easily assigned to a single category. Seeds from 44 lines identified in the 2025 screen
were selected for controlled conditions M3 screening, with 2–4 siblings per line
giving 126 seeds in total.

== M3 controlled conditions screen confirms a subset of M2 phenotypes are heritable

To determine whether the phenotypes observed during the M2 greenhouse screen
reflect true genetic lesions rather than environmental variation, seeds from the
44 selected M2 lines were re-grown as M3 siblings under controlled low-nitrogen
conditions, inoculated with a single isogenic _Rhizobium leguminosarum_ bv.
_viciae_ strain carrying a GUS reporter (blue histochemical stain marking
bacterially colonised nodules). Of the 44 lines assessed, 11 showed consistent
atypical nodulation phenotypes across both the initial controlled conditions
screen and a subsequent microscopy screen, while the remainder appeared
fully wild-type in at least one of the two rounds. This was somewhat surprising,
as the phenotypes recorded during the M2 screening were in many cases quite
dramatic (large multilobed nodule clusters and/or hypernodulating root systems).
The fact that the full EMS population was not grown under conditions conducive to
accurate nodule phenotyping (monoculture of rhizobia, different substrate, etc.)
will have contributed to this difference in phenotype.

Where possible, selected lines were backcrossed to the parent JI2822 line
to begin segregating out background EMS mutations. Due to time and resource
constraints, only five crosses were successfully completed; this is not
critical in the short term, as backcrossing can be performed in the next
generation using seeds that have been self-fertilised in generation M3.

One line in particular stands out: 1205-2, represented by two siblings
(1205-2(1) and 1205-2(2); M3 siblings of M2 line 1205-2), both of which displayed a striking hypernodulation
phenotype — large numbers of small, round, white nodules closely spaced but not
fused, on a large, waxy root system with thick cable-like roots and a prominent
tap root (cf. @fig-1205 D for representative wild-type nodule morphology). The
phenotype was consistent across both siblings, suggesting a genuine genetic
basis. Multiple attempts to cross these plants to JI2822 were unsuccessful, and
the M3 plants themselves bore zero viable seeds, indicating severe sterility.
These lines will likely be re-ordered from the panel to continue investigating
their phenotype (@fig-1205).

#figure(
  image("figures/ems_1205.png", width: 100%),
  caption: [
    Hypernodulation phenotype of EMS line 1205-2.
    *(A)* 1205-2 (1) plant and root system showing dense/waxy roots
    and abundant small nodules. Scale bar = 5 cm.
    *(B, C)* Keyence close-ups of 1205-2 (1) and 1205-2 (2) nodules; round,
    small, and closely spaced, with both uncolonised (white) and colonised
    (blue-stained) forms. Scale bars = 1000 µm.
    *(D)* Wild-type JI2822 for comparison; elongated, individually spaced nodules
    with pink colouration and blue GUS staining.
  ],
) <fig-1205>

Beyond 1205-2, a number of other lines displayed notable nodule morphologies
under controlled conditions, including fused multi-lobed clusters, fan-shaped
nodules, branching and trifurcating forms, and dense arrays of elongated
cylindrical nodules in the hypernodulating line 1310-1 (2)
(@fig-ems-phenotypes). None of the 11 candidate lines showed obvious shoot phenotypes at this stage.

#figure(
  image("figures/ems_phenotype_collage.png", width: 100%),
  caption: [
    Representative atypical nodule morphologies observed in M3 EMS lines
    under controlled single-strain inoculation conditions. All roots were
    GUS-stained following inoculation with _Rhizobium leguminosarum_ bv.
    _viciae_ Rlv3841-GUS; blue staining indicates active bacterial colonisation.
    Scale bars = 1000 µm (embedded in source images).
    *(A)* Fused multi-lobed nodule cluster, line 1077-1 (2).
    *(B)* Four-lobed nodule with apical splitting, line 1132-2 (1).
    *(C)* Three-lobed nodule with apical splitting; note irregular lobation,
    line 1132-2 (1) (different individual, illustrating the range of forms
    within this line).
    *(D)* Large flat fan-shaped nodules, line 1237-1 (1).
    *(E)* Hypernodulation with densely packed elongated cylindrical nodules,
    line 1310-1 (2).
    *(F)* Multilobed club-shaped fused nodule, line 1310-1 (3).
  ],
) <fig-ems-phenotypes>

Accession 614-2 is a representative example. During the M2 JIC greenhouse screen
it displayed two atypical nodule phenotypes: large, brain-like lobed nodule clusters
and a population of elongated, distorted nodules (@fig-614 A–B). Two M3 sibling
plants — 614-2 (1) and 614-2 (2) — grown from seeds of this M2 line under
controlled laboratory conditions both showed root systems and nodule morphology
indistinguishable from the wild-type JI2822 parent line (@fig-614 C–H): compact,
individually spaced, pale-pink nodules along lateral roots, with no cluster
formation or extraordinary elongation.

// ── 614-2 figure ─────────────────────────────────────────────────────────────

#figure(
  image("figures/614-2_controlled_conditions.png", width: 100%),
  caption: [
    Comparison of EMS line 614-2 nodule phenotypes between the M2 JIC
    greenhouse screen and M3 controlled conditions, alongside the JI2822
    wild-type reference. M2 = second mutant generation, screened at JIC; M3 =
    seeds collected from M2 plants and grown in the laboratory under controlled
    single-strain inoculation conditions.
    *(A)* Brain-like lobed nodule clusters observed on the M2 614-2 plant during
    the JIC screen. *(B)* Elongated, distorted nodules also present on the same
    M2 plant. No scale bars for A–B (no calibration available for M2 screen
    images).
    *(C–D)* Complete root systems of M3 siblings 614-2 (1) and 614-2 (2), 28
    days post-inoculation.
    *(E)* JI2822 wild-type root system for comparison.
    *(F–G)* Close-up nodule morphology of 614-2 (1) and 614-2 (2); individually
    spaced, compact nodules with no cluster formation or extraordinary elongation.
    *(H)* JI2822 wild-type nodule close-up for comparison.
    Scale bars: C–E = 5 cm; F–H = 5 mm.
  ],
) <fig-614>

One fairly easy to test hypothesis was also considered: that the increase in root
density as the plants grew to fill their small pots may have both reduced nutrient
availability and compressed the surrounding earth to trigger nodulation and induce
non-wild-type phenotypes. To test this, we grew the pea cultivar Avola in four conditions: non-compressed
and compressed 1:1 F2 compost/vermiculite mix, our standard nodulation assay
1:1 Terragreen/vermiculite mix, and 1:1 sharp sand/Terragreen mix, to understand
whether substrate had an effect on nodule morphology. Large clusters of nodules,
broadly resembling the hypernodulating phenotypes observed during the EMS screen,
were reproducible in this experiment, and were most pronounced in the compressed
F2/vermiculite and 1:1 sharp sand/Terragreen conditions.

== FN mutagenesis identifies candidate knockout lines in LSH1, NOOT1, and NOOT2

Deletion and duplication calls from 197 FN lines were intersected with 488
nodulation ortholog gene models mapped to the JI2822 assembly. Of these, 57
genes (11.7%) were overlapped by at least one FN indel — 42 by deletions and
15 by duplications, with some genes hit by both (@fig-fn-overview). The large
majority of hits arose from single FN lines, reflecting the expected
stochastic distribution of FN lesions across the panel.

Among the six core nodulation genes of primary interest, three carry FN deletions:
LSH1 (PSAT_LOCUS25697; FN line FN2026/7, ~1.93 Mb deletion on chr6),
NOOT1 (PSAT_LOCUS13168; FN line FN3185/1325, ~3.19 Mb deletion on chr3),
and NOOT2 (PSAT_LOCUS26318; FN line FN2087/2, ~2.07 Mb deletion on chr6;
@fig-fn-loci).  LSH2, NIN, and NF-YA1 were not overlapped by any FN indel in
the panel.

#figure(
  image("figures/fn_dataset_overview.png", width: 100%),
  caption: [
    Overview of the fast neutron indel dataset.
    *(A)* Size distributions of deletions and duplications (log scale).
    *(B)* Distribution of total indel calls per FN line.
    *(C)* Total FN indel coverage (sequence deleted and duplicated) per chromosome (Mb).
    *(D)* Genome-wide distribution of the 488 nodulation ortholog genes coloured
    by FN hit status (red: deletion; blue: duplication; purple: both; grey: no hit).
  ],
) <fig-fn-overview>

#figure(
  image("figures/fn_priority_composite.png", width: 100%),
  caption: [
    FN deletion locus plots for three core nodulation genes with hits.
    Each bar represents the homozygous deletion call drawn to genomic scale;
    bar label gives deletion size. Gene body indicated by solid black bar.
    *LSH1* (PSAT_LOCUS25697, chr6): FN line FN2026/7, ~1.93 Mb deletion.
    *NOOT1* (PSAT_LOCUS13168, chr3): FN line FN3185/1325, ~3.19 Mb deletion.
    *NOOT2* (PSAT_LOCUS26318, chr6): FN line FN2087/2, ~2.07 Mb deletion.
  ],
) <fig-fn-loci>


== Hairy root transformation in _P. sativum_ cv. Avola enables LSH overexpression and yields an altered root morphology phenotype <sec-hairy-root>

To assess the feasibility of hairy root transformation as a tool for functional
characterisation of LSH1 and LSH2 in _P. sativum_, transformation of
_P. sativum_ cv. Avola was attempted following the protocol of Leppyanen et al.
@leppyanen2019agrobacterium. The JI2822 background used for the EMS and FN panels
was reported by collaborators at the JIC to be recalcitrant to transformation,
including by hairy root approaches; cv. Avola was therefore selected as a more
tractable starting point. Three _A. rhizogenes_ strains (AR1193, ARqua1, and K599)
were evaluated with a reporter-only construct carrying DsRed under the _Lotus
japonicus_ ubiquitin promoter (pLjUBI::DsRed), with _n_ = 10 plants per strain.
Concurrently, ~18 plants each were transformed with constructs additionally
carrying LSH1 or LSH2 (pLjUBI::DsRed-LSH1 and pLjUBI::DsRed-LSH2), intended to
generate overexpressor composite plants (see Methods for full protocol).

A total of six DsRed-expressing hairy roots were recovered over a 30-day period
from inoculation to screening: four from K599, one from AR1193, and one from
ARqua1 (@fig-dsred). All successful transformants carried
the reporter-only construct (pLjUBI::DsRed), suggesting that hairy root
transformation is feasible in cv. Avola, though transformation efficiency was low.

// ── Setup figure ─────────────────────────────────────────────────────────────

#figure(
  image("figures/setup_composed.png", width: 100%),
  caption: [
    Hairy root transformation setup for _P. sativum_ cv. Avola.
    *(A)* Etiolated seedlings germinated in BNM in a ventilated tip box prior
    to hypocotyl sectioning and inoculation.
    *(B–C)* Seedling support using windmill-folded filter paper in glass jars
    containing modFP agar: seedling hypocotyl threaded through the central
    aperture of the filter paper (B); seedling in final position with root tip
    in contact with modFP medium (C).
  ],
) <fig-setup>

// ── DsRed reporter figure ─────────────────────────────────────────────────────

#figure(
  image("figures/dsred_avola_hairy_roots.jpeg", width: 95%),
  caption: [
    DsRed fluorescence in reporter-only (pLjUBI::DsRed) transformed hairy roots
    of _P. sativum_ cv. Avola. *(A–D)* Four transformants recovered from strain
    K599. *(E)* Transformant from AR1193. *(F)* Transformant from ARqua1.
    Scale bars = 5 mm.
  ],
) <fig-dsred>

In contrast, the LSH overexpressor constructs, introduced via strain AR1193,
achieved transformation efficiencies of 16/18 and 14/18 for the LSH1 and LSH2
overexpressor constructs, respectively, over a 23-day period from inoculation to
final screening. While the transformed roots were still mostly in early stages of
growth (@fig-lsh1; @fig-lsh2), the LSH2 overexpressor line showed a pronounced helical root curvature and a
marked reduction in lateral root number (@fig-curly).
Seven lines with the LSH2 overexpressor construct and five lines with the LSH1
overexpressor construct were taken forward for attempted stable line cultivation by
surgically removing transformed calli to transformation medium to promote
hormone-mediated re-differentiation; meanwhile, seven lines with the LSH2
overexpressor and eleven with the LSH1 overexpressor were kept to grow in 500 ml
pots in 1:1 Terragreen/vermiculite mix supplemented with ~100 ml modFP medium,
along with the six successfully transformed pLjUBI::DsRed control plants.

// ── Curled root figure ────────────────────────────────────────────────────────

#figure(
  image("figures/curly_lsh2oe.jpg", width: 65%),
  caption: [
    Altered root morphology in a pLjUBI::DsRed-LSH2 composite plant (_P. sativum_
    cv. Avola, line 43224), showing pronounced helical curvature and reduced
    lateral root formation. Scale bar = 5 cm.
  ],
) <fig-curly>

// ── LSH2 composite figure ─────────────────────────────────────────────────────

#figure(
  image("figures/hairy_root_lsh2oe_v2.png", width: 100%),
  caption: [
    Brightfield and fluorescence imaging of pLjUBI::DsRed-LSH2 transformed hairy
    roots in _P. sativum_ cv. Avola (line 43224). *(A–C)* Brightfield
    stereomicroscopy images of three independent transformed hairy roots recovered
    following _A. rhizogenes_ (AR1193) inoculation. *(D–F)* Corresponding DsRed
    fluorescence images confirming stable transgene expression in transformed root
    tissue. Helical root curvature and reduced lateral root formation are visible
    in panels C and F. Scale bars = 5 mm.
  ],
) <fig-lsh2>

// ── LSH1 composite figure ─────────────────────────────────────────────────────

#figure(
  image("figures/hairy_root_lsh1oe_v2.png", width: 100%),
  caption: [
    Brightfield and fluorescence imaging of pLjUBI::DsRed-LSH1 transformed hairy
    roots in _P. sativum_ cv. Avola (line 43225). *(A–C)* Brightfield
    stereomicroscopy images of three representative transformed hairy roots
    recovered following _A. rhizogenes_ (AR1193) inoculation. *(D–F)*
    Corresponding DsRed fluorescence images confirming stable transgene expression.
    Scale bars = 5 mm.
  ],
) <fig-lsh1>

// ── kmer diversity panel section ─────────────────────────────────────────────

== _k_-mer clustering organises natural diversity across pea nodulation gene loci <sec-kmer>

To explore natural sequence variation at nodulation gene loci across the
diversity panel, hierarchical clustering was performed on the binary _k_-mer
presence/absence matrices as described in Methods, using pairwise Euclidean
distances across all 488 ortholog gene loci with Ward linkage; the dendrogram
was cut to yield 10 discrete clusters (@fig-kmer-global).

The per-gene clustermaps for NOOT2 and NIN illustrate the resolution available
at individual loci (@fig-kmer-noot2; @fig-kmer-nin). Both genes show clear
differences in mutation profiles across accessions, with distinct groups lacking
_k_-mer coverage at specific exonic and intronic positions, consistent with
structural variation or SNP-driven _k_-mer dropout at those sites. Notably,
NIN displays considerably greater variance than NOOT2, including in flanking
regions outside the gene body itself, whereas NOOT2 _k_-mer dropout is more
localised. The same lower-variance profile is shared by NOOT1 (not shown),
suggesting this is a property of the NOOT loci rather than an artefact of any
individual gene window.

// ── Global kmer clustermap ────────────────────────────────────────────────────

#figure(
  image("figures/kmer_clustering/global_clustermap_full.png", width: 100%),
  caption: [
    Global hierarchical clustering of 656 accessions by _k_-mer presence/absence
    across 488 nodulation gene orthologs. Accessions are ordered by Ward linkage
    dendrogram (top); annotation bars indicate cluster assignment (10 clusters)
    and species identity. Heatmap shows mean _k_-mer presence across all 488
    ortholog genes, sorted by ascending mean presence; core nodulation genes
    (LSH1, LSH2, NIN, NOOT1, NOOT2, NF-YA1) are labelled.
  ],
) <fig-kmer-global>

// ── NOOT2 per-gene clustermap ─────────────────────────────────────────────────

#figure(
  image("figures/kmer_clustering/MtNOOT2.png", width: 100%),
  caption: [
    Per-accession _k_-mer presence/absence across the NOOT2 gene window. Each
    row is one accession (ordered by Ward linkage); each column is a variable
    _k_-mer position. The annotation bar indicates genomic feature context
    (promoter, exon, intron, downstream). Blue = _k_-mer present; white = absent.
  ],
) <fig-kmer-noot2>

// ── NIN per-gene clustermap ───────────────────────────────────────────────────

#figure(
  image("figures/kmer_clustering/MtNIN.png", width: 100%),
  caption: [
    Per-accession _k_-mer presence/absence across the NIN gene window.
    Layout as in @fig-kmer-noot2.
  ],
) <fig-kmer-nin>

// ── Methods ──────────────────────────────────────────────────────────────────

= Methods

== EMS Panel Screening

During the screening of the EMS mutagenesis panel, qualitative shoot (internode
length, height, leaf size, stipule size, petiole length, peduncle length) and root
characteristics (nodule abundance, nodule morphology, root density, root
complexity, nodule size) were recorded for each plant. In total, seeds from 44 unique segregating M2 lines
were requested, with 2–4 siblings per line, giving 126 seeds.

== Controlled Conditions Experiments

Seeds from 44 flagged M2 lines (126 siblings) were
germinated and cultivated under controlled conditions across three separate growth
batches. Prior to germination, seeds were surface-sterilised (1 min 70% EtOH, 5 min
2% NaOCl, three distilled water rinses, 20-min imbibition) and planted in
clear 1-litre plastic pots without drainage holes containing 1:1
vermiculite/Terragreen substrate watered with 400ml of
nitrogen-free rooting solution, inoculated with 1 ml of _Rhizobium leguminosarum_ bv. _viciae_ strain Rlv3841
carrying a chromosomally integrated, constitutively expressed β-glucuronidase
(GUS) reporter (construct MK4J18, streptomycin-selected), and wrapped in cling
film and left to grow for 28 days in an 18h light cycle growth chamber at 21ºC.
The plants were then removed from the
pots, roots washed, and plants phenotyped (weight, shoot and root
characteristics as above, and height). A representative sample of the root
system was taken for histochemical GUS staining (as detailed in Lee et al. 2024
@lee2024light); stained roots were subsequently imaged under a Keyence VHX 7000
digital microscope for close-up morphological assessment. Plants were then
repotted in 1:1 vermiculite/F2 compost (Levington F2 general purpose compost)
to set seed.

For the substrate comparison experiment, the pea cultivar Avola was grown in
500ml pots in non-compressed 1:1 F2 compost/vermiculite mix, compressed 1:1 F2
compost/vermiculite mix, standard 1:1 Terragreen/vermiculite nodulation assay
mix, and 1:1 sharp sand/Terragreen mix. Seeds were sterilised and inoculated as
described above. Compression was achieved by filling the pots to the brim and
using a stamp tool to compact the compost before planting - this was only
performed for the compost-based substrates.

== Hairy Root Transformation

Hairy root transformation of _P. sativum_ cv. Avola was performed following
@leppyanen2019agrobacterium with minimal modifications. In brief, etiolated
seedlings were germinated in buffered nodulation medium (BNM;
@wickens2026spatiotemporal) in a ventilated tip box (@fig-setup A), sectioned
in the hypocotyl region, and wounded root tips were inoculated by contact with
_A. rhizogenes_ growing on plate. Seedlings were then supported through
windmill-folded filter paper in glass jars containing modified Fahraeus (modFP)
agar, without antibiotic selection, and incubated for 10–14 days at 21°C under
a 16 h/8 h light/dark cycle (@fig-setup B–C). Roots were assessed for DsRed
fluorescence by epifluorescence stereomicroscopy.

== Identification of _P. sativum_ Nodulation Gene Orthologs

A curated gene list of nodulation-relevant genes was compiled from the _Medicago
truncatula_ A17 ecotype (v5r1.7) annotation, comprising genes with established
or candidate roles in nodulation signalling, organogenesis, and regulation,
including LSH1, LSH2, NOOT1, NOOT2, NF-YA1, and related family members
@magne2018mtnodule @laporte2014ccaat. For each A17 gene, a corresponding
ortholog in the R108 ecotype had been manually curated from the _M. truncatula_
R108 Hi-C (high-throughput chromosome conformation capture) assembly
(medtr.R108.gnmHiC_1.ann1; referred to as R108 throughout) annotation, stored as
the `R108_orthologs_final` field in the gene list.

To identify _P. sativum_ JI2822 orthologs, OrthoFinder @emms2019orthofinder was
run on predicted protein sequences from all three species: _M. truncatula_ A17,
_M. truncatula_ R108, and _P. sativum_ JI2822. _P. sativum_ orthologs were
identified using a two-step bridging strategy via R108: for each gene in the gene
list, the pre-curated R108 ortholog ID was used to look up its OrthoFinder
orthogroup, and all _P. sativum_ JI2822 proteins within that orthogroup were
extracted as candidate orthologs. This was supplemented with pairwise
R108-to-_P. sativum_ ortholog assignments from the OrthoFinder pairwise Orthologues
output. Orthogroups containing
exactly one gene from each species (1:1:1 single-copy orthogroups) were flagged as
highest-confidence assignments. The resulting _P. sativum_ locus IDs (PSAT_LOCUS
identifiers) were mapped to genomic coordinates in the JI2822 v1.3 assembly
@rayner2025rebalancing.

== Fast Neutron Indel Intersection

FN deletion and duplication calls (_n_ = 1,603 deletions and 451 duplications
across 197 FN lines) were obtained from skim-sequencing coverage analysis of the
JI2822 FN panel @ellis2026mapped. Deletions were
called at two stringency thresholds per line: a hemizygous call (coverage ~50% of
genome average, the less confident call) and a homozygous call
(coverage ~0%, the inner high-confidence call). These were intersected with the
488 ortholog gene models in JI2822 coordinates using a minimum 1 bp overlap
criterion.

== Natural Diversity Panel: _k_-mer Analysis

To identify sequence variation at nodulation gene loci across the natural
diversity panel, a _k_-mer-based presence/absence approach was applied to the
656-accession dataset with the help of collaborators with access to the data
(B. Steuernagel, JIC, unpublished method). Full variant calling across all
accessions would in principle yield higher-resolution SNP and indel calls, but
is prohibitively computationally intensive at this panel size; the _k_-mer
approach offers a practical alternative that captures presence/absence
variation across gene windows without per-accession alignment and genotyping.
Individual _k_-mer libraries were
generated for each accession using KMC @deorowicz2015kmc (_k_ = 31).  A
corresponding _k_-mer library was generated from the input gene sequences (the
_P. sativum_ ortholog set described above). KMC was then used to compute the
intersection of each accession library with the gene library, yielding, for
each gene and each accession, the subset of _k_-mers present in both. These
intersections were merged into a binary presence/absence matrix ([0, 1], where
0 is presence of a kmer and 1 is absence) indexed by position (columns) and
accession (rows).

== Hierarchical Clustering of Accessions by _k_-mer Diversity

Hierarchical clustering was performed on binary _k_-mer presence/absence
matrices as follows.

For the per-gene analysis, the presence/absence matrix for each gene (positions
× accessions) was clustered over accessions using Ward linkage on Euclidean
distances, implemented via the `scipy.cluster.hierarchy` module (SciPy v1.16.3). Only variable
positions (those not fixed at 0 or 1 across all accessions) were retained.
Results were visualised as clustermaps with a genomic feature annotation bar
(promoter, exon, intron, downstream) derived from the GFF3 annotation.

For the global analysis, pairwise Euclidean distances between all accessions
were computed across the full set of variable _k_-mer positions in all 488
ortholog genes, and Ward linkage was applied to the resulting distance matrix
using `scipy.cluster.hierarchy`. The dendrogram was cut into
discrete clusters using `fcluster`. The global heatmap displays mean _k_-mer presence per gene across all 488
orthologs, with rows sorted by ascending mean presence to place the most
variable genes at the top; the six core nodulation genes are labelled.


// ── Discussion ────────────────────────────────────────────────────────────────

= Discussion

Here we have exploited both natural and induced sequence variation in _Pisum
sativum_ to identify candidate lines for the functional characterisation of key
nodulation regulators. Across the three resources surveyed, 57 FN indel lines
were identified as reverse genetics candidates (including knockouts in NOOT1,
NOOT2, and LSH1), 11 EMS lines showed heritable aberrant nodulation phenotypes
(including 1205-2, 1132-2, 1310-1, and 8 further lines), with 86 additional
candidates from 524 plants screened in 2026 awaiting M3 follow-up, and a
structured diversity panel of 656 accessions clustered by _k_-mer variation
across 488 nodulation-related gene loci, ready to be mined for natural
variation.

In addition, we demonstrated a successful transient transformation experiment
in cv. Avola, allowin us to use _Agrobacterium rhizogenes_-mediated
transformation to assess gene function directly in composite plants after
candidate genes are identified. One notable discrepancy warrants comment: the
reporter-only (pLjUBI::DsRed) experiments yielded only six transformants across
~30 plants, whereas the LSH overexpressor constructs
introduced via AR1193 achieved efficiencies of 16/18 and 14/18 for LSH1 and
LSH2 respectively. The basis for this difference is unclear — slight variation
in bacterial culture preparation or inoculum density between experiments may
have contributed, and more controlled side-by-side comparisons will be needed
to resolve it.  Nonetheless, the result demonstrates that AR1193 is a capable
transformation strain in this system, and that the protocol provides a workable
route to functional validation of candidate genes emerging from the EMS and FN
screens.

The LSH2 overexpressor lines also produced a notable root phenotype: pronounced
helical curvature and a reduction in lateral root number. A comparable reduction
in lateral root formation has been reported in LSH overexpressor lines in
_M. truncatula_ @lee2024light, suggesting that this aspect of LSH2 function
may be conserved between the two species. Whether the helical curvature is a
_P. sativum_-specific response or reflects overexpression strength remains to
be determined as lines mature.

As detailed in the results (Controlled Conditions), the dramatic nodule cluster
phenotypes observed during the M2 greenhouse screen largely did not persist
when lines were re-grown as M3 siblings under controlled conditions. The
substrate comparison experiment in cv. Avola provides a plausible explanation:
physical compression of the substrate appears sufficient to drive cluster
formation, raising the possibility that many of the phenotypes flagged during
the EMS screen were environmentally rather than genetically determined. The
experiment was conducted in cv. Avola rather than JI2822 and with only four
plants per condition, so the result should be interpreted cautiously;
nonetheless it is an important caveat when prioritising EMS lines for further
investigation. A further confound is plant age: several of the lines that
retained atypical morphologies under controlled conditions (including those
shown in @fig-ems-phenotypes) were still relatively young at the time of
phenotyping, and it is possible that more pronounced phenotypes had not yet
fully manifested. These lines will be re-examined at end of life to determine
whether aberrant morphologies become more apparent with nodule maturation, and
whether the phenotypes are reproducible across siblings.

The _k_-mer diversity analysis reveals a striking difference between NIN and the
NOOT loci: NIN shows considerably greater inter-accession variance, including
in flanking regions outside the gene body, whereas NOOT1 and NOOT2 show low
variance across the panel. The basis for this is not immediately clear. NIN is
an ancient regulator of nodulation @liu2020evolution, but both NIN and the NOOT
genes would have been present in the common ancestor of extant pea accessions
and should have had equivalent time to accumulate intra-clade variation.
Whether the contrast reflects differences in selection pressure, local
recombination rate, or gene-family dynamics remains an open question.

A full characterisation of the population structure visible in @fig-kmer-global
is outside the scope of this report, but a few features are worth noting.
Cluster 7 corresponds entirely to _P. fulvum_, and cluster 8 contains no _P.
sativum_, being composed exclusively of _P. abyssinicum_ and _P. sativum_
subsp. _elatius_ accessions — which are themselves distinguishable within the
cluster by eye in the heatmap. _P. sativum_ subsp. _elatius_ is also present
across several other clusters, notably cluster 10, while a single _P.
abyssinicum_ accession appears in cluster 1. The remaining clusters are not
obviously explained by species identity and likely reflect genuine
within-_P. sativum_ haplotypic diversity. This partitioning provides a
principled basis for selecting representative accessions to request for wet lab
follow-up, and opens the door to a more rigorous investigation of natural
variation in nodulation gene content across the genus — including potential
geographic associations that would require the full variant-calling resolution
this _k_-mer approach intentionally trades away.

While the _k_-mer approach facilitated the analysis of a huge number of
accessions' sequencing data, we necessarily lose out on variant resolution and
the specific nature of each variant.  An alternative approach to diversity
mining at these loci would be to lift over the existing variant call tables
produced for the ZW6 cultivar into the JI2822 coordinate space via a chain
file, which would give direct SNP and indel calls rather than _k_-mer dropout
as a proxy for variation. With more time this would have been worth pursuing;
the main obstacle is generating the chain file itself.  A genome-to-genome
alignment between ZW6 and JI2822 was attempted for this purpose, but the
alignment failed to produce a usable result — the pea genome is highly
repetitive and the minimap2 @li2018minimap2 aligner was not able to allocate enough memory for the
indexes, even given ~500GB of RAM. A more sophisticated masking strategy prior
to alignment, for example using a pea-specific repeat library, would likely be
necessary to make this work.

In summary, the forward and reverse genetics approaches described form a strong
basis for investigating nodulation pathways in _P. sativum_ more deeply. The
most obvious barrier preventing researchers from adopting it as a model species
is its resistance to classical transformation methods, but the hairy root
transformation system demonstrated in cv. Avola, while requiring further
optimisation, shows that _Agrobacterium rhizogenes_-mediated transformation is
achievable in pea - albeit for root-related traits alone.  It should be noted
that one original motivation for working in this species was the prospect of
identifying correlated shoot phenotypes that betray the co-option of general
developmental regulators; this has not yet been realised in the current
candidate set, but remains a compelling reason to pursue thorough phenotyping
as lines are stabilised and advanced. How we will validate shoot-root
co-regulation in nodulation pathways is not clear, but stable transformation
methods may be investigated - what is clear is that fast neutron knockouts
covered a significant proportion of the nodulation-related genes we collated,
which will help in categorising knockout phenotypes in both shoot and root.

// ── References ────────────────────────────────────────────────────────────────

#bibliography("refs.bib", title: "References")
