#set document(
  title: "Understanding the Evolutionary Landscape of Centromere Evolution in Arabidopsis Using Monte Carlo Simulations",
  author: "Michael Beavitt",
  date: datetime(year: 2026, month: 1, day: 5),
)

#set page(
  paper: "a4",
  margin: (x: 2.5cm, y: 2.5cm),
  numbering: "1",
)

#set text(
  font: "New Computer Modern",
  size: 11pt,
)

#set par(justify: true)

#set heading(numbering: "1.1")

// Title page
#align(center)[
  #text(size: 16pt, weight: "bold")[
    Understanding the Evolutionary Landscape of Centromere Evolution in Arabidopsis Using Monte Carlo Simulations
  ]

  #v(1em)

  #text(size: 12pt)[
    Name: Michael Beavitt \
    Email: #link("mailto:mab282@cam.ac.uk")[mab282\@cam.ac.uk] \
    Department: Department of Plant Sciences \
    Supervisor: Professor Ian Henderson \

    #v(0.5em)
    Date: 05/01/2026
  ]
]

#v(2em)

*Word count:*

*Statement:*

#pagebreak()

= Introduction

Eukaryotic centromeres are essential components which allow for accurate segregation of chromosomes and chromatids during cell division. The relationship between the DNA composition of centromere sequences and their function has been long questioned -- so much so, that it has been described as the "centromere paradox" @henikoff2001centromere. While it's clear that the inheritance of centromeric function is primarily epigenetic, there still remains the enigma of the abundance of satellite repeat arrays across the eukaryotes. Do the satellite repeats themselves confer some kind of fitness to the centromere, making them more likely to be inherited during meiosis? If satellite repeats appear at centromeric loci, are the molecular components of cell division themselves -- the kinetochore -- responsible for generating satellite repeats somehow? With recent advancements in long read sequencing technologies, gap-free centromere sequences in Arabidopsis have been resolved for the first time, giving us insights into their composition @naish2021genetic @wlodzimierz2023cycles. By using these resources to turn our attention to satellite array composition and understand how they might have come to be, perhaps we can begin to understand why they exist.

While the mechanisms and drive of centromere evolution are uncertain, it is undisputed that turnover is much faster in centromere arrays than the rest of the genome, not only in the case of unusual frame-preserving INDELs but also single nucleotide polymorphisms (SNPs) which are nine times more frequent than the arms of the chromosomes, occurring at a rate of 6.12 × 10#super[−8] per bp per generation. Large inversions and tandem duplications have also been observed @wlodzimierz2023cycles. These dynamics are not just limited to Arabidopsis, but occur across Eukaryotic lineages @henderson2025cyclical. Transposable element (TE) insertion is also highly frequent in centromeric arrays, with concrete evidence for direct targeting to centromere-specific proteins @tsukahara2025centrophilic.

Many purport that the appearance of higher order repeat (HOR) sequences and centromere identity is a direct result of feedback mechanisms involving the kinetochore machinery, but the recent simulation work of Dong et al challenges this, showing that high frequency randomly inserted/deleted sequences combined with low frequency random single nucleotide polymorphisms is sufficient to generate centromere-like repeat arrays. This finding is surprising, given that we know that at least one driver of evolution is directly related to kinetochore machinery given Tsukahara et al's finding that transposon insertions are targeted to CENH3 chromatin specifically, which exclusively occupies subsets of centromeric DNA.

What I wanted to do in this rotation project was to understand how the model proposed by Dong et al truly captures the structures and sequence diversity of repeat arrays in real genomes, and additionally modify the existing model to capture a hypothesised mechanism in eukaryote genomes -- kinetochore-associated recombination. While the original model produced intriguing results, only five centromeres were simulated. Large scale comparisons of hundreds of simulated centromere satellite arrays to their real counterparts would give new insights into the strengths and limitations of the model, and potentially help us to understand the true mechanisms at play. The recent publication of 66 _Arabidopsis thaliana_ genomes, containing 330 well-resolved centromere sequences @wlodzimierz2023cycles provided a clear route forward to extend the analysis.

As touched upon already, this project also aimed to build upon the existing simulation framework by applying mechanisms inspired by the "KARMA" model (Kinetochore-associated Recombination Machinery in Arabidopsis) in an attempt to both understand its feasibility and explore the consequences of accepting the assumptions it represents.

In brief, the KARMA model hypothesises that kinetochore protein complexes attract recombination factors which increase the turnover of satellite DNA in close proximity to the kinetochore itself. This, in theory, homogenises the array, purging non-similar repeats or invading transposon sequences and producing the highly homogeneous satellite arrays seen in nature. The motivation behind applying it to the simulation framework was to understand the implications of accepting the assumptions of this popular model -- would they produce unnatural or incorrect centromeric repeats, and if so in what manner? On the other hand, would adding this positive feedback loop, reinforcing self-similar regions of satellite DNA produce more "centromeric" sequences than the existing model?

= Methods

== Code Optimisation

The original Dong et al paper's simulation code, published on GitHub (#link("https://github.com/schneebergerlab/replicated-assemblies-centromere-study")) provided a comprehensive starting point for simulating the three desired events -- point mutations, insertions and deletions. The original program unfortunately took ~3 days to complete a single 6 million generation run, which meant running the desired ~600 simulations would use approximately 43200 CPU hours or almost 5 years in real time on a single core. Multithreading would go some way to making this time constraint more manageable, but a better approach would involve optimising the core program as much as possible before applying multithreading.

Optimisations initially involved removing unneeded checks and calculations. In particular, the original program performed complex bounds checking to ensure each of the 178 bp repeats was indexed correctly, tracking repeat positions in a separate table. This was only necessary if repeat length could change, but all mutations modelled (SNPs, insertions and deletions in multiples of 178 bp) were frame-preserving, so this situation never arose. I therefore replaced the bounds checking with simple modulo arithmetic to index each repeat. Combined with replacing the most expensive array calculations with vectorized NumPy @harris2020array operations, this reduced runtime to ~20 minutes.

Tracking each base individually remained a major bottleneck, with the array size scaling as 178 × repeat number (e.g. 2.67 million characters for 15,000 repeats). This made O(n) insertion and deletion operations extremely costly, so an abstracted repeat-resolution representation was implemented to enable more efficient list operations. This reduced runtime further to ~3 minutes, making hundreds of simulations feasible.

To implement the kinetochore-associated recombination model, it was necessary to calculate array self-similarity at each generation to identify self-similar regions for targeted INDEL placement. Array self-similarity was initially calculated using the fractal dimension D via box counting on a pairwise distance matrix @liebovitch1989fast. This approach suffered from two major drawbacks: it was difficult to optimise without extensive refactoring, and it required binarising the distance matrix using an arbitrary similarity threshold. For this reason, the correlation dimension D2 was implemented using the Grassberger-Procaccia algorithm @grassberger1983characterization, producing similar self-similarity results with much faster runtime and no free parameters. Whereas box counting measures how the number of points scales with box size, D2 measures how the number of point pairs within a given distance scales with increasing radius. Even so, simulations still required 2-3 days per run, necessitating a different approach.

A novel self-similarity algorithm was therefore developed using a k-mer-based Hamming distance approximation. Each repeat was converted into a 256-bit presence-absence fingerprint of all possible 4-mers, and pairwise distances were computed using the `__builtin_popcountll()` intrinsic on the XOR of fingerprint pairs. This approach substantially outperformed both previous methods, reducing simulation time to ~3 hours with self-similarity calculated at each generation. The software was packaged as a Python library and published on GitHub (#link("https://github.com/mbeavitt/kmer-variance/")).

To further optimise performance, an AVX2-compatible implementation was used, collapsing four 64-bit XOR operations into a single 256-bit vector XOR. Although AVX2 lacks a native 256-bit popcount instruction, recent versions of Clang automatically apply Muła's PSHUFB-based algorithm @mula2018faster, computing the popcount entirely within vector registers. These optimisations yielded a modest 1.67× speed-up, reducing runtime to ~2 hours. While the actual time savings here were not substantial, the experience of learning to apply hardware-level optimisation was valuable, particularly given the dramatic speed-ups achieved by AVX-accelerated bioinformatics tools such as BWA-MEM @8820962.

The optimised simulation framework, incorporating all efficiency improvements and the novel k-mer-based self-similarity algorithm, was packaged as a Python tool called censim for ease of use and reproducibility (#link("https://github.com/mbeavitt/censim")).

== Model implementation

Two simulation models were implemented and compared: the Uniform Recombination model (baseline) and the Kinetochore-associated Recombination model (KARMA). Both models share the same mutation rates but differ in how mutation positions are selected.

In the *Uniform Recombination model*, all mutation events occur uniformly across the array with probabilities per generation:

$ "SNPs" &arrow.r "Poisson"(lambda=0.1) $
$ "Duplications" &arrow.r "Poisson"(lambda=0.25), "size" = "Poisson"(lambda=7.6) times 178 $
$ "Deletions" &arrow.r "Poisson"(lambda=0.25), "size" = "Poisson"(lambda=7.6) times 178 $

All mutation positions are drawn uniformly from the array:

$ "Position" tilde "Uniform"(0, L-1) $

where $L$ is the current array length in repeat units.

In the *Kinetochore-associated Recombination model*, SNPs continue to occur uniformly, but insertions and deletions are weighted by the local self-similarity of the array. The self-similarity at position $i$ is quantified using the correlation dimension $D_2 (i)$, and mutation positions are sampled from a categorical distribution with weights proportional to $D_2 (i)$:

$ "SNPs" &arrow.r "Uniform"(0, L-1) $
$ "Duplications" &arrow.r "Categorical"({0, ..., L-1}, w_i = D_2 (i)) $
$ "Deletions" &arrow.r "Categorical"({0, ..., L-1}, w_i = D_2 (i)) $

The mutation rates per generation remain identical to the Uniform model:

$ "SNPs" &arrow.r "Poisson"(lambda=0.1) $
$ "Duplications" &arrow.r "Poisson"(lambda=0.25), "size" = "Poisson"(lambda=7.6) times 178 $
$ "Deletions" &arrow.r "Poisson"(lambda=0.25), "size" = "Poisson"(lambda=7.6) times 178 $

This modification captures the KARMA hypothesis that kinetochore-associated recombination machinery preferentially targets self-similar regions, leading to increased turnover in these areas and potentially driving homogenisation of satellite arrays.

== Running the simulations

400 simulations were run for each group (uniform recombination and kinetochore-associated recombination) on a Linux computer with a 30 core CPU (Ryzen 9 5950x) and 64GB RAM. As we were expecting roughly 15% of the runs to collapse due to random contraction and expansion dynamics, the target number of 330 uncollapsed repeat arrays of each group necessitated approximately 400 runs. The simulations were run in parallel, with one core each and 2GB of RAM. There were no memory constraints to speak of as the main memory bottleneck was the repeat array, taking up only 2-10MB of RAM.

The starting parameters of the simulation were kept essentially the same as in Dong et al's model, using 15,000×178bp repeats as the starting array state, the starting repeat template being the same _A. thaliana_ CEN178 consensus sequence -- "AGTATAAGAACTTAAACCGCAACCCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAAGACACAAAGCCAAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAGAAGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGATGTCATGTGTATGATTG" -- and the mutation frequency per generation parameters remaining the same for each mutation type (_SNP λ=0.1, INDEL λ=0.5; Poisson_) with the INDEL size similarly drawn from a Poisson distribution (_λ=7.6_). The array simulations proceeded with checkpointing every million generations until they reached the six-millionth generation, or until the array collapsed to a threshold of 2000 repeats at which point the simulation terminated.

TRASH @wlodzimierz2023trash was used to predict HORs in the simulated data (HOR metrics already existed for the 66 Arabidopsis accessions and were obtained from the authors of the publication). Since we already knew the location and size of all of the repeats in each simulated centromere, a table of repeats was composed from the simulated centromere sequences and used as input to only the HOR identification module of TRASH: HORT.R. The module was modified to use mafft in low memory mode, and run in parallel using GNU parallel @tange_2025_17692695 on a machine with a 30 core CPU (Ryzen 9 5950x) and 64GB RAM with flags `--hor_threshold=3`, `--min_hor_len=3`.

In comparing HORs/kb, the simulated data groups were approximately log-normal in distribution with different variances, so Welch's t-test was applied on the log-transformed data to compare means between the Uniform Recombination and Kinetochore-associated Recombination groups. The _A. thaliana_ group exhibited substantially higher variance and greater deviation from log-normality compared to the simulated groups, precluding direct parametric comparison.

== Metric calculations

Block size was measured as the number of repeat units contained within each block. For non-overlapping HORs, block gap was computed as the distance between the end of block A and the start of block B; overlapping configurations were assigned a gap value of zero.

I assessed the organizational quality of each block by calculating the ratio of the number of unique monomer sequences to block size. Blocks with higher ratios indicate greater organizational diversity, where more distinct sequences contribute to the block structure. The average block quality (unique monomers per unit) was computed as:

$ Q_"avg" = 1/2 (U_A/n_A + U_B/n_B) $ <eq:block-quality>

where $n_X$ is the block size and $U_X$ is the number of unique monomer sequences in the corresponding block.

Inter-block similarity ($S_"HOR"$) was quantified by computing position-wise sequence distances between aligned repeat units in blocks A and B. For each position $i in [1, n]$, the Levenshtein edit distance $d(s_A^i, s_B^i)$ between the sequences at position $i$ in block A and block B was calculated. The mean positional distance was then transformed to a similarity score:

$ S_"HOR" = 1/(1 + macron(d)_"pos") $ <eq:hor-similarity>

where

$ macron(d)_"pos" = 1/n sum_(i=1)^n d(s_A^i, s_B^i) $ <eq:mean-distance>

This metric ranges from 0 to 1, where 1 indicates perfectly identical blocks, and the metric approaches 0 as the blocks' pairwise repeats diverge in terms of their sequences.

Each of the above metrics plus relevant metrics extracted from the TRASH analysis (block size, block gap) was compiled per centromere into a table, giving 1045 HOR tables. Due to the large total size of the tables (~500 GB), distributions were computed by first calculating fixed-width histogram bins per HOR table to facilitate plotting in Python without loading the entire contents of each group's tables into memory:

#figure(
  table(
    columns: 5,
    align: (left, center, left, center, left),
    [*Variable*], [*Number of Bins*], [*Range*], [*Bin Width*], [*Binning Type*],
    [Internal Diversity], [2000], [0 - 104.0], [0.052], [Continuous],
    [Block Size (units)], [2532], [3 - 2534], [1], [Discrete],
    [Block Gap (units)], [400], [0 - 35,477.4], [88.69], [Continuous],
    [Unique Monomers per Unit], [200], [0 - 1.05], [0.0053], [Continuous],
    [HOR similarity], [200], [0 - 1.05], [0.0053], [Continuous],
    [Composite HOR metric], [400], [0 - 57,499,686.6], [143,749.2], [Continuous],
  ),
  caption: [Binning parameters for HOR metric distributions]
) <tab:binning>

= Results

Out of the 400 simulations run for each of the uniform recombination and kinetochore-linked recombination models, 352 and 363 simulations, respectively, reached the six millionth generation without collapsing to zero by random chance. As a result, the distributions of repeat number counts for the simulated centromeres are right-skewed, reflecting a hard lower bound at zero repeats and the early termination of simulations that approach this boundary. Using similarity heatmap plots inspired by those in use in the tool ModDotPlot @sweeten2024moddotplot to inspect the repeat array structures at a coarse-grained level at each checkpoint of a million generations, we see self-similar patches resembling those observed in real centromere arrays @naish2021genetic @wlodzimierz2023cycles @sweeten2024moddotplot as previously observed by Dong et al (@fig:supp-heatmaps, Supplementary).

The true centromeric repeat arrays (@fig:faceted-histograms) had both a higher number of unique repeats per kilobase on average, as well as higher variance in the number of unique repeats, with each of the distributions appearing approximately normally distributed. This suggests that neither of the models fully capture the observed mean/variance structure of unique repeats per kilobase in real repeat arrays. The mean sizes of the repeat arrays were approximately similar, but the variance was much greater in both of the simulated arrays. Due to the nature of the simulations (random expansion/contraction) it's unlikely that the higher mean of the Kinetochore-associated Recombination has any relation to the modifications made to this model, but the differing variances are of interest, suggesting that the real arrays' evolution has a size-constraining mechanism.

#figure(
  image("faceted_histograms.png", width: 100%),
  caption: [
    *Centromeric repeat diversity comparison of simulated arrays and A. thaliana arrays.*
    *A:* Unique repeats per kilobase across three groups of repeat arrays: Kinetochore-associated Recombination model (blue, n=363), Uniform Recombination model (grey, n=352), and _A. thaliana_ centromeric arrays (red, n=330). _A. thaliana_ arrays are more diverse with higher variance (μ=2.3 unique repeats/kb, σ²=0.080) compared to both kinetochore-associated (μ=1.1, σ²=0.014) and uniform recombination models (μ=1.5, σ²=0.006).
    *B:* Distribution of total repeat counts per array. _A. thaliana_ repeat arrays contain fewer total repeats and lower variance (μ=15,461, σ²=2.0×10⁷) than arrays in either kinetochore-associated (μ=29,640, σ²=1.7×10⁸) or uniform recombination simulations (μ=20,173, σ²=9.8×10⁷).
  ]
) <fig:faceted-histograms>

The TRASH software was used to investigate the structure of the repeat arrays by identifying higher order repeats (HORs), identifying millions of HORs per group (Kinetochore-associated model: 380,319,638, Uniform model: 160,562,751, _A. thaliana_: 340,782,660). The number of HORs per kb appeared to differ between the three groups, and the variance differed greatly between both simulations and the HORs/kb in the real array (@fig:num-hors). The mean HORs/kb in the Kinetochore-associated Recombination model was significantly greater than the Uniform Recombination model (Welch's t-test, p < 0.001) but apparently lower than the real _A. thaliana_ arrays (statistical test not applied). The high variance of the _A. thaliana_ arrays precluded meaningful statistical comparison with the simulated arrays, but was of interest in and of itself -- this suggests, as before, that the generative processes that produce the repeat arrays seen in nature are not being fully captured by either of the models.

#figure(
  image("num_hors_histogram.png", width: 100%),
  caption: [
    *Higher order repeat density between simulated and real repeat arrays (log scale).*
    HORs per kilobase were measured across three groups of repeat arrays: Kinetochore-associated Recombination model (blue, n=363), Uniform Recombination model (grey, n=352), and _A. thaliana_ centromeric arrays (red, n=330). _A. thaliana_ had more HORs/kb on average (μ=261.4, σ²=97422.33) than either of the simulated datasets as well as a wider spread of HOR density values. The Kinetochore-associated Recombination model (μ=192.6, σ²=4695.82) exhibited significantly more HORs/kb than the Uniform Recombination model (μ=124.5, σ²=1678.62), Welch's t-test p < 0.001. X-axis displayed on logarithmic scale.
  ]
) <fig:num-hors>

A HOR is defined as an alignment of two blocks of similar repeats separated by any distance, with no two repeats in alignment having an edit distance greater than 3 and at least three repeats in the block. To further investigate the structural differences between the higher order repeats in the three groups, the unique monomers per unit metric was calculated and plotted. This metric measures, across every identified pair of HOR blocks, how diverse each block is internally. Blocks composed entirely of unique repeats (e.g. ABCDEF … ABCDEF, where A-F are different unique satellite repeats) have a score of 1, and blocks composed entirely of the same repeat (e.g. AAAAAA … AAAAAA, where 'A' is the same satellite repeat) have a score of 1/n, where n is the length of the block. The total number of HORs differed between each group, and so a cumulative distribution of the frequencies of unique monomers per unit in each group was plotted. The HOR diversity in _A. thaliana_ repeat arrays was observed to be very high with ~95% of HORs having a diversity score of 1.0 -- that is, all of the repeats in the block are completely unique. In contrast, only ~65% of repeats had a diversity score of 1.0 in the Uniform Recombination model and ~35% in the Kinetochore-associated Recombination model.

#figure(
  image("unique_monomers_per_unit_cumulative.png", width: 100%),
  caption: [
    *Cumulative distributions of simulated and real repeat array HOR diversity frequencies.*
    Cumulative frequency distributions were plotted for each group of repeat arrays (Kinetochore-associated Recombination, Uniform Recombination, _A. thaliana_), showing the fraction of HORs with a ratio of unique monomers per repeat unit exceeding a given value.
  ]
) <fig:unique-monomers>

Comparing HOR block sizes directly using a log-log plot, we recapitulate the power law describing the prevalence of different HOR sizes which was identified in a recent publication, and identify that the simulated datasets exhibit approximate power-law scaling but are truncated, exhibiting fewer extremely high values than expected under a pure power law (@fig:block-size). In particular, the distribution describing the prevalence of block sizes in the simulated repeat arrays is typical of a truncated power law. Where a power law is described as $P(x) prop x^(-a)$, truncated power laws are described as $P(x) prop (x^(-alpha))(e^(-x/x_c))$ where $x_c$ is the cutoff scale and $e$ is euler's constant. In particular, this suggests that where the simulated models are the result of additive growth, the real HOR sizes are produced by multiplicative growth. As the simulation uses a Poisson model with a fixed mean and variance to describe the expansion of a given repeat, this could explain why we observe that HOR size appears to follow a truncated power law. A multiplicative effect (for example, that large HORs grow more quickly) might be an alternative model to consider in order to produce a distribution of HOR sizes closer to that observed in nature. The Kinetochore-associated Recombination model only linked the position of INDELs to the self-similarity at a given location in the array -- it would be of interest to understand if additionally linking the INDEL size to self-similarity would reproduce a non-truncated power law distribution.

#figure(
  image("block.size.in.units_combined.png", width: 100%),
  caption: [
    *Comparison of distributions of HOR block size between the simulated and real repeat arrays (log-log scale).*
    *A:* Density plots describing the distribution of HOR block sizes in each of the groups. Dotted lines indicate median values of each group. Both the X-axis and Y-axis displayed on logarithmic scale.
    *B, C, D:* Histograms displaying the distribution of HOR block sizes in each group independently. Each have similar median values (either 3 or 4 blocks), but _A. thaliana_ has fewer very large HORs in general (99th percentile=10) while the simulated datasets are more likely to have very large HORs (Kinetochore-associated Recombination: 99th percentile = 111, Uniform Recombination: 99th percentile=85). Both the X-axis and Y-axis displayed on logarithmic scale.
  ]
) <fig:block-size>

In contrast, the HOR block gap distributions (@fig:block-gap) appear to follow a typical logarithmic distribution where large gaps become increasingly rare. Comparing the simulations and the _A. thaliana_ repeat arrays, it's clear that _A. thaliana_ has a higher frequency of HOR blocks with a moderately large gap, as well as more HORs with extremely large gaps as indicated by the 90th/99th percentiles (@fig:block-gap, D). In addition, it has fewer HORs with a small gap than either of the simulated datasets. Interestingly, the Kinetochore-associated Recombination model arrays show fewer small-gap HORs and more large-gap HORs (@fig:block-gap, B) than the Uniform Recombination model (@fig:block-gap, D).

#figure(
  image("block_gap_units_combined.png", width: 100%),
  caption: [
    *Comparison of distributions of HOR block gaps between the simulated and real repeat arrays (log scale).*
    *A:* Density plots describing the distribution of HOR block gaps in each of the groups. Dotted lines indicate median values of each group. X-axis displayed on logarithmic scale.
    *B, C, D:* Histograms displaying the distribution of HOR block gaps in each group independently. X-axis displayed on logarithmic scale.
  ]
) <fig:block-gap>

The final metric compared was HOR similarity, which compared the degree of concordance between the two blocks of the HOR using the mean Levenshtein distance between all pairwise repeats in the alignment -- a score of 1 indicates total identity between the two blocks, and the score approaches zero as they diverge (up to a cutoff determined by the value of `--hor_threshold` in the HOR identification program TRASH). Plotting the density distribution of these values in each group highlighted that they had remarkable similarity to each other, with the Uniform Recombination model and _A. thaliana_ displaying identical median and 90/99th percentile values as well as a highly similar histogram. The Kinetochore-associated Recombination model displayed a notable excess of high HOR similarity values and a paucity of low ones compared to the other two.

#figure(
  image("hor_similarity_combined.png", width: 100%),
  caption: [
    *Comparison of distributions of HOR similarity scores between the simulated and real repeat arrays.*
    *A:* Density plots describing the density of HOR similarity values across each of the groups. Dotted lines indicate median values of each group.
    *B, C, D:* Histograms displaying the distribution of HOR similarity values in each group independently. _A. thaliana_ and the Uniform Recombination model's data have exactly the same Median, 90th percentile and 99th percentile (0.16, 0.23, 0.40 respectively).
  ]
) <fig:hor-similarity>

Finally, a composite metric consisting of the product of the block gap, HOR similarity, block size and HOR diversity was calculated in an attempt to capture and showcase the concept that real _A. thaliana_ repeat arrays exhibit large, highly similar, internally diverse and widely spaced higher order repeats -- each of the models are able to move closer to real repeat arrays in terms of one or more of these metrics, but neither of them approach _A. thaliana_ repeat array organisation in terms of all four, making it a useful metric to gauge overall model performance.

#figure(
  image("composite_hor_metric_overlaid.png", width: 100%),
  caption: [
    *Composite HOR metric.*
    The composite HOR metric density was plotted for each group, representing the multiplicative combined influences of HOR block gap, diversity, similarity and size. The two simulated datasets are similar to each other in terms of the shape of their distribution, while the _A. thaliana_ dataset has orders of magnitude higher scores in this metric, highlighting its potential use as an overall quality metric in future simulations.
  ]
) <fig:composite-metric>

A useful feature of this composite metric is to be able to filter real centromeric data to discover large and unexpected long-range repeat organisation, likely due to the occurrence of long range and high magnitude duplications during crossover mutation events.

= Discussion

The results of the HOR similarity and HOR diversity calculations clearly indicate that the differences in block gap distributions and HORs per kb in the Kinetochore-associated Recombination model compared with the Uniform Recombination model are in no small part due to a simple excess of repeat similarity as a result of targeting INDEL events to self-similar regions. It's clear that the fundamental generative processes underpinning repeat diversity are correctly captured by a simple INDEL/SNP model, at least in terms of HOR metrics, but that the appearance of HORs with high intra-block diversity yet low inter-block diversity, along with a large block gap and block size cannot be explained by this simple model. In simpler terms, the model goes some way to explain accurately the evolution of repeat sequences themselves but not the organisation of these repeats along the centromere.

The other primary phenomenon that goes unexplained by the simple Uniform Recombination model is the deviation of HOR unit sizes from a power law distribution, forming instead a truncated power law distribution, suggesting there is a missing multiplicative mechanism tying together the size of HOR units and the magnitude of expansion events. Future investigations modelling expansion events as a function of self-similarity may go some way in testing this hypothesis definitively.

Both a random mutation model with a balance of diversifying and homogenising mutations and a directed model with recombination machinery targeting likely kinetochore attachment regions are able to produce centromere-like arrays, but it's not possible to conclusively support one model over the other in this simulation-based approach. What is supported, however, is the validity of the kinetochore-associated recombination model in comparison to the uniform recombination model, as both produce arrays with centromere-like characteristics. That said, it's certainly true that the driving mechanism that produces self-similar satellite "patches" in the simulations (and apparently in biological centromere arrays) is balanced, frame-preserving insertion-deletion events combined with single nucleotide polymorphisms.

Given the high rates of transposon insertion events and their targeting to centromere-specific proteins @tsukahara2025centrophilic, it's clear that their contribution to the evolution of centromeres cannot be ignored. Further simulation-based experiments would benefit from including high frequency transposon insertions. If their addition precludes the emergence of HORs in the uniform recombination or kinetochore-associated recombination models, this would suggest a biological mechanism for purging these transposable elements from the arrays. Unfortunately, testing this hypothesis was not possible within the time constraints of the rotation project due to the added complexity of modelling frame-disrupting insertions in our wholly in-frame mutation paradigm.

Another failing of the current approach is that we ignore the effects of centromere drive, a mechanistic theory of centromere evolution popularised by Henikoff et al in 2001 @henikoff2001centromere. This model proposes that asymmetric female meiosis might impose a selection pressure on centromeres composed of tandemly-repeated satellite arrays, with repeat arrays having characteristics that allow them to recruit more kinetochore proteins such as greater homogeneity or larger repeat arrays or some combination of both of these characteristics result in their selection over another competing centromere array. It's conceivable that such a model could explain the surprising homogeneity of centromere arrays in isolation, as well as in addition to the molecular turnover models presented here and in Dong et al.

Finally, the usability of a Monte Carlo simulation involving the shuffling and editing of DNA represented as ASCII strings in memory is limited in exploring more creative mechanistic models. The necessary overheads of inserting/deleting elements from large arrays millions of times per simulation impose a substantial calculation burden, limiting the types of mechanisms that can be modelled in reasonable computational timescales. Modelling and testing other moving parts in combination with the current model for example CENH3 positions, methylation, centromere drive, or including the rest of the genome in the evolution (for example it was shown that a centromere allele linked to centromere drive was not fixed in the population due to linked deleterious effects inherited with the allele @fishman2015centromere) would require more and more time to be devoted to optimising the simulation. By instead capturing the high-level mechanistic effects of centromere turnover rather than trying to discretely model individual mutation events would resolve this issue -- either by using an existing, heavily optimised evolutionary framework such as SLiM @haller2025slim or by turning to numerical dynamical systems as representations of sequences and the relationships between them. Due to the complex nature of representing the mutation and inheritance patterns of centromere arrays as discrete alleles, it's likely that the latter approach would be better suited to further more complex simulation experiments.

In planning the project and the order of the experiments, in hindsight it would have been better to start with comparing the results of the uniform recombination model with the existing natural centromeres in the 66 assembled Arabidopsis accessions to understand its shortcomings more completely before implementing modifications to the model. On the other hand, the kinetochore-associated recombination model is widely cited by numerous groups studying centromere evolution @gao2025global @naish2024structure, and a computational model of centromeric repeats capturing this had not yet been implemented, making this an interesting angle to consider in its own right.

Overall, the results of this work demonstrate that a model capturing the dynamics and end results of centromere evolution should be based on comparatively frequent INDELs and low frequency SNPs, compared to chromosome arms, but that it is clear that long-range and large-scale duplications are an inherent feature of natural centromere arrays that are identifiable as higher order repeats. The models proposed here and by Dong et al. produce higher order repeats but the repeat blocks are either too self-similar, not far enough apart, too small, or not diverse enough internally. In other words, large tandem duplications occur locally but are quickly broken up by other duplications and mutations. A model either modelling a much higher frequency of tandem duplications over long ranges or a mechanism for reducing the effects of their degradation would go some way in explaining the mutational dynamics of eukaryotic centromere satellite repeats.

The effect of transposon invasions also cannot be ignored, as it is a strong and potentially driving mutational force in centromere arrays -- the introduction of these events in the current simulation paradigm is not straightforward and may necessitate alternate methods.

= Supplementary Figures

#pagebreak()

// Landscape page for supplementary figure
#set page(
  flipped: true,
  margin: (x: 2cm, y: 2cm),
)

#figure(
  image("combined_facet.png", width: 100%),
  caption: [
    *Self-similarity heatmaps.*
    The self-similarity heatmaps of each generation were plotted where dark blue represents high similarity between two points and bright yellow represents low similarity between two points. The plot below shows the value of the kmer similarity index metric, used to add the "KARMA" feedback loop addition to the simulation model and as an indicator of self-similarity hotspots.
  ]
) <fig:supp-heatmaps>

#pagebreak()

// Reset to portrait
#set page(
  flipped: false,
  margin: (x: 2.5cm, y: 2.5cm),
)

#bibliography("references.bib")
