# Understanding the evolutionary landscape of centromere evolution in Arabidopsis using monte carlo simulations

Name: Michael Beavitt  
Email: [mab282@cam.ac.uk](mailto:mab282@cam.ac.uk)  
Department: Department of Plant Sciences  
Supervisor: Professor Ian Henderson

Word count:

Statement:

Date: 05/01/2026

# Introduction

Eukaryotic centromeres are essential components which allow for accurate segregation of chromosomes and chromatids during cell division. The relationship between the DNA composition of centromere sequences and their function has been long questioned \- so much so, that it has been described as the “centromere paradox” (henikoff2001centromere). While it’s clear that the inheritance of centromeric function is primarily epigenetic, there still remains the enigma of the abundance of satellite repeat arrays across the eukaryotes. Do the satellite repeats themselves confer some kind of fitness to the centromere, making them more likely to be inherited during meiosis? If satellite repeats appear at centromeric loci, are the molecular components of cell division themselves \- the kinetochore \- responsible for generating satellite repeats somehow? With recent advancements in long read sequencing technologies, gap-free centromere sequences in Arabidopsis have been resolved for the first time, giving us insights into their composition (Naish 2021, Wlodzimierz 2023). By using these resources to turn our attention to satellite array composition and understand how they might have come to be, perhaps we can begin to understand why they exist.  
   
While the mechanisms and drive of centromere evolution are uncertain, it is undisputed that turnover is much faster in centromere arrays than the rest of the genome, not only in the case of unusual frame-preserving INDELs but also single nucleotide polymorphisms (SNPs) which are nine times more frequent than the arms of the chromosomes, occurring at a rate of 6.12 × 10−8 per bp per generation (Dong preprint). Large inversions and tandem duplications have also been observed (Wlodzimierz 2023). These dynamics are not just limited to Arabidopsis, but occur across Eukaryotic lineages (Henderson/Wlodzimierz 2025 preprint). Transposable element (TE) insertion is also highly frequent in centromeric arrays, with concrete evidence for direct targeting to centromere-specific proteins (Tsukahara et al 2025).  
Many purport that the appearance of higher order repeat (HOR) sequences and centromere identity is a direct result of feedback mechanisms involving the kinetochore machinery, but the recent simulation work of Dong et al challenges this, showing that high frequency randomly inserted/deleted sequences combined with low frequency random single nucleotide polymorphisms is sufficient to generate centromere-like repeat arrays. This finding is surprising, given that we know that at least one driver of evolution is directly related to kinetochore machinery given Tsukahara et al’s finding that transposon insertions are targeted to CENH3 chromatin specifically, which exclusively occupies subsets of centromeric DNA.

What I wanted to do in this rotation project was to understand how the model proposed by Dong et al truly captures the structures and sequence diversity of repeat arrays in real genomes, and additionally modify the existing model to capture a hypothesised mechanism in eukaryote genomes \- kinetochore-associated recombination. While the original model produced intriguing results, only five centromeres were simulated. Large scale comparisons of hundreds of simulated centromere satellite arrays to their real counterparts would give new insights into the strengths and limitations of the model, and potentially help us to understand the true mechanisms at play. The recent publication of 66 *Arabidopsis thaliana* genomes, containing 330 well-resolved centromere sequences (wlodzimierz 2023\) provided a clear route forward to extend the analysis.

As touch upon already, this project also aimed to build upon the existing simulation framework by applying mechanisms inspired by the “KARMA” model (Kinetochore-associated Recombination Machinery in Arabidopsis) in an attempt to both understand its feasibility and explore the consequences of accepting the assumptions it represents.  
In brief, the KARMA model hypothesises that kinetochore protein complexes attract recombination factors which increase the turnover of satellite DNA in close proximity to the kinetochore itself. This, in theory, homogenises the array, purging non-similar repeats or invading transposon sequences and producing the highly homogeneous satellite arrays seen in nature. The motivation behind applying it to the simulation framework was to understand the implications of accepting the assumptions of this popular model \- would they produce unnatural or incorrect centromeric repeats, and if so in what manner? On the other hand, would adding this positive feedback loop, reinforcing self-similar regions of satellite DNA produce more “centromeric” sequences than the existing model?

# Methods

## Code Optimisation

The original Dong et al paper’s simulation code, published on github ([https://github.com/schneebergerlab/replicated-assemblies-centromere-study](https://github.com/schneebergerlab/replicated-assemblies-centromere-study)) provided a comprehensive starting point for simulating the three desired events \- point mutations, insertions and deletions. The original program unfortunately took \~3 days to complete a single 6 million generation run, which meant running the desired \~600 simulations would use approximately 43200 CPU hours or almost 5 years in real time on a single core. Multithreading would go some way to making this time constraint more manageable, but a better approach would involve optimising the core program as much as possible before applying multithreading.

Optimisations initially involved removing unneeded checks and calculations. In particular, the original program performed complex bounds checking to ensure each of the 178 bp repeats was indexed correctly, tracking repeat positions in a separate table. This was only necessary if repeat length could change, but all mutations modelled (SNPs, insertions and deletions in multiples of 178 bp) were frame-preserving, so this situation never arose. I therefore replaced the bounds checking with simple modulo arithmetic to index each repeat. Combined with replacing the most expensive array calculations with vectorized NumPy (harris2020array) operations, this reduced runtime to \~20 minutes.

Tracking each base individually remained a major bottleneck, with the array size scaling as 178 × repeat number (e.g. 2.67 million characters for 15,000 repeats). This made O(n) insertion and deletion operations extremely costly, so an abstracted repeat-resolution representation was implemented to enable more efficient list operations. This reduced runtime further to \~3 minutes, making hundreds of simulations feasible.

Array self-similarity was initially calculated using the fractal dimension D via box counting on a pairwise distance matrix (liebovitch1989fast). This approach suffered from two major drawbacks: it was difficult to optimise without extensive refactoring, and it required binarising the distance matrix using an arbitrary similarity threshold. For this reason, the correlation dimension D2 was implemented using the Grassberger-Procaccia algorithm (grassberger1983characterization), producing similar self-similarity results with much faster runtime and no free parameters. Whereas box counting measures how the number of points scales with box size, D2 measures how the number of point pairs within a given distance scales with increasing radius. Even so, simulations still required 2-3 days per run, necessitating a different approach.

A novel self-similarity algorithm was therefore developed using a k-mer-based Hamming distance approximation. Each repeat was converted into a 256-bit presence-absence fingerprint of all possible 4-mers, and pairwise distances were computed using the \_\_builtin\_popcountll() intrinsic on the XOR of fingerprint pairs. This approach substantially outperformed both previous methods, reducing simulation time to \~3 hours with self-similarity calculated at each generation. The software was packaged as a Python library and published on GitHub ([https://github.com/mbeavitt/kmer-variance/](https://github.com/mbeavitt/kmer-variance/)).

To further optimise performance, an AVX2-compatible implementation was used, collapsing four 64-bit XOR operations into a single 256-bit vector XOR. Although AVX2 lacks a native 256-bit popcount instruction, recent versions of Clang automatically apply Muła’s PSHUFB-based algorithm (Muła, Kurz & Lemire, 2018), computing the popcount entirely within vector registers. These optimisations yielded a modest 1.67× speed-up, reducing runtime to \~2 hours. While the actual time savings here were not substantial, the experience of learning to apply hardware-level optimisation was valuable, particularly given the dramatic speed-ups achieved by AVX-accelerated bioinformatics tools such as BWA-MEM (Vasimuddin et al. 2019).

## Running the simulations

400 simulations were run for each group (uniform recombination and kinetochore-associated recombination) on a Linux computer with a 30 core CPU (Ryzen 9 5950x) and 64GB RAM. As we were expecting roughly 15% of the runs to collapse due to random contraction and expansion dynamics, the target number of 330 uncollapsed repeat arrays of each group necessitated approximately 400 runs. The simulations were run in parallel, with one core each and 2GB of RAM. There were no memory constraints to speak of as the main memory bottleneck was the repeat array, taking up only 2-10MB of RAM. 

The starting parameters of the simulation were kept essentially the same as in Dong et al’s model, using 15,000x178bp repeats as the starting array state, the starting repeat template being the same *A. thaliana* CEN178 consensus sequence \- "AGTATAAGAACTTAAACCGCAACCCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAAGACACAAAGCCAAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAGAAGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGATGTCATGTGTATGATTG" \- and the mutation frequency per generation parameters remaining the same for each mutation type (*SNP λ=0.1, INDEL λ=0.5; Poisson*) with the INDEL size similarly drawn from a Poisson distribution (*λ=7.6)* (dong et al publication). The array simulations proceeded with checkpointing every million generations until they reached the six-millionth generation, or until the array collapsed to a threshold of 2000 repeats at which point the simulation terminated.

TRASH (wlodzimierz2023trash) was used to predict HORs in the simulated data (HOR metrics already existed for the 66 Arabidopsis accessions and were obtained from the authors of the publication). Since we already knew the location and size of all of the repeats in each simulated centromere, a table of repeats was composed from the simulated centromere sequences and used as input to only the HOR identification module of TRASH: HORT.R. The module was modified to use mafft in low memory mode, and run in parallel using GNU parallel (tange\_2025\_17692695) on a machine with a 30 core CPU (Ryzen 9 5950x) and 64GB RAM with flags \--hor\_threshold=3, \--min\_hor\_len=3.

In comparing HORs/kb, the simulated data groups were approximately log-normal in distribution with different variances, so Welch’s t-test was applied on the log-transformed data to compare means between the Uniform Recombination and Kinetochore-associated Recombination groups. The *A. thaliana* group exhibited substantially higher variance and greater deviation from log-normality compared to the simulated groups, precluding direct parametric comparison.

## Metric calculations

Block size was measured as the number of repeat units contained within each block. For non-overlapping HORs, block gap was computed as the distance between the end of block A and the start of block B; overlapping configurations were assigned a gap value of zero.

I assessed the organizational quality of each block by calculating the ratio of the number of unique monomer sequences to block size. Blocks with higher ratios indicate greater organizational diversity, where more distinct sequences contribute to the block structure. The average block quality (unique monomers per unit) was computed as:

$$Q\_{\\text{avg}} \= \\frac{1}{2}\\left(\\frac{|U\_A|}{n\_A} \+ \\frac{|U\_B|}{n\_B}\\right)$$

where $$n\_X$$ is the block size and $$|U\_X|$$ is the number of unique monomer sequences in the corresponding block.

Inter-block similarity ($$S\_{\\text{HOR}}$$) was quantified by computing position-wise sequence distances between aligned repeat units in blocks A and B. For each position $$i \\in \[1, n\]$$, the Levenshtein edit distance $$d(s\_A^i, s\_B^i)$$ between the sequences at position $$i$$ in block A and block B was calculated. The mean positional distance was then transformed to a similarity score:   
$$S\_{\\text{HOR}} \= \\frac{1}{1 \+ \\bar{d}{\\text{pos}}}$$ 

where 

$$\\bar{d}{\\text{pos}} \= \\frac{1}{n} \\sum\_{i=1}^{n} d(s\_A^i, s\_B^i)$$ 

This metric ranges from 0 to 1, where 1 indicates perfectly identical blocks, and the metric approaches 0 as the blocks’ pairwise repeats diverge in terms of their sequences.

Each of the above metrics plus relevant metrics extracted from the TRASH analysis (block size, block gap,  was compiled per centromere into a table, giving 1045 HOR tables. Due to the large total size of the tables (\~500 GB), distributions were computed by first calculating fixed-width histogram bins per HOR table to facilitate plotting in python without loading the entire contents of each group’s tables into memory:

| Variable | Number of Bins | Range | Bin Width | Binning Type |
| :---- | ----- | :---- | ----- | :---- |
| Internal Diversity | 2000 | 0 \- 104.0 | 0.052 | Continuous |
| Block Size (units) | 2532 | 3 \- 2534 | 1 | Discrete |
| Block Gap (units) | 400 | 0 \- 35,477.4 | 88.69 | Continuous |
| Unique Monomers per Unit | 200 | 0 \- 1.05 | 0.0053 | Continuous |
| HOR similarity | 200 | 0 \- 1.05 | 0.0053 | Continuous |
| Composite HOR metric | 400 | 0 \- 57,499,686.6 | 143,749.2 | Continuous |

# Results

Out of the 400 simulations run for each of the uniform recombination and kinetochore-linked recombination models, 352 and 363 simulations, respectively, reached the six millionth generation without collapsing to zero by random chance. As a result, the distributions of repeat number counts for the simulated centromeres are right-skewed, reflecting a hard lower bound at zero repeats and the early termination of simulations that approach this boundary. Using similarity heatmap plots inspired by those in use in the tool ModDotPlot (sweeten2024moddotplot) to inspect the repeat array structures at a coarse-grained level at each checkpoint of a million generations, we see self-similar patches resembling those observed in real centromere arrays (naish 2021, wlodzimeirz 2023 not the TRASH paper but the other one, sweeten2024moddotplot) as previously observed by dong et el (dong preprint) (Figure S1)

The true centromeric repeat arrays (Figure X (faceted\_histograms.png)) had both a higher number of unique repeats per kilobase on average, as well as higher variance in the number of unique repeats, with each of the distributions appearing approximately normally distributed. This suggests that neither of the models fully capture the observed mean/variance structure of unique repeats per kilobase in real repeat arrays. The mean sizes of the repeat arrays were approximately similar, but the variance was much greater in both of the simulated arrays. Due to the nature of the simulations (random expansion/contraction) it’s unlikely that the higher mean of the Kinetochore-associated Recombination has any relation to the modifications made to this model, but the differing variances are of interest, suggesting that either the real arrays’ evolution have a size-constraining mechanism.

![][image1]  
**Figure X (faceted\_histograms.png). Centromeric repeat diversity comparison of simulated arrays and A. thaliana arrays**  
**A:** Unique repeats per kilobase across three groups of repeat arrays: Kinetochore-associated Recombination model (blue, n=363), Uniform Recombination model (grey, n=352), *and A. thaliana* centromeric arrays (red, n=330). *A. thaliana* arrays are more diverse with higher variance (μ=2.3 unique repeats/kb, σ²=0.080) compared to both kinetochore-associated (μ=1.1, σ²=0.014) and uniform recombination models (μ=1.5, σ²=0.006).  
**B:** Distribution of total repeat counts per array. A. thaliana repeat arrays contain fewer total repeats and lower variance (μ=15,461, σ²=2.0×10⁷) than arrays in either kinetochore-associated (μ=29,640, σ²=1.7×10⁸) or uniform recombination simulations (μ=20,173, σ²=9.8×10⁷)

The TRASH software was used to investigate the structure of the repeat arrays by identifying higher order repeats (HORs), identifying millions of HORs per group (Kinetochore-associated model: 380,319,638, Uniform model: 160,562,751, *A. thaliana*: 340,782,660). The number of HORs per kb appeared to differ between the three groups, and the variance differed greatly between the both simulations and the HORs/kb in the real array (Figure X (num\_hors\_histogram.png)). The mean HORs/kb in the Kinetochore-associated Recombination model was significantly greater than the Uniform Recombination model (Welch’s t-test, p \< 0.001) but apparently lower than the real A. thaliana arrays (statistical test not applied). The high variance of the A. thaliana arrays precluded meaningful statistical comparison with the simulated arrays, but was of interest in and of itself \- this suggests, as before, that the generative processes that produce the repeat arrays seen in nature are not being fully captured by either of the models.

![][image2]  
**Figure X (num\_hors\_histogram.png). Higher order repeat density between simulated and real repeat arrays (log scale)**  
HORs per kilobase were measured across three groups of repeat arrays: Kinetochore-associated Recombination model (blue, n=363), Uniform Recombination model (grey, n=352), *and A. thaliana* centromeric arrays (red, n=330). A. thaliana had more HORs/kb on average (μ=261.4, σ²=97422.33) than either of the simulated datasets as well as a wider spread of HOR density values. The Kinetochore-associated Recombination model (μ=192.6, σ²=4695.82) exhibited significantly more HORs/kb than the Uniform Recombination model (μ=124.5, σ²=1678.62), Welch’s t-test p \< 0.001. X-axis displayed on logarithmic scale.

A HOR is defined as an alignment of two blocks  of similar repeats separated by any distance, with no two repeats in alignment having an edit distance greater than 3 and at least three repeats in the block. To further investigate the structural differences between the higher order repeats in the three groups, the unique monomers per unit metric was calculated and plotted. This metric measures, across every identified pair of HOR blocks, how diverse each block is internally. Blocks composed entirely of unique repeats (e.g. ABCDEF … ABCDEF, where A-F are different unique satellite repeats) have a score of 1, and blocks composed entirely of the same repeat (e.g. AAAAAA … AAAAAA, where ‘A’ is the same satellite repeat) have a score of 1/n, where n is the length of the block. The total number of HORs differed between each group, and so a cumulative distribution of the frequencies of unique monomers per unit in each group was plotted. The HOR diversity in *A. thaliana* repeat arrays was observed to be very high with \~95% of HORs having a diversity score of 1.0 \- that is, all of the repeats in the block are completely unique. In contrast, only \~65% of repeats had a diversity score of 1.0 in the Uniform Recombination model and \~35% in the Kinetochore-associated Recombination model.

![][image3]  
**Figure X (unique\_monomers\_per\_unit\_cumulative.png). Cumulative distributions of simulated and real repeat array HOR diversity frequencies**  
Cumulative frequency distributions were plotted for each group of repeat arrays (Kinetochore-associated Recombination, Uniform Recombination, A. thaliana), showing the fraction of HORs with a ratio of unique monomers per repeat unit exceeding a given value.

Comparing HOR block sizes directly using a log-log plot, we recapitulate the power law describing the prevalence of different HOR sizes which was identified in a recent publication, and identify that the simulated datasets exhibit approximate power-law scaling but are truncated, exhibiting fewer extremely high values than expected under a pure power law (Figure X (block.size.in.units\_combined.png)). In particular, the distribution describing the prevalence of block sizes in the simulated repeat arrays is typical of a truncated power law. Where a power law is described as P(x)∝x^-a, truncated power laws are described as P(x) ∝ (x^−α)(e^−x/xc) where xc is the cutoff scale and e is euler’s constant. In particular, this suggests that where the simulated models are the result of additive growth, the real HOR sizes are produced by multiplicative growth. As the simulation uses a Poisson model with a fixed mean and variance to describe the expansion of a given repeat, this could explain why we observe that HOR size appears to follow a truncated power law. A multiplicative effect (for example, that large HORs grow more quickly) might be an alternative model to consider in order to produce a distribution of HOR sizes closer to that observed in nature. The Kinetochore-associated Recombination model only linked the position of INDELs to the self-similarity at a given location in the array \- it would be of interest to understand if additionally linking the INDEL size to self-similarity would reproduce a non-truncated power law distribution.  
![][image4]  
**Figure X (block.size.in.units\_combined.png). Comparison of distributions of HOR block size between the simulated and real repeat arrays (log-log scale)**  
**A:** Density plots describing the distribution of HOR block sizes in each of the groups. Dotted lines indicate median values of each group. Both the X-axis and Y-axis displayed on logarithmic scale.  
**B, C, D:** Histograms displaying the distribution of HOR block sizes in each group independently. Each have similar median values (either 3 or 4 blocks), but *A. thaliana* has fewer very large HORs in general (99th percentile=10) while the simulated datasets are more likely to have very large HORs (Kinetochore-associated Recombination: 99th percentile \= 111, Uniform Recombination: 99th percentile=85). Both the X-axis and Y-axis displayed on logarithmic scale.

In contrast, the HOR block gap distributions (Figure X (block\_gap\_units\_combined.png)) appear to follow a typical logarithmic distribution where large gaps become increasingly rare. Comparing the simulations and the *A. thaliana* repeat arrays, it’s clear that A. thaliana has a higher frequency of HOR blocks with a moderately large gap, as well as more HORs with extremely large gaps as indicated by the 90th/99th percentiles (Figure X (block\_gap\_units\_combined.png), D). In addition, it has fewer HORs with a small gap than either of the simulated datasets. Interestingly, the Kinetochore-associated Recombination model arrays show fewer small-gap HORs and more large-gap HORs (Figure X (block\_gap\_units\_combined.png), B) than the Uniform Recombination model (Figure X (block\_gap\_units\_combined.png), D).

![][image5]  
**Figure X (block\_gap\_units\_combined.png). Comparison of distributions of HOR block gaps between the simulated and real repeat arrays (log scale)**  
**A:** Density plots describing the distribution of HOR block sizes in each of the groups. Dotted lines indicate median values of each group. X-axis displayed on logarithmic scale.  
**B, C, D:** Histograms displaying the distribution of HOR block gaps in each group independently. X-axis displayed on logarithmic scale.

The final metric compared was HOR similarity, which compared the degree of concordance between the two blocks of the HOR using the mean Levenshtein distance between all pairwise repeats in the alignment \- a score of 1 indicates total identity between the two blocks, and the score approaches zero as they diverge (up to a cutoff determined by the value of \--hor\_threshold in the HOR identification program TRASH). Plotting the density distribution of these values in each group highlighted that they had remarkable similarity to each other, with the Uniform Recombination model and A. thaliana displaying identical median and 90/99th percentile values as well as a highly similar histogram. The Kinetochore-associated Recombination model displayed a notable excess of high HOR similarity values and a paucity of low ones compared to the other two.

![][image6]  
**Figure X (hor\_similarity\_combined.png). Comparison of distributions of HOR similarity scores between the simulated and real repeat arrays**  
**A:** Density plots describing the density of HOR similarity values across each of the groups. Dotted lines indicate median values of each group.  
**B, C, D:** Histograms displaying the distribution of  HOR similarity values in each group independently. A. thaliana and the Uniform Recombination model’s data have exactly the same Median, 90th percentile and 99th percentile (0.16, 0.23, 0.40 respectively).

Finally, a composite metric consisting of the product of the block gap, HOR similarity, block size and HOR diversity was calculated in an attempt to capture and showcase the concept that real *A. thaliana* repeat arrays exhibit large, highly similar, internally diverse and widely spaced higher order repeats \- each of the models are able to move closer to real repeat arrays in terms of one or more of these metrics, but neither of them approach *A. thaliana* repeat array organisation in terms of all four, making it a useful metric to gauge overall model performance.

![][image7]
Figure X. Composite HOR metric (composite_hor_metric_overlaid.png)
The composite HOR metric density was plotted for each group, representing the multiplicative combined influences of HOR block gap, diversity, similarity and size. The two simulated datasets are similar to each other in terms of the shape of their distribution, while the A. thaliana dataset has orders of magnitude higher scores in this metric, highlighting its potential use as an overall quality metric in future simulations.


A useful feature of this composite metric is to be able to filter real centromeric data to discover large and unexpected long-range repeat organisation, likely due to the occurrence of long range and high magnitude duplications during crossover mutation events.

# Discussion

The results of the HOR similarity and HOR diversity calculations clearly indicate that the differences in block gap distributions and HORs per kb in the Kinetochore-associated Recombination model compared with the Uniform Recombination model are in no small part due to a simple excess of repeat similarity as a result of targeting INDEL events to self-similar regions. It’s clear that the fundamental generative processes underpinning repeat diversity are correctly captured by a simple INDEL/SNP model, at least in terms of HOR metrics, but that the appearance of HORs with high intra-block diversity yet low inter-block diversity, along with a large block gap and block size cannot be explained by this simple model. In simpler terms, the model goes some way to explain accurately the evolution of repeat sequences themselves but not the organisation of these repeats along the centromere.  
The other primary phenomenon that goes unexplained by the simple Uniform Recombination model is the deviation of HOR unit sizes from a power law distribution, forming instead a truncated power law distribution, suggesting there is a missing multiplicative mechanism tying together the size of HOR units and the magnitude of expansion events. Future investigations modelling expansion events as a function of self-similarity may go some way in testing this hypothesis definitively.

Both a random mutation model with a balance of diversifying and homogenising mutations and a directed model with recombination machinery targeting likely kinetochore attachment regions are able to produce centromere-like arrays, but it’s not possible to conclusively support one model over the other in this simulation-based approach. What is supported, however, is the validity of the kinetochore-associated recombination model in comparison to the uniform recombination model, as both produce arrays with centromere-like characteristics. That said, it’s certainly true that the driving mechanism that produces self-similar satellite “patches” in the simulations (and apparently in biological centromere arrays) is balanced, frame-preserving insertion-deletion events combined with single nucleotide polymorphisms.  
Given the high rates of transposon insertion events and their targeting to centromere-specific proteins (Tsukahara paper), it’s clear that their contribution to the evolution of centromeres cannot be ignored. Further simulation-based experiments would benefit from including high frequency transposon insertions. If their addition precludes the emergence of HORs in the uniform recombination or kinetochore-associated recombination models, this would suggest a biological mechanism for purging these transposable elements from the arrays. Unfortunately, testing this hypothesis was not possible within the time constraints of the rotation project due to the added complexity of modelling frame-disrupting insertions in our wholly in-frame mutation paradigm.

Another failing of the current approach is that we ignore the effects of centromere drive, a mechanistic theory of centromere evolution popularised by Henikoff et al in 2001 (henikoff 2001 paper). This model proposes that asymmetric female meiosis might impose a selection pressure on centromeres composed of tandemly-repeated satellite arrays, with repeat arrays having characteristics that allow them to recruit more kinetochore proteins such as greater homogeneity or larger repeat arrays or some combination of both of these characteristics results in their selection over another competing centromere array. It’s conceivable that such a model could explain the surprising homogeneity of centromere arrays in isolation, as well as in addition to the molecular turnover models presented here and in (Dong et al preprint).

Finally, the usability of a monte carlo simulation involving the shuffling and editing of DNA represented as ASCII strings in memory is limited in exploring more creative mechanistic models. The necessary overheads of inserting/deleting elements from large arrays millions of times per simulation impose a substantial calculation burden, limiting the types of mechanisms that can be modelled in reasonable computational timescales. Modelling and testing other moving parts in combination with the current model for example CENH3 positions, methylation, centromere drive, or including the rest of the genome in the evolution (for example it was shown that a centromere allele linked to centromere drive was not fixed in the population due to linked deleterious effects inherited with the allele, fishman and kelly paper 2015\) would require more and more time to be devoted to optimising the simulation. By instead capturing the high-level mechanistic effects of centromere turnover rather than trying to discretely model individual mutation events would resolve this issue \- either by using an existing, heavily optimised evolutionary framework such as SLiM (slim 5 paper) or by turning to numerical dynamical systems as representations of sequences and the relationships between them. Due to the complex nature of representing the mutation and inheritance patterns of centromere arrays as discrete alleles, it’s likely that the latter approach would be better suited to further more complex simulation experiments.

In planning the project and the order of the experiments, in hindsight it would have been better to start with comparing the results of the uniform recombination model with the existing natural centromeres in the 66 assembled Arabidopsis accessions to understand its shortcomings more completely before implementing modifications to the model. On the other hand, the kinetochore-associated recombination model is widely cited by numerous groups studying centromere evolution (gao2025global, naish2024structure), and a computational model of centromeric repeats capturing this had not yet been implemented, making this an interesting angle to consider in its own right.

Overall, the results of this work demonstrate that a model capturing the dynamics and end results of centromere evolution should be based on comparatively frequent INDELs and low frequency SNPs, compared to chromosome arms, but that it is clear that long-range and large-scale duplications are an inherent feature of natural centromere arrays that are identifiable as higher order repeats. The models proposed here and by Dong et al. produce higher order repeats but the repeat blocks are either too self-similar, not far enough apart, too small, or not diverse enough internally. In other words, large tandem duplications occur locally but are quickly broken up by other duplications and mutations. A model either modelling a much higher frequency of tandem duplications over long ranges or a mechanism for reducing the effects of their degradation would go some way in explaining the mutational dynamics of eukaryotic centromere satellite repeats.  
The effect of transposon invasions also cannot be ignored, as it is a strong and potentially driving mutational force in centromere arrays \- the introduction of these events in the current simulation paradigm is not straightforward and may necessitate alternate methods.

# Supplementary Figures

# References

@article{mula2018faster,  
  title={Faster population counts using AVX2 instructions},  
  author={Mu{\\l}a, Wojciech and Kurz, Nathan and Lemire, Daniel},  
  journal={The Computer Journal},  
  volume={61},  
  number={1},  
  pages={111--120},  
  year={2018},  
  publisher={Oxford University Press}  
}

@INPROCEEDINGS{8820962,  
  author={Vasimuddin, Md. and Misra, Sanchit and Li, Heng and Aluru, Srinivas},  
  booktitle={2019 IEEE International Parallel and Distributed Processing Symposium (IPDPS)},   
  title={Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems},   
  year={2019},  
  volume={},  
  number={},  
  pages={314-324},  
  keywords={Kernel;Bioinformatics;Sequential analysis;Genomics;Acceleration;Program processors;Multicore processing;sequence mapping;BWA;acceleration;architecture aware implementation;multicore CPU},  
  doi={10.1109/IPDPS.2019.00041}}

@article{naish2021genetic,  
  title={The genetic and epigenetic landscape of the Arabidopsis centromeres},  
  author={Naish, Matthew and Alonge, Michael and Wlodzimierz, Piotr and Tock, Andrew J and Abramson, Bradley W and Schm{\\"u}cker, Anna and Mand{\\'a}kov{\\'a}, Terezie and Jamge, Bhagyshree and Lambing, Christophe and Kuo, Pallas and others},  
  journal={Science},  
  volume={374},  
  number={6569},  
  pages={eabi7489},  
  year={2021},  
  publisher={American Association for the Advancement of Science}  
}

@article{wlodzimierz2023cycles,  
  title={Cycles of satellite and transposon evolution in Arabidopsis centromeres},  
  author={Wlodzimierz, Piotr and Rabanal, Fernando A and Burns, Robin and Naish, Matthew and Primetis, Elias and Scott, Alison and Mand{\\'a}kov{\\'a}, Terezie and Gorringe, Nicola and Tock, Andrew J and Holland, Daniel and others},  
  journal={Nature},  
  volume={618},  
  number={7965},  
  pages={557--565},  
  year={2023},  
  publisher={Nature Publishing Group UK London}  
}

@article{henderson2025cyclical,  
  title={Cyclical evolution of centromere architecture across 193 eukaryote species},  
  author={Henderson, Ian and W{\\l}odzimierz, Piotr and Perez-Roman, Estela and Hong, Michael and Zhang, Meng and Oliveira, Ludmila and Gonzalez-Isa, Jacob and Jenike, Katharine and Burns, Robin and Zhou, Chenxi and others},  
  year={2025}  
}

@article{tsukahara2025centrophilic,  
  title={Centrophilic retrotransposon integration via CENH3 chromatin in Arabidopsis},  
  author={Tsukahara, Sayuri and Bousios, Alexandros and Perez-Roman, Estela and Yamaguchi, Sota and Leduque, Basile and Nakano, Aimi and Naish, Matthew and Osakabe, Akihisa and Toyoda, Atsushi and Ito, Hidetaka and others},  
  journal={Nature},  
  volume={637},  
  number={8046},  
  pages={744--748},  
  year={2025},  
  publisher={Nature Publishing Group UK London}  
}

@article{henikoff2001centromere,  
  title={The centromere paradox: stable inheritance with rapidly evolving DNA},  
  author={Henikoff, Steven and Ahmad, Kami and Malik, Harmit S},  
  journal={Science},  
  volume={293},  
  number={5532},  
  pages={1098--1102},  
  year={2001},  
  publisher={American Association for the Advancement of Science}  
}

@article{fishman2015centromere,  
  title={Centromere-associated meiotic drive and female fitness variation in Mimulus},  
  author={Fishman, Lila and Kelly, John K},  
  journal={Evolution},  
  volume={69},  
  number={5},  
  pages={1208--1218},  
  year={2015},  
  publisher={Blackwell Publishing Inc Malden, USA}  
}

@article{haller2025slim,  
  title={SLiM 5: Eco-evolutionary simulations across multiple chromosomes and full genomes},  
  author={Haller, Benjamin C and Ralph, Peter L and Messer, Philipp W},  
  journal={Molecular Biology and Evolution},  
  pages={msaf313},  
  year={2025},  
  publisher={Oxford University Press}  
}

@article{gao2025global,  
  title={A global view of human centromere variation and evolution},  
  author={Gao, Shenghan and Oshima, Keisuke K and Chuang, Shu-Cheng and Loftus, Mark and Montanari, Annalaura and Gordon, David S and Human Genome Structural Variation Consortium and Human Pangenome Reference Consortium and Hsieh, PingHsun and Konkel, Miriam K and others},  
  journal={bioRxiv},  
  pages={2025--12},  
  year={2025},  
  publisher={Cold Spring Harbor Laboratory}  
}

@article{wlodzimierz2023trash,  
  title={TRASH: tandem repeat annotation and structural hierarchy},  
  author={Wlodzimierz, Piotr and Hong, Michael and Henderson, Ian R},  
  journal={Bioinformatics},  
  volume={39},  
  number={5},  
  pages={btad308},  
  year={2023},  
  publisher={Oxford University Press}  
}

@software{tange\_2025\_17692695,  
      author       \= {Tange, Ole},  
      title        \= {GNU Parallel 20251122 ('Mamdani')},  
      month        \= Nov,  
      year         \= 2025,  
      note         \= {{GNU Parallel is a general parallelizer to run  
                       multiple serial command line programs in parallel  
                       without changing them.}},  
      publisher    \= {Zenodo},  
      doi          \= {10.5281/zenodo.17692695},  
      url          \= {https://doi.org/10.5281/zenodo.17692695}  
}  
@Article{         harris2020array,  
 title         \= {Array programming with {NumPy}},  
 author        \= {Charles R. Harris and K. Jarrod Millman and St{\\'{e}}fan J.  
                 van der Walt and Ralf Gommers and Pauli Virtanen and David  
                 Cournapeau and Eric Wieser and Julian Taylor and Sebastian  
                 Berg and Nathaniel J. Smith and Robert Kern and Matti Picus  
                 and Stephan Hoyer and Marten H. van Kerkwijk and Matthew  
                 Brett and Allan Haldane and Jaime Fern{\\'{a}}ndez del  
                 R{\\'{i}}o and Mark Wiebe and Pearu Peterson and Pierre  
                 G{\\'{e}}rard-Marchant and Kevin Sheppard and Tyler Reddy and  
                 Warren Weckesser and Hameer Abbasi and Christoph Gohlke and  
                 Travis E. Oliphant},  
 year          \= {2020},  
 month         \= sep,  
 journal       \= {Nature},  
 volume        \= {585},  
 number        \= {7825},  
 pages         \= {357--362},  
 doi           \= {10.1038/s41586-020-2649-2},  
 publisher     \= {Springer Science and Business Media {LLC}},  
 url           \= {https://doi.org/10.1038/s41586-020-2649-2}  
}

@article{liebovitch1989fast,  
  title={A fast algorithm to determine fractal dimensions by box counting},  
  author={Liebovitch, Larry S and Toth, Tibor},  
  journal={physics Letters A},  
  volume={141},  
  number={8-9},  
  pages={386--390},  
  year={1989},  
  publisher={Elsevier}  
}

@article{grassberger1983characterization,  
  title={Characterization of strange attractors},  
  author={Grassberger, Peter and Procaccia, Itamar},  
  journal={Physical review letters},  
  volume={50},  
  number={5},  
  pages={346},  
  year={1983},  
  publisher={APS}  
}

@article{naish2024structure,  
  title={The structure, function, and evolution of plant centromeres},  
  author={Naish, Matthew and Henderson, Ian R},  
  journal={Genome Research},  
  volume={34},  
  number={2},  
  pages={161--178},  
  year={2024},  
  publisher={Cold Spring Harbor Lab}  
}  
@article{sweeten2024moddotplot,  
  title={ModDotPlot—rapid and interactive visualization of tandem repeats},  
  author={Sweeten, Alexander P and Schatz, Michael C and Phillippy, Adam M},  
  journal={Bioinformatics},  
  volume={40},  
  number={8},  
  pages={btae493},  
  year={2024},  
  publisher={Oxford University Press}  
}
