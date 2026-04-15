Title: Real time genomics - the custom hardware revolution 

Aim: To sell the reader on the concept of real time genomics through on-device hardware accelerated genomics analysis via case studies and benchmarks

Case Studies:
    - Apple released a custom M-series chip with unified memory architecture and fixed-function accelerators, significantly disrupting the personal computer market. The tradeoff is that you can't upgrade your memory, and they are betting that a significant proportion of their userbase will make good use of their fixed-function accelerators (video transcoding, AI inference, security, image processing). (cite this plz)
    - Illumina predicted that hardware would pose the biggest bottleneck in variant calling pipelines, and bought Edico Genome for $100 million. Their DRAGEN pipeline now completes a typical whole genome sequencing variant calling analysis in under 30 minutes, compared to the 20+ hours it usually takes. (cite this plz)
    - CERN's LHC generates up to 1 billion proton-proton collisions per second, and the trigger system needs to decide whether or not to keep the resulting collision data in under a microsecond to avoid explosions of data and block downstream analysis. They designed an ingenious custom approach where a neural network trained to recognise significant collision events is compiled directly to FPGA hardware - traditional computer chips are too slow. (cite this plz)

Why not genomics:
Genomics is lagging behind - the vast majority of our computational work
happens on slow general purpose computer hardware, and the industry is in
desperate need for a hardware revolution. The NGS data analysis market is worth
$1.1B as of 2025 and big players like Illumina and cloud compute platforms such
as AWS and GCP are gobbling up more than their fair share. As sequencing costs
come down, the portion of cost allocated to data analysis and processing will
become a primary concern - doubly so, now that the explosion of AI compute is
raising prices.

Resilience/Dependence on external players:
The past 12 months has seen more disruption to core compute services (AWS,
Cloudflare) than any other time period since the start of the information age.
Critical cloud service disruptions from the top three providers (AWS, Azure,
and GCP) have increased by 52% since 2022, with an 18% annual increase in 2024
alone. Total critical outage duration rose to 221 hours in 2024, up 51% since
2022 [1]. Between August 2024 and August 2025, the three major providers
together experienced more than 100 service outages [2]. 2025 was particularly
severe: the October 2025 AWS outage alone generated over 17 million
Downdetector reports and lasted over 15 hours, disrupting services across more
than 60 countries [3].  

Many are beginning to question the wisdom of dependence
on the big players for their core compute infrastructure. Prices are subject to
increase, blackouts and lack of availability of instances leading to spot
reclaims are more frequent than ever, and geopolitical and regulatory pressures
are driving organisations to reconsider dependence on a handful of
US-headquartered hyperscalers for critical infrastructure.

[1] Parametrix Cloud Outage Risk Report 2024, via CXO Today: https://cxotoday.com/news-analysis/have-cloud-outages-become-more-frequent-or-is-the-internet-more-dependent/
[2] Rest of World: https://restofworld.org/2026/cloud-outages-2025-global-business-impact/
[3] TechTarget: https://www.techtarget.com/searchcloudcomputing/feature/Cloud-outages-expected-to-be-the-new-normal-in-2026

The Solution:
What the genomics industry needs is a roadmap to resilient, fast compute
infrastructure. Buying custom hardware and developing it in house is an insane
proposal, with R&D costs dwarfing any analysis or sequencing costs. Relying on
sequencing providers like Illumina to foot the cost and develop fast and
scalable compute infrastructure is a gamble - not only will their products be
necessarily constrained by the requirements of their biggest customers (no
flexibility), their solutions will be proprietary and subject to vendor
lock-in.  

We propose to develop hardware-accelerated computational primitives for genomics: licensed
algorithm IP and FPGA-optimised firmware cores that sequencing platform
companies, cloud providers, and diagnostic labs embed into their own products.
Rather than competing with existing pipeline providers like Illumina's DRAGEN
or NVIDIA's Parabricks at the application layer, we operate one level down,
providing the fundamental building blocks (sequence alignment, k-mer indexing,
variant calling kernels) as licensable, silicon-ready modules. The model is
analogous to ARM's relationship to the smartphone industry: ARM does not sell
phones, but its architecture is inside nearly every phone sold.

The genomics compute stack is fragmenting. Sequencing platforms are
diversifying (Illumina, PacBio, Oxford Nanopore each demand different
algorithmic optimisations), deployment environments are splitting between
on-prem, cloud, and edge, and regulatory pressure is pushing clinical labs
toward validated, auditable compute components rather than opaque end-to-end
platforms. This creates demand for a neutral, hardware-aware algorithmic layer
that can be integrated across platforms without tying customers to a single
vendor's ecosystem. A lab analysing methylation calling data on-prem and Illumina
data in the cloud should not need two entirely separate software stacks to do
fundamentally the same linear-time string matching and probabilistic variant
calling underneath.

Our revenue model combines per-unit IP licensing fees (for FPGA and ASIC
integrations), royalty-based arrangements with platform vendors who ship our
cores inside their instruments, and annual license fees for cloud-deployable
firmware images. The initial target is the secondary analysis market, currently
valued at approximately $1.1 billion and growing at around 15% annually, where
the core computational bottlenecks (alignment, sorting, duplicate marking,
haplotype assembly) are well-defined, performance-critical, and ripe for
hardware acceleration. Early partnerships with sequencing hardware OEMs and
clinical cloud platforms establish us as the default compute substrate for
genomics, the same way ARM became the default instruction set for mobile - not
by building the end product, but by making every end product better.
