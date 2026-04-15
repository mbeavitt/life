# Real-Time Genomics: The Custom Hardware Revolution
### Hardware IP Licensing for the Genomics Industry

---

## Executive Summary

The genomics industry is generating data faster than it can analyse it. Next-generation sequencing throughput has scaled by orders of magnitude over the past decade, yet the computational infrastructure used to process that data has barely changed: slow, general-purpose CPUs running software pipelines on cloud infrastructure that is increasingly unreliable, expensive, and geopolitically fraught.

This proposal outlines a plan to solve that problem not by building another end-to-end analysis platform, but by licensing the fundamental computational building blocks – hardware-accelerated algorithmic cores for sequence alignment, k-mer indexing, variant calling, and methylation calling – directly to the companies that build sequencing instruments, cloud platforms, and clinical diagnostic systems. The model is deliberately analogous to ARM Holdings: we do not sell phones, but the goal is for this architecture to be inside every sequencing instrument that matters.

The secondary genomics analysis market is valued at approximately $1.1 billion (2025) and growing at 15% annually. The core computational bottlenecks are well understood, performance-critical, and structurally under-served by existing solutions. We are seeking £500,000 in seed funding to bring the first generation of licensed IP cores to commercial pilot stage within 18 months.

---

## 1. The Problem: Genomics Has a Compute Crisis

### 1.1 The Bottleneck Is Not the Sequencer

The cost of sequencing a human genome has fallen from ~$100 million in 2001 to under $200 today, a decline that outpaces Moore's Law by a significant margin. The cost of *analysing* that genome, however, has not followed the same curve. A standard whole-genome sequencing (WGS) variant calling pipeline typically takes 20 or more hours on general-purpose hardware. A single genome, from raw reads to clinical variant report, may cost more in compute than the sequencing itself – and as sequencing costs continue to fall, compute will become the dominant line item.

The algorithms are mature. The bottleneck is that the industry is running fundamentally parallel, data-intensive workloads – string matching, graph traversal, probabilistic inference – on hardware designed for serial, general-purpose computation. The mismatch is structural; no amount of software optimisation closes that gap.

### 1.2 Cloud Dependence Is Becoming a Liability

The genomics industry has responded to its compute problem largely by offloading to cloud infrastructure. This is rational in the short term but increasingly untenable as a long-term strategy. Critical cloud service disruptions from AWS, Azure, and GCP have increased by 52% since 2022, with an 18% annual increase in 2024 alone. Total critical outage duration reached 221 hours across the major providers in 2024, up 51% since 2022 [1]. Between August 2024 and August 2025, the three major providers together experienced more than 100 service outages [2]. The October 2025 AWS outage alone generated over 17 million Downdetector reports, lasted over 15 hours, and disrupted services across more than 60 countries [3].

For clinical genomics – where a delayed variant call can mean a delayed diagnosis – this is not an acceptable reliability profile. Regulatory pressure in the EU and UK is already pushing clinical laboratories toward on-premises, auditable compute for sensitive genomic data. Geopolitical tensions are prompting research institutions to reconsider dependence on US-headquartered hyperscalers for critical scientific infrastructure. Prices are rising as the AI compute boom strains cloud capacity, and spot instance reclaims are increasingly disrupting long-running genomics jobs.

The industry is reaching a structural inflection point. Genomics compute is moving closer to the hardware. The question is who owns the IP when it does.

> [1] Parametrix Cloud Outage Risk Report 2024, via CXO Today: https://cxotoday.com/news-analysis/have-cloud-outages-become-more-frequent-or-is-the-internet-more-dependent/
> [2] Rest of World, "Cloud Outages 2025: Global Business Impact": https://restofworld.org/2026/cloud-outages-2025-global-business-impact/
> [3] TechTarget, "Cloud Outages Expected to Be the New Normal in 2026": https://www.techtarget.com/searchcloudcomputing/feature/Cloud-outages-expected-to-be-the-new-normal-in-2026

---

## 2. Proof of Concept: Hardware Acceleration Works

Three precedents from adjacent industries demonstrate that custom hardware is the proven solution to this class of problem.

### 2.1 Apple Silicon: The Unified Memory Revolution

When Apple introduced the M1 chip in November 2020, it embedded dedicated fixed-function accelerators – for video transcoding, neural inference, image signal processing, and security – directly alongside CPU and GPU cores in a unified memory architecture [4]. The result was a step-change in performance-per-watt that redefined the laptop category. Committing specific silicon area to specific workloads, even at the cost of flexibility, yields returns that software optimisation alone cannot approach.

Genomics analysis presents exactly the same profile: well-defined, recurring, high-throughput workloads that run millions of times per day across thousands of laboratories worldwide.

> [4] Apple Inc., "Apple unleashes M1," Apple Newsroom, November 2020

### 2.2 Illumina DRAGEN: The $100 Million Proof Point

Illumina, the dominant sequencing platform company, recognised the compute bottleneck early and in 2018 acquired Edico Genome for approximately $100 million [5]. The resulting DRAGEN Bio-IT Platform uses a dedicated FPGA-based processing unit to complete a WGS variant calling analysis in under 30 minutes – a reduction from the 20+ hours required on conventional hardware. DRAGEN is now embedded in Illumina's flagship sequencing instruments and has become the industry's de facto performance benchmark.

This validates both the technical approach and the commercial model: a major industry player was willing to pay nine figures for hardware-accelerated genomics IP. Our proposition is to provide that capability as a licensable layer, not a proprietary product locked to a single vendor.

> [5] Illumina Inc., "Illumina Acquires Edico Genome," Press Release, 2018

### 2.3 CERN LHC: When Software Simply Cannot Keep Up

CERN's Large Hadron Collider generates up to one billion proton-proton collisions per second. The trigger system must decide, in under one microsecond, whether to retain or discard the resulting collision data – a latency requirement that conventional computing architectures cannot meet [6]. CERN's solution was to train a neural network to recognise scientifically significant collision events and compile it directly to FPGA fabric, bypassing the CPU entirely. The silicon *is* the algorithm.

This represents the logical endpoint of hardware-software co-design: when the data rate is high enough and the latency requirement strict enough, the algorithm must become silicon. Genomic sequencing instruments are approaching this regime.

> [6] Duarte, J. et al., "Fast inference of deep neural networks in FPGAs for particle physics," Journal of Instrumentation 13 (2018) P07027. https://arxiv.org/abs/1804.06913

---

## 3. Market Analysis

### 3.1 Addressable Market

The NGS secondary analysis market – covering the computational processing of raw sequencing reads into variant calls, methylation calls, and structural annotations – was valued at approximately $1.1 billion in 2025 and is growing at approximately 15% CAGR [7], driven by the continued decline of sequencing costs, the expansion of clinical genomics programmes, and the growth of large-scale population genomics initiatives. At this growth rate, the market reaches approximately $2.2 billion by 2030.

This is the initial target. The broader genomics informatics market, which includes tertiary analysis, clinical decision support, and data management, is a natural adjacency and is valued at over $5 billion. The FPGA and custom silicon in life sciences market is nascent but growing rapidly, accelerated by regulatory pressure for validated compute and the AI-driven demand for inference-optimised hardware.

### 3.2 Customer Segments

**Tier 1 – Sequencing Platform OEMs (Primary Target)**
Illumina, Oxford Nanopore Technologies, PacBio, and emerging players such as Element Biosciences and Singular Genomics. These companies have both the engineering capability to integrate hardware IP and the commercial motivation to differentiate on analysis speed and on-device processing. The Illumina/Edico precedent demonstrates willingness to pay.

**Tier 2 – Clinical Cloud Platforms**
Companies such as DNAnexus, Sentieon, and Seqera Labs provide managed genomics analysis infrastructure to hospital systems and clinical labs. As regulatory pressure mounts for auditable, deterministic compute, hardware-validated IP cores become a compliance asset.

**Tier 3 – Clinical Diagnostic Laboratories (Indirect)**
Hospital genetics departments, cancer genomics centres, and population health programmes that run on-premises compute. These customers are served indirectly through licensing to OEMs and platform providers, but represent the end-demand signal that makes the licensing market real.

### 3.3 Market Dynamics

The genomics compute stack is fragmenting in ways that create demand for a neutral, hardware-aware algorithmic layer:

- Sequencing platforms are diverging. Illumina short-read, Oxford Nanopore long-read, and PacBio HiFi each demand distinct algorithmic approaches (seed-extension alignment, signal-level basecalling, CCS consensus generation). No single software stack serves all three well.
- Deployment environments are splitting. On-premises clinical systems, cloud-based population genomics, and edge-deployed point-of-care diagnostics each carry different latency, cost, and regulatory constraints.
- Regulatory pressure is increasing. The EU's In Vitro Diagnostic Regulation (IVDR) and equivalent UK frameworks are pushing clinical labs toward validated, auditable compute components rather than opaque end-to-end software platforms.

A clinical lab analysing methylation data on-premises and short-read Illumina data in the cloud should not require two entirely different software stacks to perform what is, at its core, the same linear-time string matching and probabilistic variant inference.

> [7] MarketsandMarkets, "Next Generation Sequencing (NGS) Data Analysis Market – Global Forecast to 2029," 2024

---

## 4. The Solution

### 4.1 Hardware-Accelerated Genomics Primitives

We propose to develop and license a library of hardware-accelerated computational primitives for genomics: modular, silicon-ready algorithmic cores that can be embedded by OEMs and platform providers into their own products.

The initial core library targets the five highest-value bottlenecks in the secondary analysis pipeline:

| Core | Function | Performance Gain (vs. CPU) |
|------|----------|---------------------------|
| SeqAlign-HW | Seed-extension sequence alignment (BWA-MEM2 equivalent) | ~20–40× |
| KmerIndex-HW | K-mer indexing and lookup for de Bruijn graph assembly | ~15–30× |
| VariantCall-HW | Haplotype-based variant calling (HaplotypeCaller equivalent) | ~10–25× |
| MethylCall-HW | CpG methylation calling from Nanopore signal data | ~20–50× |
| DupMark-HW | Duplicate marking and read sorting (Picard-equivalent) | ~8–15× |

Each core is delivered as:
- A portable HDL implementation (Verilog/VHDL) for FPGA integration
- An annotated specification for ASIC tape-out integration
- A software-layer API enabling seamless integration into existing bioinformatics pipelines (GATK-compatible, CWL/Nextflow-compatible)

### 4.2 The ARM Analogy

Our strategic positioning is deliberately analogous to ARM Holdings. ARM does not manufacture chips. It designs instruction set architectures and processor microarchitectures and licenses them to companies that do. Every major smartphone processor in the world – Apple A-series, Qualcomm Snapdragon, Samsung Exynos – is built on ARM IP. ARM earns a royalty on every chip shipped, without bearing the capital cost of fabrication.

The aim is to occupy the same position in the genomics compute stack: providing the algorithmic substrate that makes instruments faster and analysis cheaper, without manufacturing instruments, operating platforms, or delivering reports. The goal is to become the default compute layer for genomics by making every end product better.

---

## 5. Innovation & Competitive Advantage

### 5.1 What Exists Today and Why It Falls Short

| Competitor | Approach | Limitation |
|---|---|---|
| Illumina DRAGEN | FPGA acceleration for Illumina-specific pipelines | Proprietary; locked to Illumina instruments; no third-party licensing |
| NVIDIA Parabricks | GPU-accelerated GATK on NVIDIA hardware | GPU only; no FPGA/ASIC path; cloud-dependent; general-purpose hardware |
| Sentieon | Optimised CPU software reimplementations | Software only; no hardware advantage; does not address latency floor |
| AWS Graviton (LifeOmics etc.) | Cloud-optimised CPU instances | Cloud-dependent; general-purpose; no custom silicon |

None of the existing players operates at the IP licensing layer. DRAGEN is an internal Illumina product – it is not licensed to Oxford Nanopore or PacBio. Parabricks is a software product that requires NVIDIA GPUs. There is no ARM-equivalent for genomics compute: no neutral, hardware-aware, platform-agnostic primitives library that any OEM can embed.

This gap is structural, not accidental. It exists because building at this layer requires rare expertise at the intersection of bioinformatics algorithm development, digital hardware design, and semiconductor IP commercialisation. That intersection is the basis of this proposal.

### 5.2 Technical Innovation

The cores are designed from the ground up for hardware execution. The primary innovation is hardware-algorithm co-design: restructuring the computational graphs of genomics algorithms to exploit the parallelism intrinsic to FPGA fabric, rather than mapping software control flow onto hardware. Key technical differentiators include:

- Pipelined dataflow architecture for alignment kernels that sustains near-100% hardware utilisation at sustained read throughput
- Approximate nearest-neighbour seeding in hardware using configurable Bloom filters – eliminates the primary memory bottleneck in standard seed-extension alignment
- Compressed de Bruijn graph traversal implemented as a streaming hardware state machine, avoiding the graph materialisation step that dominates RAM usage in software assemblers
- Signal-level methylation inference as a recurrent hardware pipeline compatible with Oxford Nanopore's raw squiggle format, eliminating the separate basecalling step

These are novel implementations – distinct from prior art – and form the basis of an initial patent portfolio.

---

## 6. Business Model & Route to Market

### 6.1 Revenue Streams

**Stream 1: Per-Unit IP Licensing (OEM Integration)**
Sequencing instrument OEMs license individual IP cores for integration into their instruments' FPGA or ASIC fabric. Pricing is on a per-unit-shipped royalty basis, typically in the range of $1–5 per instrument unit depending on core complexity and exclusivity terms. At an OEM scale of 10,000–50,000 instruments shipped per year, this generates meaningful recurring revenue without marginal cost.

**Stream 2: Annual Platform License (Cloud Providers)**
Cloud genomics platforms and managed analysis providers license firmware images and API-compatible software bindings for deployment in virtualised FPGA instances (AWS F1, Azure NP-series, Google Cloud FPGA). Priced at £50,000–£200,000 per year per platform, depending on throughput and support tiers.

**Stream 3: Development & Integration Fees (Early Stage)**
In the pilot phase, OEM and platform partners pay a one-time integration and validation fee (£75,000–£200,000) covering HDL delivery, testbench packages, regulatory documentation support, and engineering time during integration. This provides early revenue while the royalty stream matures.

**Stream 4: ASIC Licensing (Long-Term)**
As adoption scales, customers moving to custom ASIC designs (integrating our cores into purpose-built chips) pay a per-wafer or per-die royalty analogous to ARM's model. This is a Year 3+ opportunity but will be the highest-margin, most defensible revenue stream.

### 6.2 Route to Market

**Phase 1 (Months 0–6): Anchor Partnership**
The route to market begins with a single, high-credibility anchor partnership – an OEM or clinical institution willing to co-develop and validate the first core (SeqAlign-HW) in exchange for preferential licensing terms. Three targets have been identified based on strategic fit, geographic proximity, and alignment with their existing computational challenges:

- **Azenta Life Sciences** – Azenta is an established genomics sequencing services provider processing large volumes of samples annually through standardised Illumina-based pipelines for academic and pharmaceutical research customers. Their competitive position is built on turnaround time and cost: faster analysis is margin. Unlike sequencing platform manufacturers, Azenta has no internal hardware stack and no incentive to build one – their core business is services, not silicon. This makes them an ideal early licensing partner: they provide high-volume, real-world production data for benchmarking, a direct commercial motivation to integrate faster compute primitives, and a credible case study for the broader sequencing services market. A validated deployment at Azenta demonstrates that our cores work at scale on production Illumina data – the most common sequencing substrate in the industry.

- **ARM Holdings (Cambridge)** – ARM's Cambridge headquarters makes them an obvious strategic partner, and the analogy is more than rhetorical. ARM has a direct interest in seeing the IP licensing model extended into high-value scientific domains, and their Ecosystem Acceleration Program has a history of supporting hardware IP startups. A strategic relationship with ARM – whether as investor, advisor, or licensing partner for ASIC integration – would materially validate our positioning and open their extensive OEM network.

- **Genomics England / NHS Genomics** – As the operator of one of the world's largest clinical whole-genome sequencing programmes, Genomics England is the highest-throughput clinical customer for genomics compute in the UK. Their existing relationships with sequencing OEMs and cloud providers, combined with their mandate for auditable, validated compute in clinical pipelines, make them an ideal early validation partner for the regulatory credibility our IP cores provide.

The Wellcome Sanger Institute (Hinxton) is a further natural early relationship, providing access to real, large-scale sequencing datasets for benchmarking and the scientific credibility of a co-authored application note.

**Phase 2 (Months 6–18): Pilot Deployment & Validation**
Deliver validated HDL cores to the anchor partner. Collect performance benchmarks against DRAGEN and CPU baselines. Generate the co-authored technical paper and/or application note that establishes credibility with the broader OEM community. Begin conversations with two additional OEM partners.

**Phase 3 (Months 18–36): Commercialisation**
Convert pilot agreements to commercial licensing arrangements. Expand core library. Hire dedicated business development resource. Present at ASHG, AGBT, and ISMB to establish technical credibility with the genomics engineering community.

**Phase 4 (Months 36+): Platform & ASIC**
Expand to cloud platform licensing. Begin ASIC licensing conversations with customers reaching instrument volumes where custom silicon becomes economically attractive.

---

## 7. Defensibility

### 7.1 Intellectual Property

The primary defensive moat is a structured patent portfolio covering novel FPGA and ASIC implementations of core genomics algorithms. Patent applications will be filed in the UK, US, and EU for:

- The pipelined dataflow architecture for seed-extension alignment
- Bloom-filter-based approximate k-mer seeding on FPGA fabric
- The streaming compressed de Bruijn graph traversal state machine
- Signal-level methylation inference as a recurrent hardware pipeline

Filing occurs in parallel with core development. Provisional applications are prioritised to establish priority dates before any publication or presentation.

### 7.2 Regulatory Moat

Clinical genomics is a regulated domain. Once a hardware IP core has been validated and documented for use in an IVD-compliant pipeline – a process that is expensive, time-consuming, and requires maintained quality management systems – the switching cost for a clinical OEM customer is very high. Being the first validated, hardware-accelerated primitives provider creates a durable competitive position in clinical markets that goes well beyond patent protection.

### 7.3 Deep OEM Integration

Hardware IP embedded in instrument firmware is not easily replaced. Unlike software libraries that can be swapped, HDL cores are integrated into FPGA bitstreams and ASIC designs at tape-out. Once a sequencing platform OEM has designed our cores into a product generation, replacement requires a new hardware revision cycle, typically 18–24 months. This creates strong retention without lock-in contracts.

### 7.4 Talent & Know-How

The combination of skills required – bioinformatics algorithm expertise, digital hardware design, FPGA/ASIC IP commercialisation – is rare, and its scarcity is part of why the gap exists. As the team grows, the accumulation of internal expertise and proprietary optimisations beyond what is patented constitutes a durable knowledge moat.

---

## 8. Team

### 8.1 The Founder

The founding thesis emerges from direct, lived experience at the computational bottleneck. The founder is completing a PhD in bioinformatics at the University of Cambridge, with research focused on the simulation of centromere evolutionary dynamics – work that requires large-scale sequence analysis at the boundary of what current tools can handle. Alongside the PhD, they bring three years of professional experience designing and maintaining Nextflow pipelines for large-scale genomics analysis – running the jobs, watching the queues, and directly experiencing where the compute breaks down.

The founder sits at an unusual intersection: deep bioinformatics expertise combined with a strong personal grounding in computational hardware, electronics, and systems programming. It is the specific combination of skills that allows someone to look at a genomics pipeline and see misaligned hardware rather than slow code, and to have the technical vocabulary to do something about it.

The problem was identified from inside the pipeline.

### 8.2 Building the Team

There is currently no co-founder or formal team – this is a pre-seed concept at the stage of seeking initial investment and partnerships. The roles needed are specific enough that filling them with the wrong people early would be worse than taking the time to find the right ones.

The critical hires are: a technical co-founder with hardware IP or FPGA engineering experience, a business development lead with OEM licensing experience, and one or two scientific advisors with clinical genomics or sequencing industry connections who can open doors to early pilot partners. Several individuals with relevant profiles have expressed interest in the idea. The goal is to use the credibility of a funded, validated venture to convert those conversations into commitments – which is part of what this pitch is for.

### 8.3 Key Roles to Fill (Post-Funding)

| Role | Priority | Purpose |
|---|---|---|
| HDL/FPGA Engineer (or co-founder) | Immediate | Core development – cannot proceed without this |
| Bioinformatics Algorithm Engineer | Early | Algorithm-hardware co-design |
| Scientific Advisor (×1–2) | Early | Access to real datasets; clinical and industry connections |
| Business Development | Month 6–12 | OEM outreach and licensing negotiations |

---

## 9. Work Plan

### Phase 0 – Foundation (Pre-Funding, Current Stage)
- Conduct IP landscape analysis (freedom-to-operate searches for core algorithms on FPGA)
- Identify and approach target anchor partners (Azenta, Sanger, Genomics England)
- Begin conversations with potential technical co-founders and scientific advisors
- File provisional patent applications on primary algorithmic innovations

### Phase 1 – Core Development (Months 1–6)
- Bring on HDL/FPGA engineer (co-founder or first hire)
- Develop and validate SeqAlign-HW (primary alignment core) on FPGA dev hardware
- Develop and validate KmerIndex-HW (k-mer indexing core)
- Deliver prototype cores to anchor partner for integration testing
- File first full patent applications

### Phase 2 – Pilot Validation (Months 6–18)
- Partner integration of SeqAlign-HW and KmerIndex-HW complete
- Benchmarking against DRAGEN and CPU baselines on partner sequencing data
- Develop VariantCall-HW (haplotype-based variant calling core)
- Submit co-authored technical paper / application note
- Begin second and third OEM partnership conversations
- File patents on VariantCall-HW and MethylCall-HW

### Phase 3 – Commercialisation (Months 18–36)
- First commercial licensing agreements signed
- Expand core library to full five-core suite
- Hire Business Development Lead
- Present at AGBT and/or ASHG annual conference
- Begin Series A preparation

---

## 10. Financial Projections

### 10.1 Key Assumptions

- Seed funding of £500,000 secured at close
- First pilot agreement signed within 6 months (£100K integration fee)
- First commercial licensing deal by Month 18
- OEM royalty: $2 per unit, partner shipping 5,000 units/year (Year 2), scaling to 20,000/year (Year 4)
- Cloud platform license: £100K/year per customer; first customer Year 2, three customers by Year 4
- Staff costs based on UK market rates; team grows from 2 to 6 over 36 months
- FPGA tooling (Xilinx/Intel Quartus) and EDA licensing: £40K/year
- Patent prosecution (UK, US, EU): £15K per application, four applications over 24 months

### 10.2 Cost Projections

| Cost Category | Year 1 | Year 2 | Year 3 |
|---|---|---|---|
| Personnel (salaries) | £220,000 | £420,000 | £620,000 |
| Hardware & Tooling | £45,000 | £50,000 | £55,000 |
| IP & Legal | £60,000 | £40,000 | £30,000 |
| Business Development & Travel | £15,000 | £35,000 | £50,000 |
| Overhead & Operations | £20,000 | £30,000 | £40,000 |
| **Total Costs** | **£360,000** | **£575,000** | **£795,000** |

### 10.3 Revenue Projections

| Revenue Stream | Year 1 | Year 2 | Year 3 | Year 4 | Year 5 |
|---|---|---|---|---|---|
| Integration & Dev Fees | £100,000 | £150,000 | £100,000 | £100,000 | £100,000 |
| OEM Royalties (per-unit) | £0 | £50,000 | £150,000 | £550,000 | £1,200,000 |
| Cloud Platform Licenses | £0 | £100,000 | £200,000 | £300,000 | £500,000 |
| ASIC Licensing | £0 | £0 | £0 | £150,000 | £400,000 |
| **Total Revenue** | **£100,000** | **£300,000** | **£450,000** | **£1,100,000** | **£2,200,000** |
| **Net Position** | **–£260,000** | **–£275,000** | **–£345,000** | **+£305,000** | **+£1,405,000** |

*Cumulative cash requirement through to profitability: ~£880,000. A seed round of £500,000 is projected to fund operations through to first commercial licensing at around Month 18, at which point a Series A would be sought to fund full commercialisation. These are illustrative projections – the key variables are time-to-anchor-partnership and OEM adoption speed.*

---

## 11. Risk Analysis

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| **Technical: cores fail to meet benchmark targets** | Medium | High | Conservative targets based on published DRAGEN benchmarks; internal validation before partner delivery; phased delivery with pilot gates |
| **Commercial: OEM partners slow to adopt** | Medium | High | Multiple parallel partner conversations; integration fee model reduces perceived risk for partners; anchor partnership targeted at most innovation-hungry OEM |
| **Competitive: Illumina or NVIDIA enters IP licensing** | Low | High | Patent portfolio establishes prior art; deep OEM integration creates switching costs; regulatory validation moat takes years to replicate |
| **IP: freedom-to-operate issues** | Medium | Medium | FTO searches completed before development begins; design-around paths identified in advance; provisional filings establish early priority dates |
| **Talent: inability to hire key technical roles** | Medium | High | Equity reserved for technical co-founder; Cambridge university network as a recruitment pipeline; skills gap between bioinformatics and hardware is real and acknowledged |
| **Market: sequencing market consolidation reduces OEM diversity** | Low | Medium | Multiple OEM tiers targeted; cloud platform licensing provides alternative route; clinical lab segment grows regardless of OEM consolidation |
| **Regulatory: IVDR/MDR creates unforeseen compliance burden** | Low | Medium | ISO 13485 QMS established early; regulatory consultant engaged for clinical customer pathway; non-clinical customers de-risk this |
| **Funding: Series A unavailable at 18 months** | Medium | Medium | Revenue from pilot and licensing by Month 18 extends runway; cost structure kept lean; bridge options explored proactively |

---

## 12. The Ask

We are seeking £500,000 in seed investment to fund the first 18 months of operations, covering:

- Initial engineering team (2 FTEs)
- FPGA tooling, hardware development platforms, and EDA licensing
- Patent prosecution across UK, US, and EU jurisdictions
- Anchor partnership development and technical delivery
- Legal, compliance, and QMS establishment

This investment funds us to the point of a validated commercial pilot with a named OEM or platform partner, two or more active licensing negotiations, a defensible patent portfolio, and benchmarked performance data showing material improvement over existing solutions. That is the package required to raise a Series A on strong terms.

Hardware acceleration works – Apple, Illumina, and CERN have each demonstrated it independently in adjacent domains. The gap in the genomics compute stack is well-documented throughout this proposal. What is needed now is the capital and early partnerships to turn a well-grounded concept into a working prototype and a first commercial pilot.

Beyond capital, the most valuable thing an investor or partner could bring at this stage is:
- Introductions to sequencing services companies or OEM engineering teams willing to co-develop a pilot
- Experience in IP licensing or semiconductor go-to-market
- A technical network that helps identify the right hardware engineering co-founder

This is an early-stage proposal from a single founder with deep domain expertise and a clear thesis. The next step is finding the right people to build it with.

---

*Document prepared March 2026. All financial projections are forward-looking estimates based on publicly available market data and disclosed competitor benchmarks. They do not constitute a guarantee of performance.*
