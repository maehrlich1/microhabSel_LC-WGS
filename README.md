# Selection among microhabitats using low-coverage WGS

This repository contains the scripts used in the analysis of within-generation selection among *Fundulus heteroclitus* individuals inhabiting distinct microhabitats. Low-coverage whole genome sequencing (LC-WGS) was performed on ~1K individuals sampled in distinct microhabitatss and at different time points to estimate allele frequency changes within a single generation as a result of microhabitat selection. Allele frequency changes were also associated with previously documented phenotypic divergence among microhabitat residents. Phenotypic divergence among microhabitat residents was previously published here (LINK) whereas the results of this LC-WGS project are published in XXXX.

The repository contains markdown files in `/markdowns` detailing the analysis process and parameter justifications. `/hps_scripts` contain code with minimal annotation used to run the analysis on HPC servers at the University of Miami.

### Raw Data

Libraries were prepared from DNA isolated from 963 *Fundulus heteroclitus* individuals. Each individual was uniquely indexed and the library pool was sequenced in a two-stage manner in order to normalize sequencing depth across individuals. A total of 19 lanes were sequenced on the Illumina HiSeq 4000 yielding 2.1 Tbases of raw read data.

## Analysis Pipeline

See `/markdowns` for more detail.

### Raw read processing

* QC of raw reads
* Adapter and base-quality trimming
* QC post-trimming

Raw reads were quality controlled using `FastQC` and sample statistics were aggregated using `MultiQC`.
Adapters were trimmed from raw reads and low quality bases clipped using `Trimmomatic-v0.39`.
A second round of QC was performed to check proper trimming behavior.

### Alignment

* Aligner Testing
* Alignment
* Mapping QC for each sequencing batch

Different alignment software was benchmarked using `Teaser`. `bwa mem` was shown to give the best compromise between percent mapping reads, mapping quality and runtime.
`bwa mem` was used to align reads from all samples to the *Fundulus heteroclitus* reference genome version 4.1.
`Qualimap2` was used to assess alignment statistics for each sequencing batch.

### Alignment Polishing

* Merging reads across sequencing batches
* Removing duplicate reads
* Clipping overlapping paired-end reads
* Post-filter mapping QC

Reads were given read-group information and merged across sequencing batches resulting in one `bam` file per sample.
Duplicate reads were removed using `Picard Tools 1.103` `MarkDuplicates` and overlapping paired-end reads were clipped using the `bamUtil 1.0.14` `clipOverlap` tool to avoid overestimating sequence support at the SNP calling stage.
Final alignment QC was performed with `Qualimap2` to generate per-sample depth and coverage statistics.

Note:
* BQSR? - GATK recommends it but papers don't do it. Need known variants file, which needs to be bootstrapped. Could do a second round if needed, keep going for now.
* Indel realignment? - No, because will use BAQ option in ANGSD (use -baq 2 option!, higher sensitivity)
* Mapping Quality Filter? - No, do it at ANGSD variant calling step. Use cutoff depending on histogram of MQs (cutoff of 20 or 25 is good)
* Remove high depth locations? - No, do it at ANGSD variant calling step

### Variant Calling

Variant calling was performed using 

* ANGSD: dont forget to: Use BAQ option 2, filter for MQ, filter for high depth
