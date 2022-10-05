# Selection among microhabitats using low-coverage WGS

This repository contains the scripts used in the analysis of within-generation selection among *Fundulus heteroclitus* individuals inhabiting distinct microhabitats. Low-coverage whole genome sequencing (LC-WGS) was performed on ~1K individuals sampled in distinct microhabitatss and at different time points to estimate allele frequency changes within a single generation as a result of microhabitat selection. Allele frequency changes were also associated with previously documented phenotypic divergence among microhabitat residents. Phenotypic divergence among microhabitat residents was previously published here (LINK) whereas the results of this LC-WGS project are published in XXXX.

The repository contains markdown files in `/markdowns` detailing the analysis process and parameter justifications. `/hps_scripts` contain code with minimal annotation used to run the analysis on HPC servers at the University of Miami.

## Raw Data

Libraries were prepared from DNA isolated from 963 *Fundulus heteroclitus* individuals. Each individual was uniquely indexed and the library pool was sequenced in a two-stage manner in order to normalize sequencing depth across individuals. A total of 19 lanes were sequenced on the Illumina HiSeq 4000 yielding 2.1 Tbases of raw read data.

## Brief Summary of Analysis Pipeline

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

Different alignment software was benchmarked using `Teaser`. `Bowtie2` in `--sensitive-local` mode was shown to give the best compromise between percent mapping reads and runtime.

### BAM Filtering

* Non-unique mapping
* Mapping Quality
* Dedup

### Seq Run Aggregation

* Merging
* Dedup again

### BAM Cleaning

* Clip overlap
* Indel realignment
* Remove high depth
