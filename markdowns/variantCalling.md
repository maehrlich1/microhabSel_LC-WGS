# Variant Calling

## Depth Distribution

For variant calling several filters are applied to the BAM files. Some of these filters take into account low or excess sequencing depth.
In order to inform the filtering approach for variant calling the depth distribution across all samples and genomic positions is generated using `ANGSD`.
Here, we only assess depth statistics on the 24 *F. heteroclitus* chromosomes to avoid bias in the unplaced scaffolds.
In order to speed up processing, the analysis is split by chromosome and one instance of `ANGSD` is run per chromosome:
```
CHROM=$(awk -v jindex=$LSB_JOBINDEX 'NR==jindex {print $0}' chr_only.txt)

angsd -bam master_bamlist.txt -out master -r $CHROM: \  # Input and output files and the chromosome to be processed.
-doCounts 1 -dumpCounts 2 -doDepth 1 -maxDepth 3500 \   # Calls to generate read counts and depth statistics. dumpCounts 2 produces per individual read counts which is needed for subsequent plotting. maxDepth of 2x expected depth seems like a good target
-minMapQ 20 -baq 2 -minQ 20 \                           # Filters on BAM files: minimum mapping quality and base quality and adjusted base quality scrores.
-P 4                                                    # ANGSD only takes a maximum of 8 threads. Due to I/O operations being the bottleneck.
```
A custom `mawk` script was used to parse the output of `-dumpCounts 2` (`$CHROM.counts.gz` files) to obtain the number of individuals with at least one read for each site. These `$CHROM.counts.gz` files contain read counts per individual. Header lines were ignored, and the number of individuals with at least 1 read were counted for every genomic position.
```
CHROM=$(awk -v jindex=$LSB_JOBINDEX 'NR==jindex {print $0}' chr_only.txt)

zcat $CHROM'.counts.gz' | tail -n +2 | mawk '{c=0;for (i = 1; i <= NF; ++i) if($i>0){++c}{print c}}' | gzip > $CHROM'.indPerSite.gz'
```
Both the global depth distribution (`XXX.depthGlobal`) and the distribution of the number of individuals with at least 1 read were plotted in `R`. By inspecting the inflection points of the distribution, the maximum depth per site and the minimum number of individuals per site for variant calling were extracted. 

## SNP Calling

The following `ANGSD` script was run for calling SNPs across the entire *F. heteroclitus* genome.
```
angsd -bam ../master_bamlist.txt -ref $REF -out $CHROM'.raw' -r $CHROM':' \
-doCounts 1 -doMajorMinor 1 -doMaf 1 -GL 1 -doGlf 1 \
-minMapQ 20 -baq 2 -minQ 20 \
-minInd 478 -setMaxDepth 1600 \
-SNP_pval 1e-6 -minMaf 0.05 -rmTriallelic 0.01 \
-P 4 #ANGSD only takes a maximum of 8 threads. Due to I/O operations being the bottleneck.
```

See below for a more thorough description/explanation of the parameters chosen:

### Processing Parameters:

* `-GL 1` (`SAMTools` genotype likelihood model)
* `-doGlf 1` (Outputs GLs in binary, log10-scaled genotype likelihood format. Best to use for downstream analysis.)
* `-doMajorMinor 1` (Infers major and minor alleles from data using GLs)
* `-doMaf 1` (Calculates MAF from GLs)
* `-doCounts 1` (Counts reads per site. Required for `-setMaxDepth`)

### Alignment Filters:

* `-minQ 20` (Base quality filter)
* `-minMapQ 20` (Mapping quality filter)
* `-baq 2` (`SAMTools` model for recalibrating base quality scores around indels. Similar to indel realignment.)

### Depth Filters:

* `-setMaxDepth 1600` (Inferred from the global depth distribution by fitting a normal curve to the data and calculating +3 standard deviations.)
* `-minInd 478` (Inferred from the distribution of individual coverage i.e. the number of individuals with at least 1 read per site. Inflection point of the lower tail used as a cutoff. This was almost exactly 50% of individuals which was used as the cutoff.)

**NOTE:** `-setMinDepth` was **NOT** used as a depth filter since `-minInd` already sets a lower bound.

### SNP-level Filters:

* `-SNP_pval 1e-6` (To define polymorphic sites.)
* `-minMaf 0.05` (Would like at least 5 reads confirming an alternative allele. Given `-minInd 478` this allows for a MAF filter of 1%. However, a 5% cutoff was chosen as a more conservative measure.)
* `-rmTriallelic 0.01` (Remove sites with a significant chance of being triallelic.)

## Variant Bias Filtering

After the raw set of SNPs has been called it is further filtered to remove any allelic biases. These biases include:

* Strand bias           - whether one or the other allele is preferentially called on the forward or reverse strand.
* Base Quality Bias     - whether one or the other allele has systematically higher base qualities.
* Mapping Quality Bias  - whether the reads of one or the other allele have systematically higher mapping qualities.
* Position Bias         - whether one or the other allele is systematically found near the ends of reads.

Bias statistics were calculated using the `-doSnpStat` option in `ANGSD` using the following code. **NOTE:** Could and probably should use GL files as input for faster processing:

```
angsd -bam ../master_bamlist.txt -rf $CHROM'.raw.regions' -sites $CHROM'.raw.sites' -out $CHROM'.raw' \
-doMajorMinor 3 -GL 1 -doHWE 1 -doSnpStat 1 \
-P 4 #ANGSD only takes a maximum of 8 threads. Due to I/O operations being the bottleneck.
```
The distributions of bias statistics were subsequently plotted in `R` and filtering cutoffs chosen according to the tails. Here the specific filtering cutoffs were to keep SNPs with:

* Strand Bias (GATK Implementation)         < 0.6
* Strand Bias (Fisher Score, Phred-scaled)  < 50
* Base Quality Bias (Phred-scaled)          < 50
* Mapping Quality Bias (Phred-scaled)       < 20
* Position Bias (Phred-scaled)              < 75

The filtered SNP set was then used to produce filtered GL and MAF files in `ANGSD` for further downstream processing.

## LD Pruning

Pairwise LD statistics between SNPs were obtained using `ngsLD` with GLs as input.
