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
-minMapQ 20 -minQ 20 \                                  # Filters on BAM files: minimum mapping quality and base quality.
-P 4                                                    # ANGSD only takes a maximum of 8 threads. Due to I/O operations being the bottleneck.
```
A custom `mawk` script was used to parse the output of `-dumpCounts 2` (`$CHROM.counts.gz` files) to obtain the number of individuals with at least one read for each site. These `$CHROM.counts.gz` files contain read counts per individual. Header lines were ignored, and the number of individuals with at least 1 read were counted for every genomic position.
```
CHROM=$(awk -v jindex=$LSB_JOBINDEX 'NR==jindex {print $0}' chr_only.txt)

zcat $CHROM'.counts.gz' | tail -n +2 | mawk '{c=0;for (i = 1; i <= NF; ++i) if($i>0){++c}{print c}}' | gzip > $CHROM'.indPerSite.gz'
```
Both the global depth distribution (`XXX.depthGlobal`) and the distribution of the number of individuals with at least 1 read were plotted in `R`. By inspecting the inflection points of the distribution, the maximum depth per site and the minimum number of individuals per site for variant calling were extracted. 

## SNP Calling






## Do:

-GL 1 (in most publications)

-doMajorMinor 1 (infers from data itself using GLs)

-doGlf 2 (output Beagle format GL file, most common)

-doMaf 1 (calculate MAF. use 1 just because Therkildsen uses it, could also do 2, a combination of both or even add 8)

-doCounts 1 (required for -setMaxDepth)

## Alignment Filters:

-minQ 20 (from literature)

-minMapQ 20 (from literature)

-baq 2 (equivalent to SAMTools)

## Depth Filters

**NOT** -setMinDepth (think about it, maybe use depth distribution? Probably not necessary since I have minInd which already sets a lower bound!)

-setMaxDepth 1600 (Used depth distribution, fit normal curve to uppper end then +3 s.d.)

-minInd 478 (This already defines the minDepth since there is 1 read per ind. Used depth distribution inflection point for cutoff = almost exactly half of individuals. Used 50% of individuals as cutoff.)

## SNP Filters

-SNP_pval 1e-6 (from literature)

-minMaf 0.05 (think about it. Maybe not do minMaf but rather minimum number of alternative allele counts! i.e. a MAF of 0.01 is fine if it is supported by 10 individuals of 1000. Change according to minInd!! MAF = 2/MinInd (to have at least 2 alternative alleles...). Later saw that there was sufficient coverage even fo 1% MAF but still went with 5% to be conservative. Also 1% SNPs do not contain much information and swamp the dataset.)

-no HWE filter due to strong selection!

DOWNSTREAM

DONT USE REGIONS FILE FORMAT FOR SITE SELECTION - SLOWWWW

Use the -sites argument and supply the relevant chromosome using -r. Should be faster.
