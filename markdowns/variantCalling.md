# Variant Calling

## Depth Distribution

For variant calling several filters are applied to the BAM files. Some of these filters take into account low or excess sequencing depth.
In order to inform the filtering approach for variant calling the depth distribution across all samples and genomic positions is generated using `ANGSD`.
Here, we only assess depth statistics on the 24 *F. heteroclitus* chromosomes to avoid bias in the unplaced scaffolds.
In order to speed up processing, the analysis is split by chromosome and one instance of `ANGSD` is run per chromosome:
```
angsd -bam master_bamlist.txt -out master -r $CHROM: \  # Input and output files and the chromosome to be processed.
-doCounts 1 -doDepth 1 -maxDepth 100000 -dumpCounts 1 \ # Calls to geenrate allele counts and depth statistics.
-minMapQ 20 -minQ 20 \                                  # Filters on BAM files: minimum mapping quality and base quality.
-P 4                                                    # ANGSD only takes a maximum of 8 threads. Due to I/O operations being the bottleneck.
```
Depth output files were subsequently plotted in `R`.


## SNP Calling






## To Do:

-GL 1 (in most publications)

-doMajorMinor 1 (infers from data itself using GLs)

-doGlf 2 (output Beagle format GL file, most common)

-doMaf 1 (calculate MAF. use 1 just because Therkildsen uses it, could also do 2, a combination of both or even add 8)

-doCounts 1 (required for doDepth and dumpCounts)

-dumpCounts 1 (prints overall depth per locus - maybe use for fine-tuning parameters?!)

-doDepth 1 (gives depth distribution for every sample)

## Filters:

-setMinDepth (think about it, maybe use depth distribution?)

-setMaxDepth (think about it. depth distribution??)

-minInd (think about it)

-minQ (probably 20)

-minMapQ (also probably 20)

-SNP_pval 1e-6 (seems like what everybody uses)

-minMaf (think about it. Maybe not do minMaf but rather minimum number of alternative allele counts! i.e. a MAF of 0.01 is fine if it is supported by 10 individuals of 1000)

-baq 2 (I decided on 2)
