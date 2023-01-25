# Variant Calling

## Depth Distribution

For variant calling several filters are applied to the BAM files. Some of these filters take into account low or excess sequencing depth.
In order to inform the filtering approach for variant calling the depth distribution across all samples and genomic positions is generated using `ANGSD`.
Here, we only assess depth statistics on the 24 *F. heteroclitus* chromosomes to avoid bias in the unplaced scaffolds.
In order to speed up processing, the analysis is split by chromosome and one instance of `ANGSD` is run per chromosome:
```
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

## Raw SNP Calling

The following `ANGSD` script was run for calling SNPs across the entire *F. heteroclitus* genome.
```
angsd -bam ../master_bamlist.txt -ref $REF -out $CHROM'.raw' -r $CHROM':' \
-doCounts 1 -doMajorMinor 1 -doMaf 1 -GL 1 -doGlf 2 -doSnpStat 1 -doHWE 1 \
-minMapQ 20 -baq 2 -minQ 20 \
-minInd 478 -setMaxDepth 1600 \
-SNP_pval 1e-6 -minMaf 0.05 -rmTriallelic 1e-4 \
-P 4 #ANGSD only takes a maximum of 8 threads. Due to I/O operations being the bottleneck.
```

See below for a more thorough description/explanation of the parameters chosen:

### Processing Parameters:

* `-GL 1` (`SAMTools` genotype likelihood model)
* `-doGlf 2` (Outputs GLs in `BEAGLE` genotype likelihood format. Best to use for downstream analysis.)
* `-doMajorMinor 1` (Infers major and minor alleles from data using GLs)
* `-doMaf 1` (Calculates MAF from GLs)
* `-doCounts 1` (Counts reads per site. Required for `-setMaxDepth`)
* `-doSnpStat 1` (Produces bias statistics per SNP which will be used for subsequent filtering.)
* `-doHWE 1` (Calculates HWE stats. Required for `-doSnpStat` and can be used for later filtering.)

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
* `-rmTriallelic 1e-4` (Remove sites with a significant chance of being triallelic. Threshold commonly used in the literature.)

## Variant Bias Filtering

After the raw set of SNPs has been called it is further filtered to remove any allelic biases. These biases include:

* Strand bias           - whether one or the other allele is preferentially called on the forward or reverse strand.
* Base Quality Bias     - whether one or the other allele has systematically higher base qualities.
* Mapping Quality Bias  - whether the reads of one or the other allele have systematically higher mapping qualities.
* Position Bias         - whether one or the other allele is systematically found near the ends of reads.

Bias statistics were printed in the `.snpStat.gz` file. The distributions of bias statistics were plotted in `R` and filtering cutoffs chosen according to the tails. Here the specific filtering cutoffs were to keep SNPs with:

* Strand Bias (GATK Implementation)         < 0.6
* Strand Bias (Fisher Score, Phred-scaled)  < 50
* Base Quality Bias (Phred-scaled)          < 50
* Mapping Quality Bias (Phred-scaled)       < 20
* Position Bias (Phred-scaled)              < 75

The list of filtered SNPs was then used to filter the raw beagle GL and MAF files for further downstream processing using the following `mawk` script:
```
cat <(zcat $CHROM'.raw.beagle.gz' | head -n 1) <(zcat $CHROM'.raw.beagle.gz' | mawk 'NR==FNR{array[$1"_"$2];next} $1 in array' $CHROM'.filt.sites' -) | gzip > $CHROM'.filt.beagle.gz'
cat <(zcat $CHROM'.raw.mafs.gz' | head -n 1) <(zcat $CHROM'.raw.mafs.gz' | mawk 'NR==FNR{array[$1,$2];next} ($1,$2) in array' $CHROM'.filt.sites' -) | gzip > $CHROM'.filt.mafs.gz'
```

## LD Pruning

Pairwise LD statistics between SNPs were obtained using `ngsLD` with GLs as input:
```
N_SITES=$(wc -l $CHROM'.filt.sites')

ngsLD --geno $CHROM'.filt.beagle.gz' --probs --n_ind 956 --n_sites $N_SITES --pos $CHROM'.filt.sites' \
--max_kb_dist 10 --min_maf 0 \
--n_threads 10 \
| cut -f 1,4,7,8,9,10,11 | gzip > $CHROM'.filt.ld.gz'
```
A `--max_kb_dist 10` was chosen since previous *F. heteroclitus* dataset did not show significant linkage beyond 10Kb.
A `--min_maf 0` was chosen to speed up computation since sites had been MAF filtered previously already.
The final `cut` call removes columns containing base information that interferes with downstream analysis.

Before plotting LD decay, LD output files were randomply downsampled to ~100 million pairwise comparisons in order to speed up calculations. Anything above 100,000 pairwise comparisons gives a decent distribution. To downsample the following `mawk` script was used:
```
zcat chr.filt.ld.gz | mawk '{if rand() <= 0.01) print $0}' | gzip > chr.filt.sample1p.ld.gz
```
Next LD decay was plotted using the `fit_LDdecay.R` script supplied with `ngsLD`:
```
Rscript --vanilla --slave ~/software/local/ngsLD/scripts/fit_LDdecay.py --ld_files input.txt --out chr.filt.sample1p.ld.plot \
--n_ind 956 --ld r2 --recomb_rate 2.34 --fit_boot 1000 --fit_bin_size 50 --fit_level 2
```
The LD decay curve informed reasonable cutoff values for LD pruning. More specifically, SNPs with r2 values above 0.1 were considered to be in linkage. Strongest linked sites were identified using a `Python` script supplied with `ngsLD`:
```
python3 ~/software/local/ngsLD/scripts/prune_ngsLD.py \
--input $CHROM'.filt.ld.gz' --output $CHROM'.filt.lnkd' \
--max_dist 10000 --min_weight 0.1 --keep_heavy
```
 **NOTE:** The `--keep_heavy` option was used to output strongly linked sites rather than a pruned dataset. This is due to the `.ls.gz` files not containing all SNP information since only SNPs in proximity were considered for LD calculation. Distant SNPs do not feature in the `.ld.gz` files and will therefore not be retained although they are unlinked. For the same reason, scaffolds with only 1 SNP have empty LD files. These SNPs should however be included in the final SNP set.
 
 To generate a pruned SNP set, the linked sites were removed from the filtered dataset:
