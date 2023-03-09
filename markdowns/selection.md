# Selection
To detect loci under selection in both pond and basin environments, samples form two time points are compared and loci are tested for a significant difference in allele frequency/counts.

## Genetic Association
Here we use the genetic association tool in the `ANGSD` software suite to test for significant allelic differences between time points. By considering the two time points as binary "phenotypes" we can detect loci that are significantly associated with one or the other time point and must have therefore undergone a significant shift from one to the next.
Furthermore, we can include relevant covariates such as sampling year and the individual pond ID in the analysis.

We use the `ANGSD` association tool and test only high-quality, filtered variant sites previously determined.
We employ a genotype dosage model in a generalized linear framework (`-doAsso 6`) which has the advantage of avoiding calling hard genotypes that would be problematic with low-coverage data. Instead, the expected genotype "dosage" is calculated from the genotype probabilities. A dosage can therefore fall "between" two genotypes.
Genotype dosage is calculated as follows:
```math
E[G|X] = p(G=1|X) + 2p(G=2|X)
```
Dosage genotypes are therefore continuous variables that can be used in logistic regression. It also has the advantage that effect sizes are calculated.

The code for running the generalized linear dosage model in `ANGSD` to calculate significant allelic shifts in the ponds is given below:
```
angsd -bam ../masterFunhe.bamlist -ref $REF -r $CHROM -sites '../variants/'$CHROM'.filt.sites' -out $CHROM'.pond.AFC' \
-GL 1 -doMajorMinor 3 -doMaf 1 -doPost 1 -doAsso 6 \
-sampleFile afc.sample -whichPhe pondBin -whichCov year,pop -Pvalue 1 \
-minMapQ 20 -baq 2 -minQ 20 \
-P 4 #ANGSD only takes a maximum of 8 threads. Due to I/O operations being the bottleneck.
```
The time points as well as covariates are supplied via a `.sample` file which has the specific format:
```
ID basinBin pondBin year pop
0 B B D D
16FBaX_000001 1 -999 0 0
18FBaP_G2O4P6 1 -999 1 0
18FP3P_P2R3R5 -999 1 1 3
18SBaX_G5R7X0 0 -999 1 0
18SP4X_000028 -999 0 1 4
.
.
.
```
The header line contains the sample ID, phenotype names (in this case binary data for each time point separated by ponds and basin) and covariate names. The second line encodes the type of data (ID, continuous, discrete). Discrete data are encoded numerically and missing values marked as -999.

**NOTE: The genetic association test among timepoints yielded inflated p-values and was not used. Instead the CMH approach (detailed below) was employed.**

## Cochran-Mantel-Haenszel tests
The CMH test runs on stratified contingency tables and detects SNPs where the direction of allele frequency change is consistent across replicates. Here we use a series of bash scripts as well as the `mantelhaenszel.test` in `R` to calculate the chi-squared distributed test statistic. We then correct the test statistic according to the genomic inflation factor to yield corrected p-values. In a final step the p-values are FDR corrected.

First calculate allele frequency estimates for each population/time point using `ANGSD`:
```
#Extract population bamlist and chromosome from file using job index.
POP=$(awk -v jindex=$LSB_JOBINDEX 'NR==jindex {print $1}' ../../pop.chr.list)
CHROM=$(awk -v jindex=$LSB_JOBINDEX 'NR==jindex {print $2}' ../../pop.chr.list)

#CODE
angsd -bam '../../'$POP'.bamlist' -ref $REF -r $CHROM -sites '../../variants/'$CHROM'.filt.sites' -out $POP'.'$CHROM \
-GL 1 -doMajorMinor 3 -doMaf 1 \
-minMapQ 20 -baq 2 -minQ 20 \
-P 4 #ANGSD only takes a maximum of 8 threads. Due to I/O operations being the bottleneck.
```
Next, we calculate the most-likely estimate of allele counts from the MAF for both the major and minor allele in all populations and time points. Note that these estimates are not integer values since they are based on genotype likelihoods. The allele counts are then passed to an Rscript for running the `mantelhaenszel.test`.
```
#Extract chromosome from file using job index.
CHROM=$(awk -v jindex=$LSB_JOBINDEX 'NR==jindex {print $0}' chrScaff.chr)

#Set number of strata (replicates) in the CMH tests
NUMSTRATA=6

#Estimate allele counts from GL based MAFs
zcat '16SP1.'$CHROM'.mafs.gz' '16FP1.'$CHROM'.mafs.gz' '16SP2.'$CHROM'.mafs.gz' '16FP2.'$CHROM'.mafs.gz' '16SP3.'$CHROM'.mafs.gz' '16FP3.'$CHROM'.mafs.gz' '18SP1.'$CHROM'.mafs.gz' '18FP1.'$CHROM'.mafs.gz' '18SP3.'$CHROM'.mafs.gz' '18FP3.'$CHROM'.mafs.gz' '18SP4.'$CHROM'.mafs.gz' '18FP4.'$CHROM'.mafs.gz' | \
mawk -v OFS="\t" '!/chromo/ {printf("%.3g\n%.3g\n", $7*(1-$6), $7*$6)}' | gzip > 'pond.'$CHROM'.mm.ac.gz'

#Get SNP count from positions file
SNPCOUNT=$(wc -l '../../variants/'$CHROM'.filt.pos' | cut -d ' ' -f 1)

#Run Rscript to calculate CMH statistics
#Requires input variables (snpcount, numstrata, infile, outfile) IN ORDER!
Rscript /home/mae120/scripts/microhabSel_LC-WGS/clusterCMH.R $SNPCOUNT $NUMSTRATA 'pond.'$CHROM'.mm.ac.gz' 'pond.'$CHROM'.cmh.gz'

#Merge position and CMH test data
paste <(cat <(echo -e 'chrom\tpos') '../../variants/'$CHROM'.filt.pos') <(zcat 'pond.'$CHROM'.cmh.gz') | gzip > 'pond.'$CHROM'.cmh.pos.gz'
rm 'pond.'$CHROM'.cmh.gz'
```
The Rscript used to calculate the CMH statistic can be found below. A Woolf test is run prior to calculating the CMH statistic in order to evaluate the homogeneity of odds ratios i.e. if the odds ratios (effects) among strata/replicates are consistent. Since the CMH test assumes homogeneity, any SNPs failing the Woolf test are masked.
```
#CMH-test including check for homogeneity of odds-ratios
#Required inputs are chromosome, number of snps, number of strata/replicates for CMH, input and output files.
#Input file is a list of interleaved allele counts (major/minor) for all SNPs, both timepoints and all replicates in this order.
#Load libraries
library(data.table)

#Parse arguments
args <- commandArgs(trailingOnly = TRUE) 

snpcount <- args[1]
numstrata <- args[2]
infile <- args[3]
outfile <- args[4]

message("\nRscript input parameters:\nNumber of SNPs: ", snpcount,
        "\nNumber of replicates: ", numstrata,
        "\nInput file: ", infile,
        "\nOutput file: ", outfile,
        "\n")

#Define Functions
woolf.test <- function(x) {
  DNAME <- deparse(substitute(x))
  if (any(x == 0))
    x <- x + 1 / 2
  k <- dim(x)[3]
  or <- apply(x, 3, function(x) (x[1,1] * x[2,2]) / (x[1,2] * x[2,1]))
  w <-  apply(x, 3, function(x) 1 / sum(1 / x))
  o <- log(or)
  e <- weighted.mean(log(or), w)
  STATISTIC <- sum(w * (o - e)^2)
  PARAMETER <- k - 1
  PVAL <- 1 - pchisq(STATISTIC, PARAMETER)
  METHOD <- "Woolf-test on Homogeneity of Odds Ratios (no 3-Way assoc.)"
  names(STATISTIC) <- "X-squared"
  names(PARAMETER) <- "df"
  structure(list(statistic = STATISTIC, parameter = PARAMETER,
                 p.value = PVAL, method = METHOD, data.name = DNAME, observed = o,
                 expected = e), class = "htest")
}

cmh.test <- function(a)
{
  woolf <- apply(a, 2, function(snp) woolf.test(snp)$p.value) #function applied across 2nd dimension (SNP dimension) of array

  cmh <- ifelse(woolf < 0.05,
                NA,
                apply(a, 2, function(snp) mantelhaen.test(snp, correct = F)) #function applied across SNP dimension (2nd)
                )

  out <- cbind(matrix(unlist(lapply(cmh, `[`, c("statistic","p.value","estimate"))),
                      ncol = 3,
                      byrow = T),
               woolf)
  
  colnames(out) <- c("chisq",
                     "p.value",
                     "commonOR",
                     "woolf.p")
  
  return(out)
}

#Calculation
message("Reading allele count data into array.")
ac.array <- array(scan(infile),
                     dim = c(2, snpcount, 2, numstrata)) #Dimensions are: major/minor, snp, fall/spring, pop/strata

message("\nStarting CMH tests.\n")
cmh <- cmh.test(ac.array)

message("Finished CMH calculations. Writing to file.\n")
fwrite(data.table(cmh), outfile, na = "NA", sep = "\t", quote = F)
message("Completed!")
```
Genomic correction using the inflation factor.


## Phasing allele frequency shifts
```
paste <(zcat 16SP3.[Ns]*.mafs.gz) <(zcat 16FP3.[Ns]*.mafs.gz) | mawk -v OFS="\t" 'NR==1 {print $1,$2,$3,$4,$5,"springMAF","springN","fallMAF","fallN","delta","year","pop"} !/chromo/ {print $1,$2,$3,$4,$5,$6,$7,$13,$14,$13-$6,"2016","Pond3"}' | gzip > 16P3.delta.gz
```
```

```
