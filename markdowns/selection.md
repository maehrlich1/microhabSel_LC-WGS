# Selection
To detect loci under selection in both pond and basin environments, samples form two time points are compared and loci are tested for a significant difference in allele frequency/counts.
Here we use the genetic association tool in the `ANGSD` software suite to test for significant allelic differences between time points. By considering the two time points as binary "phenotypes" we can detect loci that are significantly associated with one or the other time point and must have therefore undergone a significant shift from one to the next.
Furthermore, we can include relevant covariates such as sampling year and the individual pond ID in the analysis.

## Testing selection through temporal association
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
