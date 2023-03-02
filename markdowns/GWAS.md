# GWAS
Here we use the genetic association tool in the ANGSD software suite to test for significant association between several measured phenotypes and genotype. We test only high-quality, filtered variant sites previously determined. We employ a genotype dosage model in a generalized linear framework (-doAsso 6) which has the advantage of avoiding calling hard genotypes that would be problematic with low-coverage data. Instead, the expected genotype "dosage" is calculated from the genotype probabilities. A dosage can therefore fall "between" two genotypes. Genotype dosage is calculated as follows:
```math
E[G|X] = p(G=1|X) + 2*p(G=2|X)
```
Dosage genotypes are therefore continuous variables that can be regressed against phenotypic data. It also has the advantage that effect sizes are calculated.

The code for running the generalized linear dosage model in ANGSD to obtain loci significantly associated with phenotypes is as follows:
```
angsd -bam ../masterFunhe.bamlist -ref $REF -r $CHROM -sites '../variants/'$CHROM'.filt.sites' -out $CHROM'.uninfGWAS' \
-GL 1 -doMajorMinor 3 -doMaf 1 -doPost 1 -doAsso 6 \
-sampleFile all_residuals.txt -whichPhe WAM.resid,ASR.resid,CTmax.resid,CaMGLU.resid,CaMFA.resid,CaMLKA.resid,CaMEND.resid  -Pvalue 1 \
-minMapQ 20 -baq 2 -minQ 20 \
-P 4 #ANGSD only takes a maximum of 8 threads. Due to I/O operations being the bottleneck.
```
Phenotypic data is supplied via the `all_residuals.txt` file which contains the residuals after removing technical variance. It has the specific format:
```
ID WAM.resid ASR.resid CTmax.resid CaMGLU.resid CaMFA.resid CaMLKA.resid CaMEND.resid
0 P P P P P P P
16FBaX_000001 -999 -999 -999 -999 -999 -999 -999
18FBaP_G3G4O6 -0.231158185346304 0.307656177 0.354084883770639 8.91413212327556 0.45526361595315 7.22587671677027 -2.51510149701529
18FP1P_P1R2P4 -0.0428932675934757 -0.330315175 0.257989803536866 -999 -999 -999 -999
18FP1X_R1G2P5 -999 -999 -999 4.93308436042313 0.207800575379817 2.08098906131342 -3.16506786852758
18SP4X_000011 -999 -999 -999 -999 -999 -999 -999
.
.
.
```
The header line contains the sample ID and phenotype names. The second line encodes the type of data (ID, continuous, discrete). Discrete data are encoded numerically and missing values marked as -999. Note that even individuals with no phenotypic data are included since these heavily inform genotype likelihoods and probabilities.
