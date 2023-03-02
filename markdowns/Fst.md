# *Fst* Calculation
*Fst* calculation using `ANGSD` requires the site allele frequencies which can then be used to produce site frequency spectra which in turn can be used to calculate *Fst*. For large datasets, input data can be reduced to only polymorphic sites and the 2D-SFS should be estimated for each chromosome as to not run into memory issues.

The procedure as a whole is:
* Use `ANGSD` to call polymorphic sites. This MUST be done using ALL samples and with *lenient* SNP filtering in order to keep rare variants i.e. NO MAF FILTER!
* Use `ANGSD` to calculate `.saf` (site allele frequency) files for each population and for each chromosome using the polymorphic sites from before
* Use the `winsfs` (faster than `realSFS`) software to calculate the 2D-SFS for each population pair/comparison and for each chromosome
* Use the above calculated per-chromosome 2D-SFS as a prior jointly with the `.saf` files from step 2 to calculate chromosomal *Fst*. Chromosomal *Fst* can then be averaged or used as a distribution.

## Call Polymorphic Sites
*Fst* calculations do not require information from monomorphic sites in the SFS. Hence, processing can be sped up significantly by only considering polymorphic sites. **NOTE** Rare variants are still important! Use lenient SNP filtering criteria (e.g. -SNP_pval 0.1) to keep sites that may be polymorphic.
```
angsd -bam ../masterFunhe.bamlist -ref $REF -r $CHROM -out $CHROM'.poly' \
-doMajorMinor 1 -doMaf 1 -GL 1 \
-minMapQ 20 -baq 2 -minQ 20 \
-SNP_pval 0.1 \
-P 4 #ANGSD only takes a maximum of 8 threads. Due to I/O operations being the bottleneck.
```
After calling polymorphic sites, generate `.sites` files from the `.maf.gz` output files and index them using `ANGSD`.

## SAF Calculation
Once polymorphic sites are determined, calculate the site allele frequencies `SAF` for each population and chromosome. For this, prepare an input file with population in the first column and chromosome in the second.
```
#Extract populations and chromosomes from file using job index.
POP=$(awk -v jindex=$LSB_JOBINDEX 'NR==jindex {print $1}' pop.chr.list)
CHROM=$(awk -v jindex=$LSB_JOBINDEX 'NR==jindex {print $2}' pop.chr.list)

angsd -bam '../'$POP'.bamlist' -ref $REF -anc $REF -r $CHROM -sites $CHROM'.poly.sites' -out $POP'.'$CHROM \
-GL 1 -doSaf 1 \
-minMapQ 20 -baq 2 -minQ 20 \
-P 4 #ANGSD only takes a maximum of 8 threads. Due to I/O operations being the bottleneck.
```
## 2D-SFS and *Fst* Calculation
The per-population `.saf` files can be combined into a two-dimensional site frequency spectrum (2D-SFS) from which *Fst* is calculated. The `winsfs` tool is faster than the `realSFS` script distributed with `ANGSD`. Prepare an input file with all pairwise population comparisons and for every chromosome.
```
#Extract populations and chromosomes from file using job index.
POP1=$(awk -v jindex=$LSB_JOBINDEX 'NR==jindex {print $1}' pairwise.chr.comps)
POP2=$(awk -v jindex=$LSB_JOBINDEX 'NR==jindex {print $2}' pairwise.chr.comps)
CHROM=$(awk -v jindex=$LSB_JOBINDEX 'NR==jindex {print $3}' pairwise.chr.comps)

#Calculate the 2DSFS for each pairwise comparison and for each chromosome
winsfs -v -t $LSB_DJOB_NUMPROC $POP1'.'$CHROM'.saf.idx' $POP2'.'$CHROM'.saf.idx' > $POP1'.'$POP2'.'$CHROM'.sfs'

#Calculate Fst
winsfs stat -v -s fst $POP1'.'$POP2'.'$CHROM'.sfs' > $POP1'.'$POP2'.'$CHROM'.fst'
```
Once chromosome-level *Fst* estimates have been calculated for each population pair, they can be averaged to obtain a global *Fst* for each comparison:
```
while read pop1 pop2; do cat $pop1.$pop2.NC*.fst | awk -v pop1=$pop1 -v pop2=$pop2 '{sum+=} END {print pop1,pop2,sum/NR}'; done < pairwise.comps > pairwise.fst
```
