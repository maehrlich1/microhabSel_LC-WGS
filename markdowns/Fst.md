# *Fst* Calculation
*Fst* calculation using `ANGSD` requires the site allele frequencies which can then be used to produce site frequency spectra which in turn can be used to calculate *Fst*. For large datasets the 2D-SFS should be estimated for each chromosome as to not run into memory issues.

The procedure as a whole is:
* Use `ANGSD` for calculating `.saf` files for each population separately
* Use the `realSFS` script found in `ANGSD/misc/realSFS` to calculate 2D-SFS for each population pair and for each chromosome
* Average across chromosomes to get the global estimate of the 2D-SFS
* Use the above calculated 2D-SFS as a prior jointly with all `.saf` files from step 1 to calculate *Fst* binary files
* Extract the *Fst* values from the *Fst* file

## SAF Calculation
Since the SFS requires the relative proportion of monomorphic sites too, SAF calculation requires information from the entire genome and not just SNPs.
Therefore raw `BAM` files are once again provided as input:
```
angsd -bam '../'$POP'.bamlist' -ref $REF -anc $REF -out $POP \
-GL 1 -doSaf 1 \
-minMapQ 20 -baq 2 -minQ 20 \
-P 8 #ANGSD only takes a maximum of 8 threads. Due to I/O operations being the bottleneck.
```

## 2D-SFS and *Fst* calculation
The following script executes the remaining steps:
```
while read CHROM
do
	~/software/local/angsd/misc/realSFS -fold 1 -r $CHROM -cores $LSB_DJOB_NUMPROC $POP1'.saf.idx' $POP2'.saf.idx' >> $POP1'.'$POP2'.chr.sfs'
done < ../chr.chr

#Average over all chromosomes to get the global 2D-SFS
mawk '{for(i=1; i<=NF; i++) {a[i]+=$i}}; END {for(i=1; i<=NF; i++) printf "%s%s", a[i]/NR, (i==NF?ORS:OFS)}' $POP1'.'$POP2'.chr.sfs' > $POP1'.'$POP2'.sfs'

#Prepare the Fst output file
~/software/local/angsd/misc/realSFS fst index -fold 1 -cores $LSB_DJOB_NUMPROC $POP1'.saf.idx' $POP2'.saf.idx' -sfs $POP1'.'$POP2'.sfs' -fstout $POP1'.'$POP2

#Get the global Fst estimate
~/software/local/angsd/misc/realSFS fst stats $POP1'.'$POP2'.fst.idx' > $POP1'.'$POP2'.fst'
```

**Note:** The `while` command loops over all chromosomes since reading in the entire genome requires too much memory. One 2D-SFS is therefore produced *per chromosome*. These 2D-SFS estimates are output as flattened matrices on a single line. The output of the `while` loop is therefore a file with one line per chromosome. These chromosomal estimates are then averaged using `mawk` to give a global estimate which is used in *Fst* calculation.
