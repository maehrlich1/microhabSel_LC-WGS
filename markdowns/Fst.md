# *Fst* Calculation
*Fst* calculation using `ANGSD` requires the site allele frequencies which can then be used to produce site frequency spectra which in turn can be used to calculate *Fst*.
The procedure as a whole is:
* Use `ANGSD` for calculating `.saf` files for each population separately
* Use the `realSFS` script found in `ANGSD/misc/realSFS` to calculate 2D SFS for each population pair
* Use the above calculated 2D-SFS as priors jointly with all `.saf` files from step 1 to calculate *Fst* binary files
* Use `realSFS` to extract the *Fst* values from the *Fst* file

## SAF Calculation
Since the SFS requires the relative proportion of monomorphic sites too, SAF calculation requires information from the entire genome and not just SNPs.
Therefore raw `BAM` files are once again provided as input:
```
angsd -bam '../'$POP'.bamlist' -ref $REF -anc $REF -out $POP \
-GL 1 -doSaf 1 \
-minMapQ 20 -baq 2 -minQ 20 \
-P 8 #ANGSD only takes a maximum of 8 threads. Due to I/O operations being the bottleneck.
```
