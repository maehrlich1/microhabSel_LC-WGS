# Principal Component Analysis

Prior to running a PCA, the filtered and pruned SNPset was further reduced by only including SNPs with a MAF > 0.1. SNP filtering according to missingness would have been possible also but the distribution of "Individual Depth" was too narrow to allow for further filtering i.e. a large number of SNPs would have been discarded.

To filter the SNP set the following code was used:
```
#Filter the .mafs.gz file for SNPs with MAF >= 0.1
zcat chr.filt.prune.mafs.gz | mawk 'FNR==1 || $6 >= 0.1' | gzip > chr.filt.prune.maf10p.mafs.gz

#Use the new MAF file to filter the .beagle.gz file
zcat chr.filt.prune.beagle.gz | mawk 'NR==FNR && NR>1{array[$1"_"$2];next} $1 in array || FNR==1 && NR!=FNR' <(zcat chr.filt.prune.maf10p.mafs.gz) - | gzip > chr.filt.prune.maf10p.beagle.gz
```
A PCA was performed using the filtered, pruned chromosomal SNP set with MAF >0.1.
`PCANGSD` was used to calculate the covariance matrix among samples using genotype likelihoods as input. A neighbour-joining tree is also produced from the covariance matrix:
```
pcangsd --beagle chr.filt.prune.maf10p.beagle.gz --out chr.filt.prune.maf10p \
--tree --tree_samples samples.list \
--threads 40
```
