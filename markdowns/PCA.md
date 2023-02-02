# Principal Component Analysis

Prior to running a PCA, we concatenate the `beagle.gz` files for all chromosomes into a single file:
```
zcat NC_0463*prune.beagle.gz | mawk 'NR>1 && $1=="marker" {next} {print $0}' | gzip > chr.filt.prune.beagle.gz
```
A PCA was performed using the entire filtered, pruned chromosomal SNP set.
`PCANGSD` was used to calculate the covariance matrix among samples using genotype likelihoods as input. A neighbour-joining tree is also produced from the covariance matrix:
```
pcangsd --beagle chr.filt.prune.beagle.gz --out chr.filt.prune \
--tree --tree_samples samples.list \
--threads 40
```
