# Principal Component Analysis

A PCA was performed using the filtered, pruned chromosomal SNP set.
`PCANGSD` was used to calculate the covariance matrix among samples using genotype likelihoods as input. A neighbour-joining tree is also produced from the covariance matrix:
```
pcangsd --beagle chr.filt.prune.beagle.gz --out chr.filt.prune \
--tree --tree_samples samples.list \
--threads 40
```
