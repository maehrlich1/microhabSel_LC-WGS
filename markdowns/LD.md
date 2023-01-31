# Linkage Disequilibrium

LD calculation was performed in order to produce an LD pruned dataset for analyses where independence of data (SNPs) is required.
Specifically, LD was calculated between pairs of SNPs within a given proximity. Next, LD decay was plotted and the shape of the curve was used to inform LD pruning.

## LD Calculation
Pairwise LD statistics between SNPs were obtained using `ngsLD` with GLs as input:
```
N_SITES=$(wc -l $CHROM'.filt.pos')

ngsLD --geno $CHROM'.filt.beagle.gz' --probs --n_ind 938 --n_sites $N_SITES --pos $CHROM'.filt.pos' \
--max_kb_dist 10 --min_maf 0 \
--n_threads 10 \
| gzip > $CHROM'.filt.ld.gz'
```
A `--max_kb_dist 10` was chosen since previous *F. heteroclitus* dataset did not show significant linkage beyond 10Kb.
A `--min_maf 0` was chosen to speed up computation since sites had been MAF filtered previously already.

## LD Decay
Before plotting LD decay, LD output files were randomply downsampled to ~100 million pairwise comparisons in order to speed up calculations. Anything above 100,000 pairwise comparisons gives a decent distribution. To downsample the following `mawk` script was used:
```
zcat chr.filt.ld.gz | mawk '{if rand() <= 0.01) print $0}' | gzip > chr.filt.sample1p.ld.gz
```
Next LD decay was plotted using the `fit_LDdecay.R` script supplied with `ngsLD`:
```
Rscript --vanilla --slave ~/software/local/ngsLD/scripts/fit_LDdecay.py --ld_files input.txt --out chr.filt.sample1p.ld.plot \
--n_ind 956 --ld r2 --recomb_rate 2.34 --fit_boot 1000 --fit_bin_size 50 --fit_level 2
```
## LD Pruning
The LD decay curve informed reasonable cutoff values for LD pruning. More specifically, SNPs with r2 values above 0.1 were considered to be in linkage.
LD pruning was not performed on scaffolds since defining physical distances between scaffolds is not possible.
Strongest linked sites were identified using a `Python` script supplied with `ngsLD`:
```
python3 ~/software/local/ngsLD/scripts/prune_ngsLD.py \
--input $CHROM'.filt.ld.gz' --output $CHROM'.filt.lnkd' \
--max_dist 10000 --min_weight 0.1 --keep_heavy
```
 **NOTE:** The `--keep_heavy` option was used to output strongly linked sites rather than a pruned dataset. This is due to the `.ls.gz` files not containing all SNP information since only SNPs in proximity were considered for LD calculation. Distant SNPs do not feature in the `.ld.gz` files and will therefore not be retained although they are unlinked. For the same reason, chromosomes/scaffolds with only 1 SNP have empty LD files. These SNPs should however be included in the final SNP set.
 
 To generate the pruned, chromosomal SNP set, the linked sites were removed from the filtered dataset:
 ```
 zcat $CHROM'.filt.beagle.gz' | mawk -F '[:\t]' 'NR==FNR{array[$1"_"$2];next} !($1 in array) || FNR==1' $CHROM'.filt.lnkd' - | gzip > $CHROM'.filt.prune.beagle.gz'
 zcat $CHROM'.filt.mafs.gz' | mawk -F '[:\t]' 'NR==FNR{array[$1,$2];next} !(($1,$2) in array) || FNR==1' $CHROM'.filt.lnkd' - | gzip > $CHROM'.filt.prune.mafs.gz'
 ```
