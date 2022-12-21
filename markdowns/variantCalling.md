# Variant Calling

## To Do:

-GL 1 (in most publications)

-doMajorMinor 1 (infers from data itself using GLs)

-doGlf 2 (output Beagle format GL file, most common)

-doMaf 1 (calculate MAF. use 1 just because Therkildsen uses it, could also do 2, a combination of both or even add 8)

-doCounts 1 (required for doDepth and dumpCounts)

-dumpCounts 1 (prints overall depth per locus - maybe use for fine-tuning parameters?!)

-doDepth 1 (gives depth distribution for every sample)

## Filters:

-setMinDepth (think about it, maybe use depth distribution?)

-setMaxDepth (think about it. depth distribution??)

-minInd (think about it)

-minQ (probably 20)

-minMapQ (also probably 20)

-SNP_pval 1e-6 (seems like what everybody uses)

-minMaf (think about it. Maybe not do minMaf but rather minimum number of alternative allele counts! i.e. a MAF of 0.01 is fine if it is supported by 10 individuals of 1000)

-baq 2 (I decided on 2)
