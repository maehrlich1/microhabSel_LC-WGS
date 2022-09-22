# Raw-read Processing

This markdown details:

* Quality Control
* Read Trimming
* Post-trim QC

## Quality Control

Quality control was performed using `FastQC` and samples were jointly analyzed using `MultiQC` using only the `-m fastqc` option for faster performance:

```
fastqc -t 39 -o /projects2/rsmas/dcrawford/MAE/microhabSel_LC-WGS/fastqc/fastqc_2ndStage *fastq.gz

multiqc -m fastqc .
```

## Read Trimming

Read trimming was performed using `Trimmomatic v.0.39` with the following settings:

`ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:1:true`

* `NexteraPE-PE.fa` supplied with `Trimmomatic v.0.39` was used for adapter clipping.
* Seed Mismatches = 2
* Palindrome clip threshold = 30
* Simple clip threshold = 10
* Min Adapter length = 1
* Keep reverse read upon read-through detection = TRUE

`SLIDINGWINDOW:4:15`

* Clip reads once the average base quality in a 4-base window drops below 15

`LEADING:3`

* Clip leading bases with quality <3

`TRAILING:3`

* Clip trailing bases with quality <3


```
#SETUP
#change directory
cd /projects2/rsmas/dcrawford/MAE/microhabSel_LC-WGS/raw_reads/1stStage/

#SET VARIABLES
SAMPLETABLE=/projects2/rsmas/dcrawford/MAE/microhabSel_LC-WGS/sample_table.tsv
TRIMMOMATIC=~/software/local/Trimmomatic-0.39/trimmomatic-0.39.jar
THREADS=4  #Optimal number of threads for Trimmomatic since output compression is limiting and more threads lead to no performance increase.
ADAPTERS=~/software/local/Trimmomatic-0.39/adapters/NexteraPE-PE.fa
R1_SUFFIX=_R1_001.fastq.gz
OUTDIR=/projects2/rsmas/dcrawford/MAE/microhabSel_LC-WGS/trimmed/1stStage/
OUT_SUFFIX=_trimmed.fastq.gz

#CODE
#Parallelize Trimmomatic. Since Trimmo is optimized at 4 threads, we can take advantage of the up to 40 cores on the cluster. However, since Trimmo uses 4 threads already we reduce the jobload per core to 25% (=1/4 threads).
parallel -j 25% --header : --colsep '\t' java -jar $TRIMMOMATIC PE -threads $THREADS -basein {fastq_prefix}$R1_SUFFIX -baseout $OUTDIR{sample}$OUT_SUFFIX ILLUMINACLIP:$ADAPTERS:2:30:10:1:true SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 :::: $SAMPLETABLE
```
Only 4 threads were used for `Trimmomatic` since output compression is the bottleneck. Further threads would not improve performance.

However, in order to take advantage of the HPC cluster, the script was parallelized using GNU `parallel` and the automatic behavior of 1 job per thread suppressed using `-j 25%` since `Trimmomatic` is already using 4 threads.

## Quality Control - Round 2

A second round of QC is performed to check whether adapter trimming was performed as expected.

The code used is essentially the same as above:

```
fastqc -t 39 -o /projects2/rsmas/dcrawford/MAE/microhabSel_LC-WGS/trimmed/fastqc/fastqc_1stStage_trimmed *fastq.gz

multiqc -m fastqc .
```
