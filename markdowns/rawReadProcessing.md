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

Trimming was performed using `
