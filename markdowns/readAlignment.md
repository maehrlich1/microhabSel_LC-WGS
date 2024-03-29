# Read Alignment

## Aligner Software Benchmarking

`BWA mem` and `Bowtie2` aligners were compared using `Teaser` benchmarking software.
A subset of the 957 samples was provided to `Teaser` which further randomly subsampled reads from each sample and aligned them to the *Fundulus heteroclitus* reference genome.

`BWA mem` was tested in default mode whilst `bowtie2` was tested with various parameter settings.

`Teaser` takes a `YAML` setup file as input:

```
include:
  - base_teaser.yaml

test_mappers:
  - bowtie2
  - bwamem

test_parameters:
  - bowtie2_presets    #Tells Teaser to run bowtie2 with various settings

teaser:
   tests:
      microhabSel_LC-WGS_1st_16SP1X_000046:
         type: real  #Not simulated data
         reference: GCF_011125445.2_MU-UCD_Fhet_4.1_genomic.fna
         paired: Yes
         import_read_files: [/projects2/rsmas/dcrawford/MAE/microhabSel_LC-WGS/trimmed/teaserTest/1st/16SP1X_000046_trimmed_1P.fastq,/projects2/rsmas/dcrawford/MAE/microhabSel_LC-WGS/trimmed/teaserTest/1st/16SP1X_000046_trimmed_2P.fastq]

threads: 5

report:
    name: microhabSel_LC-WGS_1st_16SP1X_000046
```

By providing the `YAML` file, `Teaser` is simply run by:

```
./teaser.py example.yaml
```

`Teaser` generates a `HTML` report containing information on:
* Percentage mapped reads
* Percentage unmapped reads
* Runtime
* Memory Usage
* Mapping quality

Based on this, `bwa mem` was chosen as the best mapper for this data set. It showed the best response curve in a plot of read count as a function of mapping quality. Although it did not map the most reads, it mapped most reads with a mapping quality of 20 and above.
[Rplot.pdf](https://github.com/maehrlich1/microhabSel_LC-WGS/files/9718908/Rplot.pdf)

## Read Alignment

Paired reads were aligned using `bwa mem` as follows:
```
bwa mem -M -a reference.fa $SAMPLE_ID'_trimmed_1P.fastq.gz' $SAMPLE_ID'_trimmed_2P.fastq.gz' > $SAMPLE_ID'.sam'
```
The `-M -a` options:
* Mark shorter split hits as secondary (For Picard compatibility)
* Outputs all found alignments for unpaired reads also. These are marked as secondary.

Before further processing `SAM` files were converted to `BAM` format and sorted using `sambamba 0.8.2`:
```
sambamba view -t 4 -S -f bam mysample.sam | sambamba sort -t 4 -m 8G -o mysample.sorted.bam' /dev/stdin && rm mysample.sam
```
## Alignment Quality Control

`Qualimap` was used to assess alignment statistics across all samples, grouped by sequencing batch. Java memory size was expanded beyond the size of the largest `BAM` file for faster processing.
```
qualimap multi-bamqc --java-mem-size=4G -r -d qualimap_list.tsv -outdir ./
```
