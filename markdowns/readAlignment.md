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

Based on this, `bowtie2 --sensitive-local` was chosen as the best mapper for this data set. It showed the second highest percentage of reads mapped whilst being significantly faster than the setting with higher mapping percentage, `bowtie2 --very-sensitive-local`.
`BWA mem` showed fast performance but was only third-best in terms of mapping percentage.

