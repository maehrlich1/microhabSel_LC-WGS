# Alignment Processing

In this step the aligned reads are:

* Given read group information for downstream analysis
* Merged between sequencing runs
* Deduplicated
* Overhang-clipped (Overlapping paired-end reads)

## Adding Read Group Info

By adding read group information individuals reads in a `BAM` file can be assigned to a specific library, sequencing instrument, lane and sample. This information is later used to correctly deduplicate reads and adjust base quality scores.
Read group information is added using the `Picard AddOrReplaceReadGroup` tool which comes with the `GATK4`. This requires information on:

* Library - a library is defined as a specific batch of DNA which was processed for sequencing. DNA batches can come from the same or different samples.
* Platform - Illumina in this case
* Platform Unit - This comprises the unique flowcell ID and lane number.
* Sample - The biological sample from which a library was made. There can be multiple libraries made from one sample.

Note: Since the `fastq` files supplied by the sequencing center were already demultiplexed and aggregated across sequencing lanes, lane information can no longer be assigned. Only the flowcell info was added in this case. While this does not allow for differentiating (and potentially accounting for) lane variance, it at least captures flowcell variance.

```
#Add Read Group to 1stStage
gatk --java-options -Xmx4g AddOrReplaceReadGroups -I '1stStage/'$SAMPLE_ID'_1stStage.sorted.bam' -O '1stStage/'$SAMPLE_ID'_1stStage.RG.sorted.bam' -RGLB $SAMPLE_ID'_lib1' -RGPL illumina -RGPU GW2011023488th -RGSM $SAMPLE_ID && rm '1stStage/'$SAMPLE_ID'_1stStage.sorted.bam'*

#Add Read Group to 2ndStage
gatk --java-options -Xmx4g AddOrReplaceReadGroups -I '2ndStage/'$SAMPLE_ID'_2ndStage.sorted.bam' -O '2ndStage/'$SAMPLE_ID'_2ndStage.RG.sorted.bam' -RGLB $SAMPLE_ID'_lib1' -RGPL illumina -RGPU GW201223000 -RGSM $SAMPLE_ID && rm '2ndStage/'$SAMPLE_ID'_2ndStage.sorted.bam'*
```

## Merging Reads

Once read group info has been assigned, the `BAM` files resulting from sequencing the same library/sample on several lanes/flowcells, can be merged. For this we use `Picard MergeSamFiles`:
```
gatk --java-options -Xmx4g MergeSamFiles -I '1stStage/'$SAMPLE_ID'_1stStage.RG.sorted.bam' -I '2ndStage/'$SAMPLE_ID'_2ndStage.RG.sorted.bam' -O 'merged/'$SAMPLE_ID'.merged.sorted.bam' --USE_THREADING true
```

## Marking Duplicate reads

Duplicate reads are removed using `Picard MarkDuplicates`. Duplication stats are printed to a separate file.

```
java -Xmx4g -jar $PICARD_DIR/MarkDuplicates.jar I=$SAMPLE_ID'.merged.sorted.bam' O=$SAMPLE_ID'.merged.sorted.dedup.bam' M=$SAMPLE_ID'_dupstat.txt' REMOVE_DUPLICATES=true
```

## Clip Overlapping Read Pairs

Overlapping paired-end reads are clipped such that sequence information in the overlapping section is not duplicated. Read information for the higher quality read (usually the forward read) is retained and the secondary read clipped. For this, the `clipoverlap` function in `Bamutil` was used.

```
bam clipOverlap --in $SAMPLE_ID'.merged.sorted.dedup.bam' --out $SAMPLE_ID'.merged.sorted.dedup.overlapclip.bam' --stats
```
