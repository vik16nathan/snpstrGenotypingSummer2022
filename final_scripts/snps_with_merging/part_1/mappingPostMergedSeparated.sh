#!/bin/bash -i
sampleListFile=$1
mapfile -t listOfSamples < $sampleListFile
for sample in "${listOfSamples[@]}"
do
    for primer in {1..30}
    do
        for allele in 1 2
        do
        fastaPrefix="./separatedSams/$sample""_Primer$primer""_allele$allele"

        # clean up read pairing information and flag with SAMtools:
        samtools sort -n -O sam "$fastaPrefix.mapped.sam" | \
        samtools fixmate -m -O bam - "$fastaPrefix.mapped.fixmate.bam"

        # sort the bam-file into coordinate order:
        samtools sort -O bam -o "$fastaPrefix.sorted.bam" "$fastaPrefix.mapped.fixmate.bam"

        # mark duplicated
        samtools markdup -S "$fastaPrefix.sorted.bam" "$fastaPrefix.sorted.duplicates.bam"

        done
    done
done
