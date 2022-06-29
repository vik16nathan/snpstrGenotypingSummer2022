#!/bin/bash -i

#listOfSamples=("1232" "1393" "1791" "2006092" "2006174" "217" "2520669" "2599208" "34896" "34966")
#for sample in "${listOfSamples[@]}"
for sample in 1232
do
    for primer in {1..30}
    do
        for allele in 1 2
        do
        fastaPrefix="./maskedSeparatedSams/P-pyrhulla_$sample""_Primer$primer""_allele$allele""_masked_combined"
       
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
