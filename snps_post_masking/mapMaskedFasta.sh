#!/bin/bash -i

#listOfSamples=("1232" "1393" "1791" "2006092" "2006174" "217" "2520669" "2599208" "34896" "34966")
for sample in 1232
do
    #for primer in {1..30}
    for primer in {1..15}
    do
        for allele in 1 2
        do

        fastaPrefix="P-pyrhulla_$sample.sorted.duplicates_Primer$primer""_allele$allele""_masked_combined"
        bwa index "./maskedFasta/$fastaPrefix.fasta"
        SampleID="$primer""_$sample""_$allele"
        bwa mem -R "@RG\tID:${SampleID}\tLB:$sample\tPL:ILLUMINA\tPM:HISEQ\tSM:${SampleID}" Pyrhulla_pyrhulla_PrimerRef.fasta "./maskedFasta/$fastaPrefix""_R1.fasta" "./maskedFasta/$fastaPrefix""_R2.fasta" > "./maskedSeparatedSams/$fastaPrefix.mapped.sam"

        done
    done
done

