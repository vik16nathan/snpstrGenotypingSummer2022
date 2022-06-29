#!/bin/bash -i

#rename to mapMaskedFastq.sh

#listOfSamples=("1232" "1393" "1791" "2006092" "2006174" "217" "2520669" "2599208" "34896" "34966")
for sample in 1232
do
    for primer in {1..30}
    do
        for allele in 1 2
        do

        fastqPrefix="P-pyrhulla_$sample""_Primer$primer""_allele$allele""_masked_combined"
        bwa index "./maskedFastq/$fastqPrefix""_R1.fastq"
        bwa index "./maskedFastq/$fastqPrefix""_R2.fastq"
        SampleID="$sample""_$allele"
        bwa mem -R "@RG\tID:${SampleID}\tLB:$sample\tPL:ILLUMINA\tPM:HISEQ\tSM:${SampleID}" ./onePrimerFastas/Pyrhulla_pyrhulla_Primer$primer.fasta "./maskedFastq/$fastqPrefix""_R1.fastq" "./maskedFastq/$fastqPrefix""_R2.fastq" > "./maskedSeparatedSams/$fastqPrefix.mapped.sam"

        done
    done
done

