#!/bin/bash -i
sampleListFile=$1
mkdir dummyReads
mkdir separatedSams
mapfile -t listOfSamples < $sampleListFile
for sample in "${listOfSamples[@]}"
do
    for primer in {1..30}
    do
        for allele in 1 2
        do

            fastqPrefix="$sample""_Primer$primer""_allele$allele"
            bwa index "./FastqInputs/$fastqPrefix""_merged.fastq"
            SampleID="$primer""_$sample""_$allele"
            bwa mem -R "@RG\tID:${SampleID}\tLB:$sample\tPL:ILLUMINA\tPM:HISEQ\tSM:${SampleID}" \
            "./onePrimerFastas/Primer$primer.fasta" "./FastqInputs/$fastqPrefix""_merged.fastq" > \
            "./separatedSams/$fastqPrefix.mapped.sam"

        done
    done
done

