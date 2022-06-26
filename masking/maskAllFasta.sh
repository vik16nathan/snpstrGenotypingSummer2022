#!/bin/bash -i

listOfSamples=("1232" "1393" "1791" "2006092" "2006174" "217" "2520669" "2599208" "34896" "34966")
for sample in "${listOfSamples[@]}"
do
    for primer in {1..30}
    do
    
        for allele in 1 2
        do

        filePrefix="P-pyrhulla_$sample.sorted.duplicates_Primer$primer""_allele$allele"
        #use masking tool on each DNA sequencing reads for each allele from each primer/sample pair
        bedtools maskfasta -fi "./FastaInputs/$filePrefix.fasta" -bed "bedForMasking/P-pyrhulla_$sample""_Primer$primer""_allele_$allele.bed" -fo "./maskedFasta/$filePrefix""_masked.fasta"

        #create .fai file for new fastas
        samtools faidx "./maskedFasta/$filePrefix""_masked.fasta" -o "./maskedFasta/$filePrefix""_masked.fasta.fai"

        #create .dict file
        gatk CreateSequenceDictionary -R "./maskedFasta/$filePrefix""_masked.fasta"


        done
    done
done
