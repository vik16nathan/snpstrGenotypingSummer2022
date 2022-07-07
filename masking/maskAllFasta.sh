#!/bin/bash -i

listOfSamples=("1232" "1393" "1791" "2006092" "2006174" "217" "2520669" "2599208" "34896" "34966")
for sample in "${listOfSamples[@]}"
do
    for primer in {1..30}
    do
    
        for allele in 1 2
        do

        outputFilePrefix="P-pyrhulla_$sample""_Primer$primer""_allele$allele"
        #use masking tool on each DNA sequencing reads for each allele from each primer/sample pair
        bedtools maskfasta -fi "./FastaInputs/${outputFilePrefix}_R1.fasta" -bed "bedForMasking/P-pyrhulla_$sample""_Primer$primer""_allele_$allele.bed" -fo "./maskedFasta/$outputFilePrefix""_masked_R1.fasta"
        bedtools maskfasta -fi "./FastaInputs/${outputFilePrefix}_R2_reversed.fasta" -bed "bedForMasking/P-pyrhulla_$sample""_Primer$primer""_allele_$allele.bed" -fo "./maskedFasta/$outputFilePrefix""_masked_R2_reversed.fasta"

        #ensure that reads aren't split across multiple lines
        forwardFileName="./maskedFasta/$outputFilePrefix""_masked_combined_R1.fasta"
        r2reverseFileName="./maskedFasta/$outputFilePrefix""_masked_combined_R2_reversed.fasta"
        prefix="./maskedFasta/$outputFilePrefix""_masked_combined"
        (seqtk seq -a "./maskedFasta/$outputFilePrefix""_masked_R1.fasta") > $forwardFileName
        (seqtk seq -a "./maskedFasta/$outputFilePrefix""_masked_R2_reversed.fasta") > $r2reverseFileName

        #Reverse the R2 fasta lines to be back in their original 3'-5' orientation
        Rscript createReverseFasta.r $prefix $r2reverseFileName

        #create .fai file for new fastas
        samtools faidx "./maskedFasta/$outputFilePrefix""_masked_combined_R1.fasta"
        samtools faidx "./maskedFasta/$outputFilePrefix""_masked_combined_R2.fasta"

        #create .dict file
        gatk CreateSequenceDictionary -R "./maskedFasta/$outputFilePrefix""_masked_combined_R1.fasta"
        gatk CreateSequenceDictionary -R "./maskedFasta/$outputFilePrefix""_masked_combined_R2.fasta"

        done
    done
done

