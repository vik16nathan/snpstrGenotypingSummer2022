#!/bin/bash -i

for sample in 1232
do
    for primer in {1..30}
    do
        for allele in 1 2
        do
        outputFilePrefix="P-pyrhulla_$sample""_Primer$primer""_allele$allele""_masked_combined_"
        fastaFile1="./maskedFasta/$outputFilePrefix""R1.fasta"
        fastaFile2="./maskedFasta/$outputFilePrefix""R2.fasta"

        fastaq to_fake_qual $fastaFile1 - | fastaq fasta_to_fastq $fastaFile1 - "./maskedFastq/$outputFilePrefix""R1.fastq"
        fastaq to_fake_qual $fastaFile2 - | fastaq fasta_to_fastq $fastaFile2 - "./maskedFastq/$outputFilePrefix""R2.fastq"

        done
    done
done
