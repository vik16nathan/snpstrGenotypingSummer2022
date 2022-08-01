#!/bin/bash -i

#Input: Primer reference fasta file with 30 lines, each of which corresponds to a primer
#Output: 30 separate primer reference fasta files, each of which has a single header and a single line

mkdir onePrimerFastas
for i in {1..30}
do
    sed -n "$((2*i - 1)),$((2*i))p" PrimerRef.fasta > ./onePrimerFastas/Primer$i.fasta
    bwa index ./onePrimerFastas/Primer$i.fasta
    samtools faidx ./onePrimerFastas/Primer$i.fasta
    
    #make dictionary for reference fasta
    gatk CreateSequenceDictionary -R ./onePrimerFastas/Primer$i.fasta
done

