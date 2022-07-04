#!/bin/bash -i

#Input: Primer reference fasta file with 30 lines, each of which corresponds to a primer
#Output: 30 separate primer reference fasta files, each of which has a single header and a single line

for i in {1..30}
do
    sed -n "$((2*i - 1)),$((2*i))p" Pyrhulla_pyrhulla_PrimerRef.fasta > ./onePrimerFastas/Pyrhulla_pyrhulla_Primer$i.fasta
    bwa index ./onePrimerFastas/Pyrhulla_pyrhulla_Primer$i.fasta
    samtools faidx ./onePrimerFastas/Pyrhulla_pyrhulla_Primer$i.fasta
    
    #make dictionary for reference fasta
    gatk CreateSequenceDictionary -R ./onePrimerFastas/Pyrhulla_pyrhulla_Primer$i.fasta
done

