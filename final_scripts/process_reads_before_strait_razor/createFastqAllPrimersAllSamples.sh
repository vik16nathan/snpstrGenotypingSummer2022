#!/bin/bash -i
sampleListFile=$1
mapfile -t listOfSamples < $sampleListFile

for sample in "${listOfSamples[@]}"
do
    for primer in {1..30}
    do
    primerString="Primer"
    primerArgString="${primerString}${primer}"
    
    #Create a fastq file as a STRaitRazor input
    #Only create STRaitRazor input files for .bam files in the separatedFilteredBams directory
    #We only want to run STRaitRazor on inputs for which there are >10 reads
    
    samtools fastq ./separatedFilteredBams/${sample}.sorted.duplicates_${primerArgString}.bam > \
    ./FastqInputs/${sample}.sorted.duplicates_${primerArgString}.fastq
    done
done