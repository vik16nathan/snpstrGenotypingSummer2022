#!/bin/bash -i

listOfSamples=("1232" "1393" "1791" "2006092" "2006174" "217" "2520669" "2599208" "34896" "34966")
for sample in "${listOfSamples[@]}"
do
    for primer in {1..30}
    do
    primerString="Primer"
    primerArgString="${primerString}${primer}"
    
    #Create a fastq file as a STRaitRazor input
    #Only create STRaitRazor input files for .bam files in the separatedFilteredBams directory
    #We only want to run STRaitRazor on inputs for which there are >10 reads
    
    samtools fastq ./STRaitRazorGenotyping/separatedFilteredBams/P-pyrhulla_${sample}.sorted.duplicates_${primerArgString}.bam > \
    ./STRaitRazorGenotyping/FastqInputs/P-pyrhulla_${sample}.sorted.duplicates_${primerArgString}.fastq
    done
done
