#!/bin/bash -i

listOfSamples=("1232" "1393" "1791" "2006092" "2006174" "217" "2520669" "2599208" "34896" "34966")
for sample in "${listOfSamples[@]}"
do
    for primer in {1..30}
    do
    primerString="Primer"
    primerArgString="${primerString}${primer}"
    
    #First, use GATK to create a .bam file for one specific primer and one specific sample
    gatk PrintReads -I P-pyrhulla_${sample}.sorted.duplicates.bam -L ${primerArgString} -O \
    ./STRaitRazorGenotyping/separatedBams/P-pyrhulla_${sample}.sorted.duplicates_${primerArgString}.bam
    
    #Create a fastq file as a STRaitRazor input
    
    samtools fastq ./STRaitRazorGenotyping/separatedBams/P-pyrhulla_${sample}.sorted.duplicates_${primerArgString}.bam > \
    ./STRaitRazorGenotyping/FastqInputs/P-pyrhulla_${sample}.sorted.duplicates_${primerArgString}.fastq
    done
done
