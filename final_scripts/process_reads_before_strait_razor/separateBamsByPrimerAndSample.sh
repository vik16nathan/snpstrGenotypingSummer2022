#!/bin/bash -i

#Use GATK to create a .bam file for one specific primer and one specific sample
sampleListFile=$1
mapfile -t listOfSamples < $sampleListFile
for sample in "${listOfSamples[@]}"
do
    samtools index ${sample}.sorted.duplicates.bam
    for primer in {1..30}
    do
        gatk PrintReads -I ${sample}.sorted.duplicates.bam -L Primer${primer} -O \
            ./separatedBams/${sample}.sorted.duplicates_Primer${primer}.bam
    done
done