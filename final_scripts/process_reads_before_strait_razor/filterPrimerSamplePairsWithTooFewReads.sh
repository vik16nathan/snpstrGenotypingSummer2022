#!/bin/bash -i

sampleListFile=$1
mapfile -t listOfSamples < $sampleListFile
for sample in "${listOfSamples[@]}"
do
    for primer in {1..30}
    do
        if [ $(samtools view -c -F 260 ./separatedBams/${sample}.sorted.duplicates_Primer${primer}.bam) -gt 10 ]
        then
            cp ./separatedBams/${sample}.sorted.duplicates_Primer${primer}.bam ./separatedFilteredBams/
        fi
    done
done