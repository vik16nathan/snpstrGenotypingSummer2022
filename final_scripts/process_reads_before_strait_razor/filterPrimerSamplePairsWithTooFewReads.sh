#!/bin/bash -i

sampleListFile=$1
threshold=$2
mapfile -t listOfSamples < $sampleListFile
for sample in "${listOfSamples[@]}"
do
    for primer in {1..30}
    do
        if [ $(samtools view -c -F 260 ./separatedBams/${sample}.sorted.duplicates_Primer${primer}.bam) -gt ${threshold} ]
        then
            cp ./separatedBams/${sample}.sorted.duplicates_Primer${primer}.bam ./separatedFilteredBams/
        fi
    done
done
