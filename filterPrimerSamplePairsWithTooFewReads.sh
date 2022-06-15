#!/bin/bash -i

listOfSamples=("1232" "1393" "1791" "2006092" "2006174" "217" "2520669" "2599208" "34896" "34966")
for sample in "${listOfSamples[@]}"
do
    for primer in {1..30}
    do
        if $(samtools view -c -F 260 ./separatedBams/P-pyrhulla_${sample}_sorted.duplicates_Primer${primer}.bam) < 10
        then
            echo "Primer $primer and sample $sample has < 10 reads"
        fi
    done
done


