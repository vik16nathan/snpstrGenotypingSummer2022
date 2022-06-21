#!/bin/bash -i


listOfSamples=("1232" "1393" "1791" "2006092" "2006174" "217" "2520669" "2599208" "34896" "34966")
for sample in "${listOfSamples[@]}"
do
    for primer in {1..30}
    do
        #for flankLength in {10}
        #do
        flankLength=10
        primerString="Primer"
        primerArgString="${primerString}${primer}"
        filePrefix="P-pyrhulla_${sample}.sorted.duplicates_${primerArgString}"
        (str8rzr -f 10 -c "./config/Pyrhulla_pyrhulla_${primerArgString}_flank${flankLength}.config" "./FastqInputs/${filePrefix}.fastq" | sed '/SumBelowThreshold/d') > "./results/${filePrefix}_STRaitRazor_flank${flankLength}.txt"
        #done
    done
done

