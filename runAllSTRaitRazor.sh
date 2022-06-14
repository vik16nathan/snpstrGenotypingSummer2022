#!/bin/bash -i


listOfSamples=("1232" "1393" "1791" "2006092" "2006174" "217" "2520669" "2599208" "34896" "34966")
for sample in "${listOfSamples[@]}"
do
    for primer in {1..30}
    do
    primerString="Primer"
    primerArgString="${primerString}${primer}"
    filePrefix="P-pyrhulla_${sample}.sorted.duplicates_${primerArgString}"
    str8rzr -c "./config/Pyrhulla_pyrhulla_${primerArgString}_multiple_flanks.config" "./FastqInputs/${filePrefix}.fastq" > "./results/${filePrefix}_STRaitRazor_multiple_flanks.txt"

    #ensure that every line within the result file is unique (many lines will be duplicates due to multiple flanking
    #regions in the config file resulting in the same genotype due to a point mutation/SNP within a flanking region)
    #(head -n 1 ./results/${filePrefix}_STRaitRazor_multiple_flanks.txt && tail -n +2 ./results/${filePrefix}_STRaitRazor_multiple_flanks.txt | sort | uniq) > ./results/${filePrefix}_STRaitRazor_multiple_flanks.txt
    Rscript processMultipleFlanksSTRaitRazorResults.r $primer $sample
    done
done
