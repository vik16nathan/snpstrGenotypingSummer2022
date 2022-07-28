#!/bin/bash -i

sampleListFile=$1
mapfile -t listOfSamples < $sampleListFile

for sample in "${listOfSamples[@]}"
do
    for primer in {1..30}
    do
    
    primerString="Primer"
    primerArgString="${primerString}${primer}"
    filePrefix="${sample}.sorted.duplicates_${primerArgString}"

    #the f -10 argument ensures that the output files keep only the alleles that were found >10 times by STRait Razor 
    str8rzr -f 10 -c "./config/${primerArgString}_multiple_flanks.config" "./FastqInputs/${filePrefix}.fastq" > "./results/${filePrefix}_STRaitRazor_multiple_flanks.txt"

    #ensure that every line within the result file is unique (many lines will be duplicates due to multiple flanking
    #regions in the config file resulting in the same genotype due to a point mutation/SNP within a flanking region)

    #We need a script to combine all reads corresponding to the same alleles, regardless of which flanking regions in the .config file were used to generate them
    #The output needs to be in reverse sorted order by total number of forward and reverse reads, and should contain only unique alleles
    Rscript processMultipleFlanksSTRaitRazorResults.r $primer $sample
    done
done
