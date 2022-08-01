#!/bin/bash -i

sampleListFile=$1
mkdir dummyReads
mapfile -t listOfSamples < $sampleListFile
for sample in "${listOfSamples[@]}"
do
    for primer in {1..30}
    do
        for allele in 1 2
        do

            echo "Primer: $primer"
            echo "Sample: $sample"
            echo "Allele: $allele"
            prefix=${sample}_Primer${primer}_allele${allele}
            
            #Use fastp for merging
            fastp --merge -i ./FastqInputs/${prefix}_R1.fastq -I \
             ./FastqInputs/${prefix}_R2.fastq --merged_out \
             ./FastqInputs/${prefix}_merged.fastq --detect_adapter_for_pe --out1 ./dummyReads/${prefix}_out1.fastq --out2 \
             ./dummyReads/${prefix}_out2.fastq --unpaired1 ./dummyReads/${prefix}_unpaired1.fastq --unpaired2 \
             ./dummyReads/${prefix}_unpaired2.fastq
            
            echo ""
            echo "**************************************************************************************"
        done
    done
done

