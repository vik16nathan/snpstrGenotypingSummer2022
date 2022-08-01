#!/bin/bash

sampleListFile=$1
mapfile -t listOfSamples < $sampleListFile
for sample in "${listOfSamples[@]}"
do
    for primer in {1..30}
    do
        echo "Primer: $primer"
        echo "Sample: $sample"
        echo ""
        Rscript findSTRpositionsFastqEachAllele.r $primer $sample
        echo "********************************************************************************"

    done
done