#!/bin/bash -i


listOfSamples=("1232" "1393" "1791" "2006092" "2006174" "217" "2520669" "2599208" "34896" "34966")
#for sample in "${listOfSamples[@]}"
for sample in 1232
do
    for primer in {1..30}
    do
        for allele in 1 2
        do
        echo "*******************************************************************"
        echo "Primer: $primer , sample: $sample , allele: $allele"
        echo "*******************************************************************"
        prefix="P-pyrhulla_$sample""_Primer$primer""_allele$allele""_masked_combined"

        inputBam="./maskedSeparatedSams/$prefix.sorted.duplicates.bam"
        samtools index $inputBam
        numReads=$(samtools view -c $inputBam)
        
        echo "*******************************************************************"
        echo "Total number of reads: $numReads"
        echo "*******************************************************************"

        gatk --java-options "-Xmx4g" HaplotypeCaller \
        -R "./onePrimerFastas/Pyrhulla_pyrhulla_Primer$primer.fasta" \
        -I  $inputBam \
        -O "./maskedSeparatedGVCFs/$prefix.g.vcf.gz" \
        --disable-read-filter "NotDuplicateReadFilter" \
        --disable-read-filter "WellformedReadFilter" \
        -ERC GVCF
        done
    done
done