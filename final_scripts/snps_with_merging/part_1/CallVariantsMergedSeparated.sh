#!/bin/bash -i

mkdir separatedGVCFs
sampleListFile=$1
mapfile -t listOfSamples < $sampleListFile
for sample in "${listOfSamples[@]}"
do
    #option 1/2 depending on whether we mask the reference .fasta or not
    for primer in {1..30}
    do
        for allele in 1 2
        do
        echo "*******************************************************************"
        echo "Primer: $primer , sample: $sample , allele: $allele"
        echo "*******************************************************************"
        prefix="$sample""_Primer$primer""_allele$allele"
        inputBam="./separatedSams/$prefix.sorted.duplicates.bam"
        samtools index $inputBam
        numReads=$(samtools view -c $inputBam)
        
        echo "*******************************************************************"
        echo "Total number of reads: $numReads"
        echo "*******************************************************************"

        gatk --java-options "-Xmx4g" HaplotypeCaller \
        -R "./onePrimerFastas/Primer$primer.fasta" \
        -I  $inputBam \
        -O "./separatedGVCFs/$prefix.g.vcf.gz" \
        --disable-read-filter "NotDuplicateReadFilter" \
        --disable-read-filter "WellformedReadFilter" \
        -ERC GVCF
        done
    done
done