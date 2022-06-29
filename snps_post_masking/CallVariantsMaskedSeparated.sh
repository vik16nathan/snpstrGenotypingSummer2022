#!/bin/bash -i

#make dictionary for reference fasta
gatk CreateSequenceDictionary -R Pyrhulla_pyrhulla_PrimerRef.fasta

listOfSamples=("1232" "1393" "1791" "2006092" "2006174" "217" "2520669" "2599208" "34896" "34966")
#for sample in "${listOfSamples[@]}"
for sample in 1232
do
    for primer in {1..30}
    do
        for allele in 1 2
        do
        prefix="P-pyrhulla_$sample.sorted.duplicates_Primer$primer""_allele$allele""_masked_combined"
        samtools index "./maskedSeparatedSams/$prefix.sorted.duplicates.bam"
        gatk --java-options "-Xmx4g" HaplotypeCaller \
        -R Pyrhulla_pyrhulla_PrimerRef.fasta \
        -I  "./maskedSeparatedSams/$prefix.sorted.duplicates.bam" \
        -O "./maskedSeparatedGVCFs/$prefix.g.vcf.gz" \
        -ERC GVCF
        done
    done
done

