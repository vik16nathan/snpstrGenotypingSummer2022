#!/bin/bash -i

Ref="PrimerRef.fasta"
samtools faidx $Ref
gatk CreateSequenceDictionary -R $Ref

 gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R $Ref \
   -V gendb://merged-separated \
   -O MergedVariantsBeforeFiltering.vcf.gz \

