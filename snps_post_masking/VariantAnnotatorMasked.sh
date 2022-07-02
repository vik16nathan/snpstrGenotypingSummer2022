#!/bin/bash -i

gunzip P_pyrhulla_masked_filteredG.vcf.gz
inputIString=""
for file in $(ls ./maskedSeparatedSams/ | grep sorted.duplicates.bam$ )
do
    inputIString="$inputIString  ""-I ./maskedSeparatedSams/$file"
done
gatk  VariantAnnotator \
   -R Pyrhulla_pyrhulla_PrimerRef.fasta \
   $inputIString \
   -V P_pyrhulla_masked_filteredG.vcf \
   -O P_pyrhulla_masked_annotated.vcf \
   -A Coverage
