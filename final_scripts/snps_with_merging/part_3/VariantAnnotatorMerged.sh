#!/bin/bash -i

gunzip MergedVariants_filteredG.vcf.gz
inputIString=""
for file in $(ls ./separatedSams/ | grep sorted.duplicates.bam$ )
do
    inputIString="$inputIString  ""-I ./separatedSams/$file"
done

gatk  VariantAnnotator \
   -R PrimerRef.fasta \
   $inputIString \
   -V MergedVariants_filteredG.vcf \
   -O MergedVariantsFilteredAnnotated.vcf \
   -A Coverage

