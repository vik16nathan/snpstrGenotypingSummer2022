#!/bin/bash -i

sampleNames=($(sed -n '1,1p' MergedVariantsFilteredAnnotatedTableOnly.vcf | awk '{for (i=10; i<NF-1; i++) printf "\""$i"\"" " "; printf "\""$i"\""}'))
fSampleNames=""
for sampleName in "${sampleNames[@]}"; do fSampleNames="$fSampleNames -F $sampleName "; done

gatk VariantsToTable \
-V MergedVariantsFilteredAnnotated.vcf \
-F CHROM -F POS -F REF -F ALT -F FORMAT $fSampleNames \
-O MergedVariantsFilteredAnnotated.table

#Prepare directory for final tables
mkdir finalExcelOutputs