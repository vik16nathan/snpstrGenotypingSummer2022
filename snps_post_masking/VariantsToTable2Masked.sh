#!/bin/bash -i

sampleNames=($(sed -n '69,69p' P_pyrhulla_masked.vcf | awk '{for (i=10; i<NF-1; i++) printf "\""$i"\"" " "; printf "\""$i"\""}'))
fSampleNames=""
for sampleName in "${sampleNames[@]}"; do fSampleNames="$fSampleNames -F $sampleName "; done

gatk VariantsToTable \
-V P_pyrhulla_masked_annotated.vcf \
 -F CHROM -F POS -F REF -F ALT -F FORMAT $fSampleNames \
-O P_pyrhulla-SNPs_masked_1232.table

