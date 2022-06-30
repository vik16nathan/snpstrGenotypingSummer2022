#!/bin/bash -i

#figure out how to automate the column names somehow


sampleNames=($(sed -n '69,69p' P_pyrhulla_masked.vcf | awk '{for (i=10; i<NF-1; i++) printf "\""$i"\"" " "; printf "\""$i"\""}'))

fSampleNames=""
for sampleName in "${sampleNames[@]}"; do fSampleNames="$fSampleNames -F $sampleName "; done

#extract the primer/sample/allele column names from the .vcf file

gatk VariantsToTable \
    -V P_pyrhulla_masked.vcf \
    -F CHROM -F POS -F REF -F ALT -F FORMAT $fSampleNames \
    -O P_pyrhulla-1232-masked.table
